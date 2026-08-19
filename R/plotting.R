#' Compute per-assay mean and variance by sample type
#'
#' @param data A data frame with a `SampleType` column and one column per
#'   assay/feature (columns 2 and 3 are assumed to be non-feature columns
#'   and are dropped, matching the original assay layout).
#'
#' @return A long data frame with columns `SampleType`, `Assay`, `mean`,
#'   `variance`.
#' @export
graph_mean_var <- function(data) {
  # calculate the means of each assay for all samples and controls
  x_mean <-
    data |>
    dplyr::group_by(SampleType) |>
    dplyr::summarise_all(~ mean(na.omit(.))) |>
    dplyr::select(-c(2, 3)) |>
    tidyr::pivot_longer(!SampleType, names_to = "Assay", values_to = "mean")

  # calculate variance for each assay for all the samples and controls
  x_var <-
    data |>
    dplyr::group_by(SampleType) |>
    dplyr::summarise_all(~ stats::var(na.omit(.))) |>
    dplyr::select(-c(2, 3)) |>
    tidyr::pivot_longer(!SampleType, names_to = "Assay", values_to = "variance")

  # combine the two data frames together
  data_graph <-
    x_mean |>
    dplyr::left_join(x_var, by = c("SampleType", "Assay"))
  return(data_graph)
}

#' Bubble-style correlation plot with significance sizing
#'
#' Plots a correlation (or similar) matrix as a grid of points where fill
#' encodes the value and point size encodes `-log10(p-value)`.
#'
#' @param corr A numeric matrix of values to plot (e.g. correlations).
#' @param p.mat A numeric matrix of p-values, same dimensions as `corr`.
#' @param labels Optional character vector of x-axis labels; if `NULL`,
#'   labels are derived from `colnames(corr)` by stripping a trailing
#'   `stat.mean` suffix.
#'
#' @return A `ggplot` object.
#' @export
customCorrplot <- function(corr, p.mat, labels = NULL) {
  data <- corr
  dt <- data.frame()
  for (c in c(1:ncol(data))) {
    i <- data.frame(x = rep(c, nrow(data)), y = c(1:nrow(data)), var = data[, c])
    rownames(i) <- rownames(data)
    dt <- rbind(dt, i)
  }
  data.p <- p.mat
  dt.p <- data.frame()
  for (c in c(1:ncol(data.p))) {
    data.p[, c] <- as.numeric(as.character(data.p[, c]))
    i <- data.frame(x = rep(c, nrow(data.p)), y = c(1:nrow(data.p)), pval = -log10(data.p[, c]))
    rownames(i) <- rownames(data.p)
    dt.p <- rbind(dt.p, i)
  }
  dt.merge <- dplyr::left_join(dt, dt.p, by = c("x", "y"))
  dt.merge$var <- as.numeric(as.character(dt.merge$var))
  dt.merge$pval <- as.numeric(as.character(dt.merge$pval)) # numeric value
  if (!is.null(labels)) {
    labels <- labels
  } else {
    labels <- stringr::str_extract(string = colnames(data), pattern = ".*?(?=stat.mean)") # get rid of stat.mean labels
    labels <- stringr::str_replace_all(string = labels, pattern = "\\.", replacement = " ") # customized labels
  } # push in customized labels if provided

  p <- ggplot2::ggplot(
    data = dt.merge,
    ggplot2::aes(x = x, y = y, size = pval, fill = var)
  ) +
    ggplot2::geom_point(shape = 21, stroke = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_size_continuous(range = c(3, 10)) +
    ggplot2::scale_fill_gradient(low = "blue", high = "red") +
    ggplot2::scale_x_continuous(
      limits = c(0, (ncol(data) + 1)),
      breaks = c(0:(ncol(data) + 1)),
      labels = c("", labels, "")
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(1:nrow(data)),
      labels = rownames(data)
    ) +
    ggplot2::theme(legend.position = "left") +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = 12,
        hjust = 1,
        vjust = 0.5,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 12,
        hjust = 1,
        vjust = 0.5,
        face = "bold"
      )
    ) +
    ggplot2::guides()
  return(p)
}

#' Plot Olink-style QC summary (pie chart + histogram) per plate
#'
#' @param data_i A list-like object with elements `AssayOlinkQC` (columns
#'   `AssayQC`, `n`, `PlateID`) and `SampleOlinkQC` (columns `Frequency`,
#'   `PlateID`).
#'
#' @return The combined grid graphics object returned by
#'   [gridExtra::grid.arrange()].
#' @export
Plot_QC <- function(data_i) {
  plot <-
    data_i$AssayOlinkQC |>
    dplyr::arrange(n) |>
    dplyr::mutate(AssayQC = factor(AssayQC, levels = c("PASS", "WARN", "FAIL"))) |>
    dplyr::mutate(prop = n / sum(n) * 100) |>
    dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop) |>
    ggplot2::ggplot(ggplot2::aes(x = "", y = n, fill = AssayQC)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggrepel::geom_label_repel(
      ggplot2::aes(y = ypos, label = paste0(AssayQC, "=", n)),
      size = 4, nudge_x = 0.8, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = c("darkgreen", "yellow", "red")) +
    ggplot2::facet_wrap(~PlateID)

  plot_hist <-
    data_i$SampleOlinkQC |>
    ggplot2::ggplot(ggplot2::aes(x = 1 - Frequency, fill = PlateID)) +
    ggplot2::geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity") +
    ggplot2::stat_count(
      binwidth = 1,
      geom = "text",
      color = "black",
      ggplot2::aes(label = ggplot2::after_stat(count)), size = 4,
      position = ggplot2::position_stack(vjust = 1.02)
    ) +
    ggplot2::xlab("% analytes failed") +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~PlateID)

  return(gridExtra::grid.arrange(grobs = list(plot, plot_hist), nrow = 2))
}


#' Proportional Venn / Euler diagram
#'
#' @param list A named list of vectors, one per set, as expected by
#'   [gplots::venn()].
#' @param cols A named vector of colors, keyed by set-combination name (as
#'   produced by [gplots::venn()], with `":"` replaced by `"&"`).
#'
#' @return A list with `plot` (the plotted [eulerr::euler()] object) and
#'   `table` (named vector of intersection sizes).
#' @export
euler_venn <- function(list, cols) {
  set.seed(123)
  items <- gplots::venn(list, show.plot = FALSE)
  combo <- lengths(attributes(items)$intersections)
  names(combo) <- gsub(":", "&", names(combo))
  cols <- cols[names(combo)[1:length(list)]]
  plot_venn <- eulerr::euler(combo, shape = "circle")
  plot <- plot(
    plot_venn,
    labels = list(cex = 2),
    edges = list(lwd = 2, col = "black"),
    fills = list(alpha = 0.9, fill = cols),
    quantities = list(
      type = c("counts", "percent"),
      cex = 1.5,
      fontface = "bold",
      col = "black"
    )
  )
  return(list(plot = plot, table = combo))
}


#' Volcano plot with pathway-colored gene highlighting
#'
#' @param data A data frame containing at least `x`, `y`, and `id_col`.
#' @param x Character scalar, column name to use on the x-axis (typically
#'   log fold-change).
#' @param y Character scalar, column name to use on the y-axis (a p-value,
#'   or already-transformed `-log10(p)` if `raw_pvalue = FALSE`).
#' @param path_list A data frame with one row per pathway, containing
#'   `path_name`, `"Description"`, and `path_id` (comma-separated gene IDs).
#' @param path_name Column name in `path_list` giving the pathway
#'   identifier/label used for coloring. Defaults to `"Pathway.ID"`.
#' @param legend_position Passed through conceptually; retained for
#'   interface parity (legend position is set via `ggplot2::theme()`
#'   internally).
#' @param path_id Column name in `path_list` with comma-separated gene IDs.
#'   Defaults to `"Genes"`.
#' @param raw_pvalue Logical; if `TRUE`, `y` is `-log10()`-transformed.
#' @param point_size Base point size.
#' @param text_size Label text size.
#' @param line_size Cutoff line width.
#' @param x_cutoff Length-2 numeric, vertical dashed-line positions.
#' @param y_cutoff Numeric, horizontal dashed-line position.
#' @param id_col Column name in `data` with the gene/feature identifier.
#' @param xlab,ylab Axis labels.
#' @param cols Named vector of colors keyed by pathway label.
#' @param max.overlaps Passed to [ggrepel::geom_text_repel()].
#'
#' @return A `ggplot` object.
#' @export
AdvancedVolcano <- function(data, x, y, path_list, path_name = "Pathway.ID",
                             legend_position = "top", path_id = "Genes",
                             raw_pvalue = TRUE, point_size = 2.5, text_size = 2.5,
                             line_size = 0.6, x_cutoff = c(-0.2, 0.2),
                             y_cutoff = -log10(0.01), id_col = "GeneID",
                             xlab = "Log2FC", ylab = "-log10(P-Value)",
                             cols = grDevices::palette(), max.overlaps = 30) {
  # set color names to the cols
  names(cols) <- path_list[, path_name]

  # modify the raw p-values to get -log10(p-value)
  if (raw_pvalue) {
    data <-
      data |>
      dplyr::mutate(log10pval = -log10(get(y)))
  } else {
    data <-
      data |>
      dplyr::mutate(log10pval = get(y))
  }

  # plot volcano plot with gray points and cut-offs
  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = log10pval)) +
    ggplot2::geom_point(color = "#CECACA", size = point_size) +
    ggplot2::geom_hline(
      yintercept = y_cutoff,
      linetype = "dashed",
      color = "#3E3B3B",
      linewidth = line_size
    ) +
    ggplot2::geom_vline(
      xintercept = x_cutoff,
      linetype = "dashed",
      color = "#3E3B3B",
      linewidth = line_size
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(
      labels = path_list[, "Description"],
      values = cols
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(face = "bold", size = 5),
      legend.text = ggplot2::element_text(size = 5),
      legend.position = "bottom",
      axis.title = ggplot2::element_text(face = "bold", size = 5)
    )

  # iteratively add colored points to the ggplot
  for (i in 1:nrow(path_list)) {
    gene_filter <-
      path_list[i, path_id] |>
      strsplit(",") |>
      unlist()
    data_x <-
      data |>
      dplyr::filter(get(id_col) %in% gene_filter) |>
      dplyr::mutate(color = rep(path_list[i, path_name]))
    plot <-
      plot +
      ggplot2::geom_point(
        data = data_x,
        ggplot2::aes(x = .data[[x]], y = log10pval, color = color),
        size = point_size + 1
      ) +
      ggrepel::geom_text_repel(
        data = data_x,
        ggplot2::aes(label = .data[[id_col]], color = color),
        size = text_size,
        force = 5,
        max.overlaps = max.overlaps,
        show.legend = FALSE,
        fontface = "bold"
      )
  }
  plot <-
    plot +
    ggplot2::labs(color = "Pathways")
  return(plot)
}

#' Publication-style heatmap of GSEA scores, saved to file
#'
#' Selects columns matching `column_pattern` (typically `"GScore"`), orders
#' rows within each category, and renders a [pheatmap::pheatmap()] to
#' `pwd/filename`.
#'
#' @param mat.heatmap A data frame with a `pathway_name` column, a
#'   `category_name` column, and one or more score columns matching
#'   `column_pattern`.
#' @param filename Output file name. Defaults to
#'   `"Integrated_Pathway_Heatmap.pdf"`.
#' @param low,high Optional numeric color-scale bounds; default to the 10th
#'   and 90th percentiles of the score columns.
#' @param pwd Directory the file should be saved to.
#' @param pathway_name Column in `mat.heatmap` used as row names. Defaults
#'   to `"Pathway"`.
#' @param category_name Column in `mat.heatmap` used for row annotation and
#'   ordering. Defaults to `"shared"`.
#' @param column_pattern Pattern used to select score columns. Defaults to
#'   `"GScore"`.
#' @param col_gradient RColorBrewer palette name. Defaults to `"RdYlBu"`.
#'
#' @return Invisibly, `NULL`. As a side effect, writes the heatmap file.
#' @export
enhancedHeatmap <- function(
  mat.heatmap,
  filename = "Integrated_Pathway_Heatmap.pdf",
  low = NULL,
  high = NULL,
  pwd,
  pathway_name = "Pathway",
  category_name = "shared",
  column_pattern = "GScore",
  col_gradient = "RdYlBu"
) {
  # Get the GScore column from the automated analysis or just the pattern
  sel_cols <- grep(column_pattern, colnames(mat.heatmap), value = TRUE)
  group_order <- dplyr::pull(mat.heatmap, dplyr::any_of(category_name)) |> unique()
  dt_graph <-
    mat.heatmap |>
    dplyr::mutate(dplyr::across(dplyr::all_of(sel_cols), as.numeric)) |>
    dplyr::mutate(!!category_name := factor(.data[[category_name]], levels = group_order)) |>
    dplyr::group_by(.data[[category_name]]) |>
    dplyr::arrange(.data[[category_name]], dplyr::across(dplyr::all_of(sel_cols), dplyr::desc), .by_group = TRUE) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(sel_cols))

  # Set the row names to the pathway names
  rownames(dt_graph) <- mat.heatmap[, pathway_name]

  # Set the column names to be the comparison names
  label <- stringr::str_extract(
    string = colnames(mat.heatmap),
    pattern = paste0("^(.*?)(?=", column_pattern, ")")
  ) |> na.omit()
  colnames(dt_graph) <- label

  # Get the category of comparison
  annot.row <- data.frame(Group = mat.heatmap[, category_name])
  rownames(annot.row) <- mat.heatmap[, pathway_name]

  # Get the color gradient to be scaled
  if (is.null(low)) low <- stats::quantile(dt_graph, probs = 0.1, na.rm = TRUE) # 10th percentile
  if (is.null(high)) high <- stats::quantile(dt_graph, probs = 0.9, na.rm = TRUE) # 90th percentile

  gradient <- seq(low, high, 0.1)
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
    n = 10,
    name = col_gradient
  )))(length(gradient))

  if (nrow(dt_graph) > 0) {
    pheatmap::pheatmap(
      mat = as.matrix(dt_graph),
      na_col = "white",
      border_color = "black",
      cellwidth = 15,
      cellheight = 15,
      filename = paste0(pwd, "/", filename),
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      scale = "none",
      fontsize_row = 14,
      breaks = gradient,
      color = colors,
      annotation_row = annot.row
    )
  }
  return(invisible(NULL))
}

#' ROC curve plot for a bootstrapped/iterated model
#'
#' Plots every iteration's ROC curve in gray and highlights the iteration
#' whose AUC is closest to the median AUC in red.
#'
#' @param boostModel A list-like object with an `auc` element (numeric
#'   vector of per-iteration AUCs) and an element named `name_roc`
#'   containing a data frame with columns `Iter`, `fpr`, `tpr`.
#' @param name_roc Name of the element in `boostModel` holding the ROC
#'   points. Defaults to `"pROC"`.
#'
#' @return A `ggplot` object.
#' @export
AdvancedROC <- function(boostModel, name_roc = "pROC") {
  maxIter <-
    which(abs(boostModel$auc - stats::median(boostModel$auc)) == min(abs(boostModel$auc - stats::median(boostModel$auc))))[1]
  bestRoc <-
    boostModel[[name_roc]] |>
    dplyr::filter(Iter == maxIter)
  plot <-
    boostModel[[name_roc]] |>
    ggplot2::ggplot(
      ggplot2::aes(x = fpr, y = tpr, group = Iter)
    ) +
    ggplot2::geom_line(linewidth = 0.5, color = "darkgray", lty = 2) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = 0.5,
      color = "black"
    ) +
    ggplot2::geom_line(
      data = bestRoc,
      ggplot2::aes(x = fpr, y = tpr),
      linewidth = 1,
      color = "red"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_bw(base_size = 14, base_rect_size = 2) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  return(plot)
}

#' Pathway heatmap saved to file (simplified variant)
#'
#' Takes a data table containing GScore-style columns and renders a
#' [pheatmap::pheatmap()], for use with previously computed GSEA graphing
#' parameters.
#'
#' @param data A data frame with a `Pathway` column and one or more columns
#'   containing `"GScore"` in their name.
#' @param filename Output file name. Defaults to
#'   `"customized pathways.pdf"`.
#' @param pwd Directory the file should be saved to.
#'
#' @return Invisibly, `NULL`. As a side effect, writes the heatmap file (if
#'   `data` has at least one row).
#' @export
pathwayHeatmap <- function(data, filename = "customized pathways.pdf", pwd) {
  dt_graph <-
    data |>
    dplyr::select(dplyr::contains("GScore")) |>
    apply(MARGIN = 2, FUN = as.numeric)
  rownames(dt_graph) <- data$Pathway
  label <- stringr::str_extract(string = colnames(dt_graph), pattern = "^(.*?)(?=GScore)") |> na.omit()
  colnames(dt_graph) <- label
  if (nrow(dt_graph) > 0) {
    pheatmap::pheatmap(
      mat = as.matrix(dt_graph),
      na_col = "white",
      border_color = "black",
      cellwidth = 15,
      cellheight = 15,
      filename = paste0(pwd, filename),
      cluster_cols = FALSE, cluster_rows = FALSE, scale = "none"
    )
  }
  return(invisible(NULL))
}

#' Circular ("Manhattan"-style) plot across chromosomes and levels
#'
#' @param data A data frame with columns for chromosome (`chr_name`),
#'   position (`pos_name`), p-value (`pval_name`), effect size
#'   (`beta_name`), and a `Level` column used for vertical stacking.
#' @param gap Ratio of spacing between chromosomes. Defaults to `1.1`.
#' @param tile_height,box_size,tile_width Plot geometry parameters.
#' @param low_col,high_col Colors for the diverging fill scale.
#' @param chr_name,pos_name,pval_name,beta_name Column names in `data`.
#'
#' @return A `ggplot` object.
#' @export
circleManhattan <- function(data, gap = 1.1, tile_height = 0.15, box_size = 0.5,
                             tile_width = 0.01, low_col = "blue", high_col = "red",
                             chr_name = "chr", pos_name = "pos", pval_name = "P.Value",
                             beta_name = "B") {
  height_gap <- tile_height * 0.6
  height_space <- tile_height * 1.15

  # Get the chromosome information
  chr_info <-
    data |>
    dplyr::group_by(get(chr_name)) |>
    dplyr::summarise(chr_len = max(get(pos_name))) |>
    dplyr::mutate(chr_start = cumsum(dplyr::lag(chr_len * gap, default = 0)))
  colnames(chr_info) <- c(chr_name, "chr_len", "chr_start")

  # Calculate -log10(p)
  data <-
    data |>
    dplyr::mutate(logp = -log10(get(pval_name))) |>
    # Merge the chromosome information into the overall data frame
    dplyr::left_join(chr_info, by = chr_name) |>
    dplyr::mutate(pos = get(pos_name) + chr_start) |>
    # Map position to angle (0 to 2*pi)
    dplyr::mutate(angle = 2 * pi * (pos - min(pos)) / (max(pos) - min(pos))) |>
    dplyr::mutate(Bscore = get(beta_name)) |> # Get the Beta x -log(pvalue) score
    dplyr::mutate(ypos = as.numeric(as.factor(Level)) * height_space + 1)

  # Prepare chromosome outlines
  chr_bounds <-
    data |>
    dplyr::group_by(get(chr_name)) |>
    dplyr::summarise(
      start_angle = min(angle),
      end_angle = max(angle),
      center_angle = (min(angle) + max(angle)) / 2
    )

  # Prepare the level labeling
  data_label <- data.frame(
    xpos = rep(0, length(unique(data$Level))),
    ypos = unique(data$ypos),
    labels = unique(data$Level)
  )

  # Prepare the chromosome outline y-pos
  ymin <- min(data$ypos - height_gap)
  ymax <- max(data$ypos + height_gap)

  # Plot: x = angle, y = -log10(p)
  p <-
    ggplot2::ggplot(data, ggplot2::aes(x = angle, y = ypos)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = Bscore, width = abs(Bscore) * tile_width),
      height = tile_height,
      alpha = 0.8
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_col,
      mid = "white",
      high = high_col,
      midpoint = 0
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, max(data$ypos) + 0.5),
      breaks = max(data$ypos) + 0.5
    ) +
    ggplot2::scale_x_continuous(limits = c(0, max(data$angle) + 0.4)) +
    ggplot2::geom_text(
      data = chr_bounds,
      ggplot2::aes(
        x = center_angle,
        y = ymax + 0.2,
        label = paste0("Chr ", unique(data[, chr_name]))
      ),
      size = 4
    ) +
    ggplot2::geom_text(
      data = data_label,
      ggplot2::aes(x = xpos, y = ypos, label = labels),
      size = 3.5,
      hjust = 1
    ) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::theme_void()

  # Add chromosome outlines as segments
  for (i in 1:nrow(chr_bounds)) {
    p <-
      p +
      ggplot2::annotate(
        geom = "path",
        x = c(
          chr_bounds$start_angle[i],
          chr_bounds$end_angle[i],
          chr_bounds$end_angle[i],
          chr_bounds$start_angle[i],
          chr_bounds$start_angle[i]
        ),
        y = c(ymin, ymin, ymax, ymax, ymin),
        color = "black",
        linewidth = box_size
      )
  }

  # Return the plot
  return(p)
}

#' Side-by-side ("staggered") volcano plots with highlighted genes
#'
#' @param data.list A named list of data frames, each with columns `GeneID`,
#'   `log2FC`, and the column named in `name.Padjust`.
#' @param col.list Named vector of colors keyed by highlight category.
#'   Defaults to `c("Highlighted" = "red")`.
#' @param highlight.list A named list (same names as `data.list`) of data
#'   frames with columns `GeneID` and `SubLabel` used to mark genes for
#'   highlighting/labeling.
#' @param ylimit Length-2 numeric, y-axis limits. Defaults to `c(-10, 10)`.
#' @param alpha_background Alpha for non-highlighted points.
#' @param pos Position adjustment passed to [ggplot2::geom_point()].
#' @param ylab Y-axis label.
#' @param name.Padjust Column name in `data.list` entries holding the
#'   adjusted p-value.
#' @param text.size Label text size.
#' @param max.overlap Passed to [ggrepel::geom_text_repel()].
#' @param order.display Optional character vector giving the display order
#'   of `Category` facets.
#' @param log2fc.cutoff.label Minimum absolute `log2FC` required for a point
#'   to be labeled.
#' @param pt.size Point size.
#'
#' @return A `ggplot` object, faceted by `Category`.
#' @export
staggeredVolcano <- function(data.list,
                              col.list = c("Highlighted" = "red"),
                              highlight.list,
                              ylimit = c(-10, 10),
                              alpha_background = 0.2,
                              pos = "identity",
                              ylab = "log2FC",
                              name.Padjust = "Padjusted",
                              text.size = 5,
                              max.overlap = 25,
                              order.display = NULL,
                              log2fc.cutoff.label = 0.2,
                              pt.size = 3) {
  # Combine the category data into each element of the list
  for (i in names(data.list)) {
    data.list[[i]] <-
      data.list[[i]] |>
      dplyr::left_join(highlight.list[[i]], by = "GeneID") |>
      dplyr::distinct(GeneID, .keep_all = TRUE)
  }

  # Combine and process data
  combined_data <-
    dplyr::bind_rows(data.list) |>
    dplyr::mutate(
      neg_log10_p = -log10(.data[[name.Padjust]]),
      label = ifelse(!is.na(SubLabel), GeneID, ""),
      alpha = dplyr::case_when(is.na(SubLabel) ~ alpha_background, TRUE ~ 1),
      col.border = dplyr::case_when(is.na(SubLabel) ~ NA, TRUE ~ "#000000")
    ) |>
    dplyr::arrange(!is.na(SubLabel), SubLabel)

  # Label the individual dot with nth percentile of dots
  combined_data <-
    combined_data |>
    dplyr::group_by(Category) |>
    dplyr::mutate(
      label = dplyr::case_when(
        abs(log2FC) > log2fc.cutoff.label ~ label,
        TRUE ~ NA
      )
    )

  # Arrange the order to display if the call is not NA
  if (!is.null(order.display) || length(order.display) > 0) {
    combined_data <-
      combined_data |>
      dplyr::mutate(Category = factor(Category, levels = order.display))
  }

  # Volcano plot with facets (side by side)
  plot <- ggplot2::ggplot(
    combined_data,
    ggplot2::aes(x = neg_log10_p, y = log2FC)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = SubLabel, alpha = alpha),
      color = "#000000",
      size = pt.size,
      position = pos,
      shape = 21
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey30"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      force = 5,
      force_pull = 2,
      size = text.size,
      max.overlaps = max.overlap
    ) +
    ggplot2::scale_fill_manual(values = col.list) +
    ggplot2::labs(x = expression(-log[10](P - value)), y = ylab) +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~Category, nrow = 1, scales = "free_x") + # Side by side in the same graph
    ggplot2::scale_y_continuous(limits = ylimit) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", size = 12)
    )
  return(plot)
}

#' Flow-diagram plot of QC step retention
#'
#' Create a flow-style diagram summarizing the number of features retained
#' at each QC step and the associated removal percentages.
#'
#' @param flow_df A data frame with at least the columns:
#'   \describe{
#'     \item{step_id}{Integer or numeric step index (1 = top).}
#'     \item{label}{Character label shown inside each box.}
#'     \item{removed_label}{Character label for text next to arrows; use
#'       `""` for steps without an arrow label.}
#'   }
#' @param space Vertical spacing between boxes. Defaults to `0.8`.
#'
#' @return A `ggplot` object representing the QC flow diagram.
#' @examples
#' flow_df <- data.frame(
#'   step_id = 1:3,
#'   label = c(
#'     "Total of 5421 proteins",
#'     "5349 proteins remain",
#'     "5095 of 5349 proteins\n>80% above lower limit of quantification"
#'   ),
#'   removed_label = c(
#'     "1.33% removed\nQC Process",
#'     "4.75% removed\n>20% samples < LOD",
#'     ""
#'   ),
#'   stringsAsFactors = FALSE
#' )
#'
#' p <- olink_qc_graph(flow_df)
#' print(p)
#'
#' @export
olink_qc_graph <- function(flow_df, space = 0.8) {
  # get number of steps
  n_steps <- nrow(flow_df)

  # x/y position of box centers
  flow_df$x <- 0
  flow_df$y <- seq(from = 0, by = -space, length.out = n_steps)

  # data for arrows: from step i to i + 1
  arrow_df <- flow_df[-n_steps, , drop = FALSE]
  arrow_df$xend <- flow_df$x[-1]
  arrow_df$yend <- flow_df$y[-1]

  # midpoints for arrow text
  arrow_df$xtxt <- (arrow_df$x + arrow_df$xend) / 2 + 0.2
  arrow_df$ytxt <- (arrow_df$y + arrow_df$yend) / 2

  p <- ggplot2::ggplot() +
    # boxes
    ggplot2::geom_label(
      data = flow_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      label.size = 0.8,
      label.r = grid::unit(0.35, "lines"),
      label.padding = grid::unit(0.6, "lines"),
      size = 4,
      fill = "white",
      color = "black"
    ) +
    # arrows from box to box
    ggplot2::geom_segment(
      data = arrow_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y - space * 0.2, # start just below upper box
        xend = .data$xend,
        yend = .data$yend + space * 0.2 # end just above lower box
      ),
      arrow = grid::arrow(length = grid::unit(0.25, "cm")),
      linewidth = 0.6
    ) +
    # arrow labels between boxes, only where non-empty
    ggplot2::geom_text(
      data = arrow_df[arrow_df$removed_label != "", , drop = FALSE],
      ggplot2::aes(
        x = .data$xtxt,
        y = .data$ytxt,
        label = .data$removed_label
      ),
      hjust = 0,
      size = 3.5
    ) +
    ggplot2::coord_cartesian(xlim = c(-3, 3), ylim = c(-n_steps * space * 0.8, 0)) +
    ggplot2::theme_void()

  return(p)
}
