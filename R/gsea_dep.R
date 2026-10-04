# ============================================================================
# GSEA across multiple gene-set collections (KEGG, Reactome, GO:BP, Hallmark)
# BMS vs RRMS and BMS vs SPMS
# ============================================================================

# ---- Gene-set retrieval ----------------------------------------------------

#' Fetch an MSigDB collection as a TERM2GENE data frame (Entrez IDs)
get_msig_t2g <- function(
  collection,
  subcollection = NULL,
  species = "Homo sapiens"
) {
  # msigdbr >= 10 uses collection/subcollection; older versions category/subcategory
  args_new <- list(species = species, collection = collection)
  if (!is.null(subcollection)) {
    args_new$subcollection <- subcollection
  }

  db <- tryCatch(
    do.call(msigdbr::msigdbr, args_new),
    error = function(e) {
      args_old <- list(species = species, category = collection)
      if (!is.null(subcollection)) {
        args_old$subcategory <- subcollection
      }
      do.call(msigdbr::msigdbr, args_old)
    }
  )

  gs_col <- intersect(c("gs_name"), names(db))[1]
  gene_col <- intersect(c("gene_symbol"), names(db))[1]

  db |>
    dplyr::select(
      gs_name = dplyr::all_of(gs_col),
      entrez_gene = dplyr::all_of(gene_col)
    ) |>
    dplyr::mutate(entrez_gene = as.character(entrez_gene)) |>
    dplyr::filter(!is.na(entrez_gene)) |>
    dplyr::distinct()
}

#' Collections to test (name -> MSigDB collection/subcollection)
gs_collections <- list(
  Hallmark = list(collection = "H", subcollection = NULL),
  KEGG = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  GO_BP = list(collection = "C5", subcollection = "GO:BP")
)

tg <- get_msig_t2g(collection = "C2", subcollection = "CP:KEGG_LEGACY")

t2g_list <- purrr::map(gs_collections, function(x) {
  t2g <- tryCatch(
    get_msig_t2g(x$collection, x$subcollection),
    error = function(e) NULL
  )
  # KEGG_LEGACY only exists in newer MSigDB releases; fall back to "CP:KEGG"
  if (is.null(t2g) && identical(x$subcollection, "CP:KEGG_LEGACY")) {
    t2g <- get_msig_t2g("C2", "CP:KEGG")
  }
  t2g
})

# ---- Ranking metric --------------------------------------------------------

#' Named, strictly sorted ranking vector: -log10(Pvalue) * log2FC
build_ranked_list <- function(df) {
  ranked <- df |>
    dplyr::mutate(rank_metric = dis_coef) |>
    dplyr::filter(is.finite(rank_metric)) |>
    dplyr::arrange(dplyr::desc(rank_metric), GeneID) |>
    dplyr::distinct(GeneID, .keep_all = TRUE)

  gene_list <- ranked$rank_metric
  names(gene_list) <- as.character(ranked$GeneID)
  gene_list
}

# ---- GSEA runners ----------------------------------------------------------

#' Run GSEA for one comparison against one TERM2GENE table
run_gsea_ranked <- function(df, term2gene, adjust_method = "fdr", seed = 7419) {
  gene_list <- build_ranked_list(df)
  set.seed(seed)
  clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    eps = 0,
    seed = TRUE,
    pAdjustMethod = adjust_method,
    verbose = FALSE
  )
}

#' Run GSEA for every comparison in a list of DE tables
run_gsea_collection <- function(
  dep_list,
  term2gene,
  adjust_method = "fdr",
  seed = 7419
) {
  purrr::imap(
    dep_list,
    ~ run_gsea_ranked(.x, term2gene, adjust_method = adjust_method, seed = seed)
  )
}

#' Collapse a list of gseaResult objects into one long tibble
tidy_gsea_results <- function(gsea_list, collection) {
  purrr::imap(
    gsea_list,
    ~ tibble::as_tibble(.x@result) |>
      dplyr::mutate(Comparison = .y, Collection = collection)
  ) |>
    dplyr::bind_rows()
}

#' Shorten MSigDB set names for plotting
clean_pathway_label <- function(x, width = 55) {
  x |>
    stringr::str_remove("^(HALLMARK|KEGG|KEGG_MEDICUS|REACTOME|GOBP|GO)_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_sentence() |>
    stringr::str_trunc(width)
}

# ---- Pathway selection -----------------------------------------------------

#' Pathways significant in at least one comparison
select_sig_pathways <- function(df, fdr = 0.05) {
  df |>
    dplyr::filter(p.adjust < fdr) |>
    dplyr::pull(ID) |>
    unique()
}

#' Top-n pathways by best (smallest) FDR across comparisons, |NES| as tiebreak
select_top_pathways <- function(df, n = 25) {
  df |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      best_padj = min(p.adjust, na.rm = TRUE),
      max_nes = max(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(best_padj, dplyr::desc(max_nes)) |>
    dplyr::slice_head(n = n) |>
    dplyr::pull(ID)
}

# ---- Heatmap ---------------------------------------------------------------

#' Wide matrix of one column per comparison
gsea_matrix <- function(df, pathways, value_col, comparisons) {
  m <- df |>
    dplyr::filter(ID %in% pathways) |>
    dplyr::select(ID, Comparison, dplyr::all_of(value_col)) |>
    tidyr::pivot_wider(
      names_from = Comparison,
      values_from = dplyr::all_of(value_col)
    ) |>
    tibble::column_to_rownames("ID") |>
    as.matrix()
  m[, comparisons[comparisons %in% colnames(m)], drop = FALSE]
}

#' NES heatmap; asterisks only for cells with FDR < 0.05
plot_gsea_heatmap <- function(df, pathways, comparisons, title, fdr = 0.05) {
  nes_mat <- gsea_matrix(df, pathways, "NES", comparisons)
  padj_mat <- gsea_matrix(df, pathways, "p.adjust", comparisons)
  padj_mat <- padj_mat[rownames(nes_mat), colnames(nes_mat), drop = FALSE]

  sig_label <- matrix(
    dplyr::case_when(
      is.na(as.vector(padj_mat)) ~ "",
      as.vector(padj_mat) >= fdr ~ "", # non-significant: no asterisk
      as.vector(padj_mat) < 0.001 ~ "***",
      as.vector(padj_mat) < 0.01 ~ "**",
      TRUE ~ "*"
    ),
    nrow = nrow(padj_mat),
    dimnames = dimnames(padj_mat)
  )

  ComplexHeatmap::Heatmap(
    nes_mat,
    name = "NES",
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    col = circlize::colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027")),
    na_col = "grey85",
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "average",
    row_labels = clean_pathway_label(rownames(nes_mat)),
    row_names_gp = grid::gpar(fontsize = 8),
    cluster_columns = FALSE,
    column_labels = stringr::str_remove(colnames(nes_mat), "_$"),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- sig_label[i, j]
      if (nzchar(lab)) {
        grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 10))
      }
    },
    border = TRUE,
    border_gp = grid::gpar(col = "black", lwd = 1),
    rect_gp = grid::gpar(col = "black", lwd = 0.5),
    width = ncol(nes_mat) * grid::unit(5, "mm"), # Heatmap BODY dimensions: same per-cell size in both axes
    height = nrow(nes_mat) * grid::unit(5, "mm")
  )
}
