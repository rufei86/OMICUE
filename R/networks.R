#' Build protein co-expression networks for multiple modules
#'
#' For each specified module, performs hierarchical clustering and dynamic
#' tree cut on `prox_avg`, then builds a weighted network and returns the
#' network, node degrees, and scaled edge widths.
#'
#' @param prox_avg Numeric proximity matrix (square, symmetric, values in
#'   `[0, 1]`).
#' @param module_ids Vector of integers, modules to output networks for
#'   (e.g., `c(2, 5, 7)`).
#' @param minClusterSize Minimum cluster size for
#'   [dynamicTreeCut::cutreeDynamic()] (default: 10).
#' @param method_clust [stats::hclust()] method (default: `"ward.D2"`).
#' @param deepSplit `deepSplit` parameter for
#'   [dynamicTreeCut::cutreeDynamic()] (default: 4).
#' @param edge_width_scale Numeric scaling of edge weights for plotting
#'   (default: 5).
#'
#' @return List containing:
#'   - `hc`: clustering object
#'   - `modules`: cluster assignment vector
#'   - `data_modules`: ProteinID/module mapping
#'   - `networks`: named list of igraph objects (per module)
#'   - `degree`: named list of degree vectors (per module)
#'   - `edge_widths`: named list of scaled edge widths (per module)
#' @export
build_protein_network_multi <- function(
  prox_avg,
  module_ids,
  minClusterSize = 10,
  method_clust = "ward.D2",
  deepSplit = 4,
  edge_width_scale = 5
) {
  # Step 1: Clustering & module assignments
  distobj <- 1 - prox_avg |> stats::as.dist()
  hc <- stats::hclust(distobj, method = method_clust)
  modules <- dynamicTreeCut::cutreeDynamic(
    hc,
    distM = as.matrix(distobj),
    minClusterSize = minClusterSize,
    method = "hybrid",
    deepSplit = deepSplit
  )
  data_modules <- data.frame(
    ProteinID = rownames(prox_avg),
    Module = modules
  )

  # Step 2: Process each module
  networks <- list()
  degree <- list()
  edge_widths <- list()

  for (mid in module_ids) {
    sel_protein_ids <- data_modules |>
      dplyr::filter(Module == mid) |>
      (\(df) df$ProteinID)()
    if (length(sel_protein_ids) >= 2) { # Only build network if at least 2 proteins
      adj_mat <- prox_avg[
        rownames(prox_avg) %in% sel_protein_ids,
        colnames(prox_avg) %in% sel_protein_ids
      ]
      diag(adj_mat) <- 0

      g <- igraph::graph_from_adjacency_matrix(
        adj_mat,
        weighted = TRUE,
        diag = FALSE,
        mode = "undirected"
      )
      degs <- igraph::degree(g)
      # For plotting: scale edge widths based on weight
      sw <- igraph::E(g)$weight * edge_width_scale
      networks[[as.character(mid)]] <- g
      degree[[as.character(mid)]] <- degs
      edge_widths[[as.character(mid)]] <- sw
    } else {
      networks[[as.character(mid)]] <- NULL
      degree[[as.character(mid)]] <- NULL
      edge_widths[[as.character(mid)]] <- NULL
    }
  }

  return(list(
    hc = hc,
    modules = modules,
    data_modules = data_modules,
    networks = networks,
    degree = degree,
    edge_widths = edge_widths
  ))
}
