#' Build a gene-list summary table for a set of pathway IDs
#'
#' For each ID in `ids`, looks up the matching gene set, its description in
#' `gse_res`, and the subset of genes present in `gage_in`, returning a
#' single tidy table.
#'
#' @param ids Character vector of pathway IDs (matched against
#'   `names(gset)`); `NA` entries are dropped.
#' @param gset A named list of gene sets (e.g. `geneSets` slot from a
#'   clusterProfiler / ReactomePA result).
#' @param gse_res A data frame with columns `Pathway` and `PathID` used to
#'   look up the human-readable description for each ID.
#' @param gage_in A named numeric vector (e.g. a ranked gene list) whose
#'   names are used to restrict each pathway's genes to those present.
#' @param gene_db An annotation database object (e.g. an `OrgDb`), passed to
#'   [AnnotationDbi::mapIds()].
#'
#' @return A data frame with columns `Pathway.ID`, `Description`, `Genes`
#'   (comma-separated gene symbols).
#' @export
pathGeneList <- function(ids, gset, gse_res, gage_in, gene_db) {
  tbl <- list()
  for (k in na.omit(ids)) {
    path.name <- names(gset[grep(names(gset), pattern = k)])
    description <- gse_res$Pathway[grep(pattern = k, x = gse_res$PathID)]
    path <- gset[grep(names(gset), pattern = k)][[1]]
    pathway.genes <- path[path %in% names(gage_in)]
    path.gene.name <- suppressMessages(
      AnnotationDbi::mapIds(gene_db, pathway.genes, "SYMBOL", "ENTREZID")
    )
    tbl[[k]] <- data.frame(
      `Pathway ID` = path.name,
      Description = description,
      Genes = paste0(path.gene.name, collapse = ",")
    )
  }
  return(dplyr::bind_rows(tbl))
}


#' Gene Set Enrichment Analysis on GO Biological Process terms
#'
#' Wraps [clusterProfiler::gseGO()] (ontology `"BP"`) with a `tryCatch` that
#' retries once on warnings, then simplifies redundant terms and returns a
#' tidy summary table.
#'
#' @param sigs A named, sorted numeric vector (e.g. ranked log-fold-changes)
#'   as expected by `geneList` in [clusterProfiler::gseGO()].
#' @param gene_db An `OrgDb` annotation object.
#' @param plimit P-value cutoff passed to `pvalueCutoff`.
#' @param p_adjust_method P-value adjustment method passed to
#'   `pAdjustMethod`.
#'
#' @return A list with elements `gse_GO` (the raw enrichResult), `result_GO`
#'   (tidy summary data frame with `PathID`, `Pathway`, `Score`, `Pval`,
#'   `GScore`), `go_gset` (the underlying gene sets) and `go_ids`.
#' @export
Enrich_GO <- function(sigs, gene_db, plimit, p_adjust_method) {
  # geneset enrichment assay with GO terms
  gse_GO <-
    tryCatch(
      expr = {
        clusterProfiler::gseGO(
          geneList = sigs,
          pAdjustMethod = p_adjust_method,
          OrgDb = gene_db,
          ont = "BP",
          pvalueCutoff = plimit, # set cut-off at p < 0.05
          minGSSize = 10,
          verbose = TRUE,
          by = "fgsea"
        )
      },
      error = function(e) {
        print("There was an error generating GO GSEA")
        NULL
      },
      warning = function(w) {
        gse_GO <- clusterProfiler::gseGO(
          geneList = sigs,
          pAdjustMethod = p_adjust_method,
          OrgDb = gene_db,
          ont = "BP",
          pvalueCutoff = plimit,
          minGSSize = 10,
          verbose = TRUE,
          by = "fgsea"
        )
        print("There was a warning message but the GSEA is complete")
        return(gse_GO)
      }
    )

  print(paste(length(gse_GO$ID), " of GO pathways detected"))
  if (length(gse_GO$ID) > 0) {
    gse_GO <- clusterProfiler::simplify(x = gse_GO, cutoff = 0.6)
    go_gset <- gse_GO@geneSets
    result_GO <- data.frame(
      PathID = gse_GO$ID,
      Pathway = paste0("GO: ", gse_GO$Description),
      Score = gse_GO$enrichmentScore,
      Pval = gse_GO$pvalue,
      GScore = gse_GO$enrichmentScore * -log10(gse_GO$pvalue)
    )
    go_ids <- result_GO$PathID
  } else {
    result_GO <- data.frame(matrix(ncol = 5))
    colnames(result_GO) <- c("PathID", "Pathway", "Score", "Pval", "GScore")
    go_ids <- NULL
    go_gset <- NULL
  }
  return(list(gse_GO = gse_GO, result_GO = result_GO, go_gset = go_gset, go_ids = go_ids))
}

#' Gene Set Enrichment Analysis on KEGG terms
#'
#' Wraps [clusterProfiler::gseKEGG()] (human, `organism = "hsa"`) with a
#' `tryCatch` that retries once on warnings, and returns a tidy summary
#' table.
#'
#' @param sigs A named, sorted numeric vector as expected by `geneList` in
#'   [clusterProfiler::gseKEGG()].
#' @param plimit P-value cutoff passed to `pvalueCutoff`.
#' @param p_adjust_method P-value adjustment method passed to
#'   `pAdjustMethod`.
#'
#' @return A list with elements `gse_KEGG`, `result_Kegg`, `kegg_gset`,
#'   `kegg_ids`.
#' @export
Enrich_KEGG <- function(sigs, plimit, p_adjust_method) {
  gse_KEGG <-
    tryCatch(
      expr = {
        clusterProfiler::gseKEGG(
          geneList = sigs,
          verbose = TRUE,
          pAdjustMethod = p_adjust_method,
          organism = "hsa",
          minGSSize = 10,
          pvalueCutoff = plimit,
          by = "fgsea"
        )
      },
      error = function(e) {
        print("There was an error generating KEGG GSEA")
        NULL
      },
      warning = function(w) {
        gse_KEGG <- clusterProfiler::gseKEGG(
          geneList = sigs,
          verbose = TRUE,
          pAdjustMethod = p_adjust_method,
          organism = "hsa",
          minGSSize = 10,
          pvalueCutoff = plimit,
          by = "fgsea"
        )
        print("There was a warning message but the GSEA is complete")
        return(gse_KEGG)
      }
    )

  print(paste(length(gse_KEGG$ID), " of Kegg pathways detected"))
  if (length(gse_KEGG$ID) > 0) {
    kegg_gset <- gse_KEGG@geneSets
    result_Kegg <- data.frame(
      PathID = gse_KEGG$ID,
      Pathway = paste0("KEGG: ", gse_KEGG$Description),
      Score = gse_KEGG$enrichmentScore,
      Pval = gse_KEGG$pvalue,
      GScore = gse_KEGG$enrichmentScore * -log10(gse_KEGG$pvalue)
    )
    kegg_ids <- result_Kegg$PathID
  } else {
    result_Kegg <- data.frame(matrix(ncol = 5))
    colnames(result_Kegg) <- c("PathID", "Pathway", "Score", "Pval", "GScore")
    kegg_ids <- NULL
    kegg_gset <- NULL
  }
  return(list(gse_KEGG = gse_KEGG, result_Kegg = result_Kegg, kegg_gset = kegg_gset, kegg_ids = kegg_ids))
}

#' Gene Set Enrichment Analysis on Reactome terms
#'
#' Wraps [ReactomePA::gsePathway()] with a `tryCatch` that retries once on
#' warnings, and returns a tidy summary table.
#'
#' @param sigs A named, sorted numeric vector as expected by `geneList` in
#'   [ReactomePA::gsePathway()].
#' @param plimit P-value cutoff passed to `pvalueCutoff`.
#' @param p_adjust_method P-value adjustment method passed to
#'   `pAdjustMethod`.
#'
#' @return A list with elements `gse_PA`, `result_pa`, `pa_gset`, `pa_ids`.
#' @export
Enrich_REACTOME <- function(sigs, plimit, p_adjust_method) {
  gse_PA <-
    tryCatch(
      expr = {
        ReactomePA::gsePathway(
          geneList = sigs,
          pAdjustMethod = p_adjust_method,
          verbose = TRUE,
          pvalueCutoff = plimit,
          minGSSize = 10,
          by = "fgsea"
        )
      },
      error = function(e) {
        print("There was an error generating REACTOME GSEA")
        NULL
      },
      warning = function(w) {
        gse_PA <- ReactomePA::gsePathway(
          geneList = sigs,
          pAdjustMethod = p_adjust_method,
          verbose = TRUE,
          pvalueCutoff = plimit,
          minGSSize = 10,
          by = "fgsea"
        )
        print("There was a warning message but the GSEA is complete")
        return(gse_PA)
      }
    )

  # graph Reactome pathways
  if (is.null(gse_PA) || !is.data.frame(gse_PA@result) || nrow(gse_PA@result) == 0) {
    print("0 of Reactome pathways detected")
    result_pa <- data.frame(matrix(ncol = 5))
    colnames(result_pa) <- c("PathID", "Pathway", "Score", "Pval", "GScore")
    pa_ids <- NULL
    pa_gset <- NULL
  } else {
    print(paste(nrow(gse_PA@result), " of Reactome pathways detected"))
    pa_gset <- gse_PA@geneSets
    result_pa <- data.frame(
      PathID = gse_PA$ID,
      Pathway = paste0("Reactome: ", gse_PA$Description),
      Score = gse_PA$enrichmentScore,
      Pval = gse_PA$pvalue,
      GScore = gse_PA$enrichmentScore * -log10(gse_PA$pvalue)
    )
    pa_ids <- result_pa$PathID
  }

  return(list(gse_PA = gse_PA, result_pa = result_pa, pa_gset = pa_gset, pa_ids = pa_ids))
}
