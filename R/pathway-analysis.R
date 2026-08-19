#' Run differential-feature comparisons across groups
#'
#' For every pairwise (or reference-vs-other) combination of levels of
#' `factor`, runs the requested per-feature test (Wilcoxon, logistic,
#' linear/gaussian regression, or rank regression), annotates significant
#' features with Entrez IDs, and writes summary CSVs to disk.
#'
#' @param factor Character scalar, column name in `data$metadata` used to
#'   define groups (or the response, for regression methods).
#' @param omic_name Character scalar, name of the omic table in `data`.
#' @param level_id Optional character scalar, random-effect grouping column
#'   passed through to [LogReg()].
#' @param name_suffix Character scalar appended to the output directory
#'   name.
#' @param method One of `"wilcox"`, `"regression"`, `"rankregression"`,
#'   `"logistic"`.
#' @param empirical Logical; use permutation p-values (methods that support
#'   it).
#' @param n_perm Number of permutations when `empirical = TRUE`.
#' @param covariates Character vector of covariate column names, or `NA`.
#' @param reference Optional character scalar; if supplied, all other levels
#'   of `factor` are compared against this reference level instead of every
#'   pairwise combination.
#' @param id_name Character scalar, sample ID column shared between
#'   `data[[omic_name]]` and `data$metadata`.
#' @param data A list-like object containing `metadata` and one or more omic
#'   tables.
#' @param gene_db An annotation database object (e.g. an `OrgDb`), passed to
#'   [AnnotationDbi::select()].
#' @param p_adjust_method P-value adjustment method for the per-comparison
#'   feature tests, passed to [stats::p.adjust()].
#' @param p_gsea_adjust_method P-value adjustment method reserved for
#'   downstream GSEA steps (kept for interface parity with
#'   [Enrich_GO()]/[Enrich_KEGG()]/[Enrich_REACTOME()]).
#' @param min_sig Minimum gene set size reserved for downstream GSEA steps.
#' @param species Species code reserved for downstream GSEA steps.
#' @param plimit P-value cutoff reserved for downstream GSEA steps.
#' @param plimit_test P-value cutoff used to filter significant features
#'   from the per-comparison tests.
#' @param oid_suffix Character scalar, suffix separating a feature's gene ID
#'   from its assay/probe identifier (e.g. `"-OID"`).
#' @param n_cores Number of cores to use for regression-based methods.
#' @param pwd Character scalar, directory in which output subfolders and
#'   CSVs are written.
#'
#' @return Invisibly, a list with `sig_gene_all` and `de_all`, the per-
#'   comparison significant-feature and full differential-expression tables.
#'   As a side effect, writes `pwd/<factor><name_suffix> - Pathway
#'   Analysis/` with subfolders and summary CSVs.
#' @export
omicue_de <- function(
  factor,
  omic_name,
  level_id = NULL,
  name_suffix = "",
  method = c("wilcox", "regression", "rankregression", "logistic"),
  empirical = FALSE,
  n_perm = 1000,
  covariates = NA,
  reference = NULL,
  id_name,
  data,
  gene_db,
  p_adjust_method = "none",
  p_gsea_adjust_method = "fdr",
  min_sig = 20,
  species = "hsa",
  plimit = 0.05,
  plimit_test = 0.05,
  oid_suffix = "-OID",
  n_cores = round((parallel::detectCores() - 6) / 2),
  pwd
) {
  method <- match.arg(method)

  # running through the loop
  sig_gene_all <- list()
  de_all <- list()

  # make sub-directory for each cell type
  j <- paste0(factor, name_suffix, " - Pathway Analysis")
  dir.create(path = paste0(pwd, "/", j))
  dirs.new <- c(
    volcano = paste0(pwd, "/", j, "/Volcano plots/"),
    gonet = paste0(pwd, "/", j, "/GO BP Network/"),
    pa = paste0(pwd, "/", j, "/Reactome Pathways/"),
    keggpath = paste0(pwd, "/", j, "/Kegg Pathways/"),
    keggref = paste0(pwd, "/", j, "/Kegg Pathways/Reference Pathway/")
  )
  lapply(X = dirs.new, function(i) dir.create(path = i))

  # filter down to the samples that are actually in the dataset
  meta_refine <-
    data[[omic_name]] |>
    dplyr::select(dplyr::all_of(id_name)) |>
    dplyr::left_join(data$metadata, by = id_name) |>
    tibble::tibble()

  # get all possible combination of the factors
  # first grab the unique categories within the factor variable
  if (!method %in% c("regression", "rankregression")) {
    categories <-
      meta_refine[, factor] |>
      unique() |>
      na.omit() |>
      dplyr::pull(dplyr::all_of(factor)) |>
      sort()
    if (is.null(reference)) {
      combn <- utils::combn(categories, m = 2, simplify = TRUE)
    } else {
      combn <- cbind(
        reference,
        setdiff(categories, reference)
      ) |>
        t()
    }
  } else {
    combn <- data.frame(
      factor = c(factor, factor),
      factor1 = c(paste0(factor, "2"), factor)
    )
  }

  for (i in c(1:ncol(combn))) {
    print(paste("Analyzing ", paste(combn[, i], collapse = " vs ")))
    # find the significantly expressed genes in each cell type
    data_i <-
      data[[omic_name]] |>
      dplyr::left_join(meta_refine, by = id_name) |>
      as.data.frame()

    if (!method %in% c("regression", "rankregression")) {
      data_i <-
        data_i |>
        dplyr::filter(get(factor) %in% combn[, i])
    }

    # get the comparison first
    comparison <- combn[, i]
    protein_names <-
      data_i |>
      dplyr::select(dplyr::contains(oid_suffix)) |>
      colnames()

    # Generate the column for the comparison
    prefix <- paste0(comparison[1], "_vs_", comparison[2], "_")

    # doing wilcoxon test on the pairing, or a regression-based method
    if (method == "wilcox") {
      degs <- simpleWilcox(protein_names, data_i, factor, comparison)
    } else if (method == "logistic") {
      degs <- LogReg(
        protein_names = protein_names,
        data_i = data_i,
        factor = factor,
        oid_suffix = oid_suffix,
        comparison = comparison,
        covariates = covariates,
        empirical = empirical,
        n_perm = n_perm,
        n_cores = n_cores,
        family = stats::binomial(),
        level_id = level_id,
        type = "logistic"
      )
    } else if (method == "regression") {
      degs <- LogReg(
        protein_names,
        data_i,
        factor,
        comparison = comparison,
        oid_suffix = oid_suffix,
        covariates = covariates,
        empirical = empirical,
        n_perm = n_perm,
        n_cores = n_cores,
        family = stats::gaussian(),
        level_id = level_id,
        type = "regression"
      )
    } else if (method == "rankregression") {
      degs <- RankReg(protein_names, data_i, factor, covariates)
    }

    # Adding in the column about the comparison
    degs <- degs |> dplyr::mutate(Comparison = prefix)

    # Quick inspection
    print(utils::head(
      degs |> dplyr::filter(Pvalue < plimit_test) |> dplyr::arrange(-Log2FC)
    ))

    # filter the significant genes
    sig_gene <-
      degs |>
      dplyr::distinct(GeneID, .keep_all = TRUE) |> # filter out the duplicated geneID
      dplyr::mutate(
        Padjusted = stats::p.adjust(p = Pvalue, method = p_adjust_method)
      ) |>
      dplyr::filter(Padjusted < plimit_test) |> # filter the p-value < 0.05
      dplyr::arrange(-Log2FC)
    sig_gene <- AnnotationDbi::select(
      gene_db,
      keys = sig_gene$GeneID,
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "SYMBOL"
    ) |>
      dplyr::rename(GeneID = SYMBOL) |>
      dplyr::left_join(sig_gene, by = "GeneID")

    # Pushing the iterative results into the list
    sig_gene_all[[i]] <- sig_gene
    de_all[[i]] <- degs

    # if too few of genes to run then skip
    print(paste("Number of significant genes = ", nrow(sig_gene)))
  }
  # Output the concatenated table of significant DEGs across all comparisons
  utils::write.csv(
    x = dplyr::bind_rows(sig_gene_all),
    file = paste0(pwd, "/", j, "/sig DEGs across all comparisons.csv")
  )
  # Output the concatenated table of all DEGs (significant or not) across all comparisons
  utils::write.csv(
    x = dplyr::bind_rows(de_all),
    file = paste0(pwd, "/", j, "/all DEGs across all comparisons.csv")
  )

  return(invisible(list(sig_gene_all = sig_gene_all, de_all = de_all)))
}
