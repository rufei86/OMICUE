# OMICUE

Tools for multi-omic differential analysis, enrichment, and visualization —
covering data wrangling, differential-feature testing (Wilcoxon, GLM,
mixed-effects, rank regression), pathway/enrichment analysis (GO, KEGG,
Reactome), protein co-expression networks, and a set of publication-style
plotting functions (volcano plots, heatmaps, ROC curves, Manhattan plots,
QC diagrams, and more).

All functions call dependencies with explicit `package::function()` syntax
(no `library()`/`require()` calls anywhere in the package code), so it is
always clear which package a given function comes from.

## Installation

`OMICUE` depends on three Bioconductor-only packages
(`AnnotationDbi`, `clusterProfiler`, `ReactomePA`) that are **not** on CRAN.
`remotes::install_github()` only installs CRAN dependencies automatically,
so install the Bioconductor packages first, then install `OMICUE`:

```r
# 1. Install Bioconductor dependencies (one-time)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("AnnotationDbi", "clusterProfiler", "ReactomePA"))

# 2. Install OMICUE from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("<your-github-username>/OMICUE")
```

Replace `<your-github-username>` with wherever this repository ends up
being hosted.

## Function reference

### Data manipulation (`R/manipulation.R`)

| Function | Description |
|---|---|
| `omicFilter()` | Subsets metadata and one or more omic tables to a chosen set of factor levels (with optional per-subject de-duplication). |
| `omicDeltaByVisit()` | Computes per-subject deltas from a baseline visit across one or more omic tables. |

### Differential expression (`R/differential-expression.R`)

| Function | Description |
|---|---|
| `simpleWilcox()` | Two-group Wilcoxon rank-sum test across a set of features. |
| `glm_model()` | Single-covariate generalized linear model (e.g., logistic regression). |
| `glm_model_mcore()` | Parallelized (multi-core) version of `glm_model()` across many features. |
| `LogReg()` | Logistic regression wrapper with optional covariates. |
| `GaussianReg()` | Linear (Gaussian) regression wrapper with optional covariates. |
| `rg_model()` | Linear mixed-effects model (via `lme4`/`lmerTest`) for repeated-measures designs. |
| `RankReg()` | Rank-based regression (via `Rfit`), robust to outliers/non-normality. |

### Pathway & enrichment analysis (`R/pathway-analysis.R`, `R/enrichment.R`)

| Function | Description |
|---|---|
| `PathwayAnalysis()` | End-to-end wrapper that runs differential testing, then GO/KEGG/Reactome enrichment on the resulting hit list. |
| `pathGeneList()` | Builds a named, ranked gene list suitable for GSEA-style enrichment functions. |
| `Enrich_GO()` | Gene Ontology enrichment via `clusterProfiler`. |
| `Enrich_KEGG()` | KEGG pathway enrichment via `clusterProfiler`. |
| `Enrich_REACTOME()` | Reactome pathway enrichment via `ReactomePA`. |

### Networks (`R/networks.R`)

| Function | Description |
|---|---|
| `build_protein_network_multi()` | Clusters a protein proximity matrix (hierarchical clustering + dynamic tree cut) and builds an `igraph` network, node degree, and scaled edge widths for each requested module. |

### Plotting (`R/plotting.R`)

| Function | Description |
|---|---|
| `graph_mean_var()` | Mean-variance relationship plot, useful for QC of omic intensity data. |
| `customCorrplot()` | Annotated correlation heatmap. |
| `Plot_QC()` | Multi-panel QC diagnostic plot. |
| `euler_venn()` | Euler/Venn diagram of overlapping feature sets. |
| `AdvancedVolcano()` | Volcano plot with flexible labeling and thresholding. |
| `enhancedHeatmap()` | Annotated heatmap with row/column clustering. |
| `AdvancedROC()` | ROC curve(s) with AUC annotation. |
| `pathwayHeatmap()` | Heatmap of enrichment results across conditions/contrasts. |
| `circleManhattan()` | Circular Manhattan-style plot across genomic/feature groups. |
| `staggeredVolcano()` | Volcano plot faceted/staggered across multiple contrasts. |
| `olink_qc_graph()` | QC plot tailored to Olink proteomics output. |

## Development notes

- Documentation (`man/*.Rd`) and `NAMESPACE` are generated with
  [`roxygen2`](https://roxygen2.r-lib.org/) (v8.1.0) from the `@`-tagged
  comments in each `R/*.R` file. After editing any function's roxygen
  comments, regenerate both with:

  ```r
  roxygen2::roxygenise()
  ```

- Unit tests (`tests/testthat/`) cover the dependency-light functions
  (`omicFilter()`, `omicDeltaByVisit()`, `simpleWilcox()`, `glm_model()`)
  that only require base R plus `dplyr`/`tidyr`/`stringr`. Run with:

  ```r
  devtools::test()
  ```

## License

MIT © Rufei Lu. See [`LICENSE.md`](LICENSE.md) for the full text.
