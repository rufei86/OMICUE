#' Compute visit-specific change from baseline for multiple omics
#'
#' For each omic table in `data`, computes the change from baseline
#' (`value - baseline_value`) for every numeric feature, at every non-baseline
#' visit, and returns the result in wide format (one row per subject, one
#' column per `delta-<visit>-<feature>`).
#'
#' @param data A list-like object containing one or more omic tables plus a
#'   `metadata` element (see `metadata_name`).
#' @param omic_names Character vector of names in `data` that contain omic
#'   tables.
#' @param subjectid Column name for subject ID.
#' @param visitid Column name for visit / timepoint.
#' @param baseline_level Value of `visitid` indicating baseline (e.g.
#'   `"Baseline"` or `0`).
#' @param metadata_name Name of the element in `data` holding sample-level
#'   metadata. Defaults to `"metadata"`.
#' @param sampleid Column name for sample ID (optional; can be `NA`).
#' @param feature_cols Optional character vector of feature (gene/protein)
#'   columns. If `NULL`, all numeric columns except id/visit/sample are used.
#'
#' @return Modified `data` list; for each omic table:
#'   - baseline rows kept as-is (no deltas),
#'   - for each non-baseline visit, new columns named
#'     `delta-<visit>-<feature>` with visit - baseline,
#'   - delta columns that are all NA are removed.
#' @export
omicDeltaByVisit <- function(data,
                              omic_names,
                              subjectid,
                              visitid,
                              baseline_level,
                              metadata_name = "metadata",
                              sampleid = NA,
                              feature_cols = NULL) {

  for (omic_name in omic_names) {
    omic <- data[[omic_name]]
    if (is.null(omic)) next

    omic <-
      omic |>
      dplyr::left_join(
        dplyr::select(
          data[[metadata_name]],
          dplyr::all_of(c(sampleid, subjectid, visitid))
        ),
        by = sampleid
      ) # Merge in the column into the omic data frame

    # Determine feature columns if not supplied
    if (is.null(feature_cols)) {
      id_cols <- c(subjectid, visitid)
      if (!is.na(sampleid)) id_cols <- c(id_cols, sampleid)

      features <- setdiff(names(omic), id_cols)
      features <- features[
        vapply(omic[features], is.numeric, logical(1))
      ]
    } else {
      features <- feature_cols
    }

    # If nothing to do for this omic, skip
    if (length(features) == 0L) {
      next
    }

    # Long format: one row per subject-visit-feature
    # Keep only the ones with baseline values
    keeps <-
      omic |>
      dplyr::filter(.data[[visitid]] == baseline_level) |>
      dplyr::pull(.data[[subjectid]])

    omic_long <-
      omic |>
      dplyr::filter(.data[[subjectid]] %in% keeps) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(features),
        names_to = "feature",
        values_to = "value"
      )

    # Join baseline and compute delta for non-baseline visits
    omic_long_delta <-
      omic_long |>
      dplyr::group_by(.data[[subjectid]]) |>
      dplyr::mutate(
        delta = dplyr::case_when(
          .data[[visitid]] == baseline_level ~ value,
          TRUE ~ value - value[which(.data[[visitid]] == baseline_level)]
        )
      )

    # Create visit-specific delta variable name: delta-<visit>-<feature>
    omic_long_delta <-
      omic_long_delta |>
      dplyr::mutate(
        delta_name = dplyr::case_when(
          .data[[visitid]] == baseline_level ~ feature,
          TRUE ~ paste0(
            "delta-",
            .data[[visitid]],
            "-",
            feature
          )
        )
      )

    # Keep original values as-is, and pivot deltas wide by delta_name
    delta_wide <-
      omic_long_delta |>
      dplyr::select(
        dplyr::all_of(subjectid),
        delta_name,
        delta
      ) |>
      tidyr::pivot_wider(
        names_from = delta_name,
        values_from = delta
      )

    # Drop delta columns that ended up all NA
    all_na_cols <- vapply(delta_wide, function(x) all(is.na(x)), logical(1))
    delta_wide <- delta_wide[, !all_na_cols, drop = FALSE]

    data[[omic_name]] <- delta_wide
  }

  return(data)
}


#' Filter omic structures and metadata
#'
#' @param data A list-like object containing `metadata` and one or more omic
#'   tables.
#' @param factor Character scalar, column name in `metadata` used for
#'   filtering.
#' @param omic_names Character vector of names in `data` that contain omic
#'   tables. Defaults to `c("olink")`.
#' @param subjectid Character scalar, column name for subject ID in
#'   `metadata` (use `NA` to skip subject-level deduplication).
#' @param sampleid_metadata Character scalar, column name for sample ID in
#'   `metadata`. Defaults to `"SampleID"`.
#' @param sampleid Character scalar, column name for sample ID in the omic
#'   tables. Defaults to `"SampleID"`.
#' @param filters Vector of values used to filter `factor`.
#'
#' @return The modified `data` list with filtered `metadata` and omic tables.
#' @export
omicFilter <- function(data,
                        factor,
                        omic_names = c("olink"),
                        subjectid = NA,
                        sampleid_metadata = "SampleID",
                        sampleid = "SampleID",
                        filters) {

  # Filter metadata, with optional subject-level distinct
  if (is.na(subjectid)) {
    data$metadata <-
      data$metadata |>
      dplyr::filter(.data[[factor]] %in% filters)
  } else {
    data$metadata <-
      data$metadata |>
      dplyr::distinct(.data[[subjectid]], .keep_all = TRUE) |>
      dplyr::filter(.data[[factor]] %in% filters)
  }

  # Compute the set of sample IDs to keep based on metadata
  keeps <- data$metadata[[sampleid_metadata]]

  # Loop over all requested omic tables and filter them
  for (omic_name in omic_names) {
    if (!is.null(data[[omic_name]])) {
      data[[omic_name]] <-
        data[[omic_name]] |>
        dplyr::filter(.data[[sampleid]] %in% keeps)
    }
  }

  # Ensure metadata is also restricted to kept sample IDs
  data$metadata <-
    data$metadata |>
    dplyr::filter(.data[[sampleid_metadata]] %in% keeps)

  return(data)
}
