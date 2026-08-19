#' Wilcoxon rank-sum test across many features
#'
#' Runs a two-group Wilcoxon test independently for every feature in
#' `protein_names`, plus the difference in medians between the two groups.
#'
#' @param protein_names Character vector of feature/column names to test.
#'   Expected to contain a suffix separated by `"-"`, e.g. `"GENE-OID"`; only
#'   the part before the first `"-"` is kept as `GeneID`.
#' @param data_i Data frame containing `factor` and all of `protein_names`.
#' @param factor Character scalar, name of the grouping column in `data_i`.
#' @param comparison Length-2 character vector giving the two levels of
#'   `factor` to compare.
#'
#' @return A data frame with columns `GeneID`, `Log2FC` (difference in
#'   medians) and `Pvalue`.
#' @export
simpleWilcox <- function(protein_names, data_i, factor, comparison) {
  degs <- data.frame(
    GeneID = stringr::str_split(
      protein_names,
      pattern = "-",
      simplify = TRUE
    )[, 1],
    Log2FC = sapply(
      protein_names,
      FUN = function(x) {
        median(na.omit(data_i[data_i[, factor] == comparison[1], x])) -
          median(na.omit(data_i[data_i[, factor] == comparison[2], x]))
      }
    ),
    Pvalue = sapply(
      protein_names,
      FUN = function(x) {
        stats::wilcox.test(
          stats::as.formula(paste0("`", x, "` ~ `", factor, "`")),
          data = data_i
        )$p.value
      }
    )
  )
  return(degs)
}

#' Fit a single-predictor GLM, with optional permutation p-value
#'
#' @param factor Character scalar, name of the response column.
#' @param covariates Character vector of covariate column names, or `NA` to
#'   fit an unadjusted model.
#' @param data_i Data frame containing `factor`, `covariates`, and `x`.
#' @param x Character scalar, name of the predictor column.
#' @param family A family object as used by [stats::glm()], e.g.
#'   `binomial()` or `gaussian()`.
#' @param empirical Logical; if `TRUE`, compute an empirical p-value via
#'   label permutation instead of the asymptotic p-value.
#' @param n_perm Number of permutations used when `empirical = TRUE`.
#' @param seed Random seed for permutations.
#'
#' @return Either a named numeric vector of coefficient statistics for `x`
#'   (when `empirical = FALSE`), or a list with `coefficient` and
#'   `empirical_p_value` (when `empirical = TRUE`).
#' @export
glm_model <- function(factor,
                       covariates,
                       data_i,
                       x,
                       family = stats::binomial(),
                       empirical = FALSE,
                       n_perm = 1000,
                       seed = 123) {

  set.seed(seed)

  response <- paste0("`", factor, "`")
  predictor <- paste0("`", x, "`")

  if (length(na.omit(covariates)) == 0) {
    formula <- stats::as.formula(paste(response, "~", predictor))
  } else {
    cov_str <- paste0("`", covariates, "`", collapse = " + ")
    formula <- stats::as.formula(paste(response, "~", predictor, "+", cov_str))
  }

  fit_obs <- stats::glm(formula = formula, data = data_i, family = family)
  coefs_obs <- summary(fit_obs)$coefficients

  if (!x %in% rownames(coefs_obs) && !paste0("`", x, "`") %in% rownames(coefs_obs)) {
    warning("Predictor name not found in model coefficients or complete colinearity exist.")
    return(list(
      coefficient = data.frame(
        `Estimate` = 0,
        `Std. Error` = 1,
        `z value` = 0,
        `Pr(>|z|)` = 1
      ),
      empirical_p_value = 1
    ))
  } else {
    idx <- if (x %in% rownames(coefs_obs)) x else paste0("`", x, "`")

    if (!empirical) {
      return(coefs_obs[idx, ])
    }

    obs_z <- abs(coefs_obs[idx, "z value"])
    perm_z <- numeric(n_perm)

    response_name <- factor

    for (b in seq_len(n_perm)) {
      dat_perm <- data_i
      dat_perm[[response_name]] <- sample(dat_perm[[response_name]])

      fit_perm <- stats::glm(formula = formula, data = dat_perm, family = family)
      coefs_perm <- summary(fit_perm)$coefficients

      if (idx %in% rownames(coefs_perm)) {
        perm_z[b] <- abs(coefs_perm[idx, "z value"])
      } else {
        perm_z[b] <- NA
      }
    }

    perm_z <- perm_z[!is.na(perm_z)]
    emp_p <- (sum(perm_z >= obs_z) + 1) / (length(perm_z) + 1)

    out <- coefs_obs[idx, ]

    return(list(
      coefficient = out,
      empirical_p_value = emp_p
    ))
  }
}

#' Fit a single-predictor GLM or GLMM, optionally in parallel with a
#' permutation p-value
#'
#' Same purpose as [glm_model()], but supports an optional random-effect
#' `level_id` (fit via `lme4`/`lmerTest`) and runs permutations in parallel
#' via [parallel::mclapply()].
#'
#' @inheritParams glm_model
#' @param level_id Optional character scalar, name of a grouping column to
#'   add as a `(1 | level_id)` random intercept.
#' @param n_cores Number of cores to use for permutations.
#'
#' @return A list with `coefficient` (named numeric vector) and
#'   `empirical_p_value`.
#' @export
glm_model_mcore <- function(
  factor,
  covariates,
  data_i,
  x,
  family = stats::binomial(),
  empirical = FALSE,
  n_perm = 1000,
  seed = 123,
  level_id = NULL,
  n_cores = parallel::detectCores() - 3
) {
  set.seed(seed)

  response <- paste0("`", factor, "`")
  predictor <- paste0("`", x, "`")

  if (length(na.omit(covariates)) == 0) {
    formula_str <- paste(response, "~", predictor)
  } else {
    cov_str <- paste0("`", covariates, "`", collapse = " + ")
    formula_str <- paste(response, "~", predictor, "+", cov_str)
  }

  if (!is.null(level_id)) {
    formula <- stats::as.formula(paste0(formula_str, " + (1 | `", level_id, "`)"))
    if (family$family == "gaussian") {
      fit_obs <- lmerTest::lmer(
        formula = formula,
        data = data_i
      )
      coefs_obs <- summary(fit_obs)$coefficients
    } else if (family$family == "binomial") {
      ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
      fit_obs <- lme4::glmer(
        formula = formula,
        data = data_i,
        family = family,
        control = ctrl
      )
      coefs_obs <- summary(fit_obs)$coefficients
    }
  } else {
    formula <- stats::as.formula(formula_str)
    fit_obs <- stats::glm(formula = formula, data = data_i, family = family)
    coefs_obs <- summary(fit_obs)$coefficients
  }

  if (
    !x %in% rownames(coefs_obs) && !paste0("`", x, "`") %in% rownames(coefs_obs)
  ) {
    warning(
      "Predictor name not found in model coefficients or complete colinearity exists."
    )
    coefs <- c(0, 1, 0, 1)
    names(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    return(list(
      coefficient = coefs,
      empirical_p_value = 1
    ))
  }

  idx <- if (x %in% rownames(coefs_obs)) x else paste0("`", x, "`")

  if (!empirical) {
    return(list(
      coefficient = coefs_obs[idx, ],
      empirical_p_value = 1
    ))
  }

  obs_z <- abs(coefs_obs[idx, "z value"])

  # Parallel over permutations
  response_name <- factor

  perm_z_list <- parallel::mclapply(
    X = seq_len(n_perm),
    mc.cores = n_cores,
    FUN = function(b) {
      dat_perm <- data_i
      dat_perm[[response_name]] <- sample(dat_perm[[response_name]])

      if (!is.null(level_id)) {
        if (family$family == "gaussian") {
          fit_perm <- lmerTest::lmer(
            formula = formula,
            data = dat_perm
          )
          coefs_perm <- summary(fit_perm)$coefficients
        } else if (family$family == "binomial") {
          fit_perm <- lme4::glmer(
            formula = formula,
            data = dat_perm,
            family = family
          )
          coefs_perm <- summary(fit_perm)$coefficients
        }
      } else {
        fit_perm <- stats::glm(formula = formula, data = dat_perm, family = family)
        coefs_perm <- summary(fit_perm)$coefficients
      }

      if (idx %in% rownames(coefs_perm)) {
        return(abs(coefs_perm[idx, "z value"]))
      } else {
        return(NA_real_)
      }
    }
  )

  perm_z <- unlist(perm_z_list, use.names = FALSE)
  perm_z <- perm_z[!is.na(perm_z)]

  if (!length(perm_z)) {
    emp_p <- 1
    warning(
      "No valid permutation statistics for ",
      x,
      "; empirical p-value set to 1."
    )
  } else {
    emp_p <- (sum(perm_z >= obs_z) + 1) / (length(perm_z) + 1)
  }

  out <- coefs_obs[idx, ]

  return(list(
    coefficient = out,
    empirical_p_value = emp_p
  ))
}

#' Logistic or linear regression across many features
#'
#' Iterates [glm_model_mcore()] over every entry in `protein_names` and
#' assembles a tidy differential-expression-style table.
#'
#' @param protein_names Character vector of feature/column names to test.
#' @param data_i Data frame containing `factor`, `covariates`, and all of
#'   `protein_names`.
#' @param factor Character scalar, name of the response column.
#' @param level_id Optional character scalar, random-effect grouping column.
#' @param comparison Length-2 character vector; for `type = "logistic"`, the
#'   two levels of `factor` used to define the reference/comparison coding.
#' @param oid_suffix Character scalar, suffix used to split `protein_names`
#'   into a `GeneID`.
#' @param covariates Character vector of covariate column names, or `NA`.
#' @param family A family object as used by [stats::glm()].
#' @param type One of `"logistic"` or `"regression"`.
#' @param empirical Logical; use permutation p-values.
#' @param n_perm Number of permutations when `empirical = TRUE`.
#' @param n_cores Number of cores for permutations.
#'
#' @return A data frame with columns `GeneID`, `Log2FC`, `Pvalue`.
#' @export
LogReg <- function(
  protein_names,
  data_i,
  factor,
  level_id = NULL,
  comparison,
  oid_suffix,
  covariates,
  family,
  type = c("logistic", "regression"),
  empirical = FALSE,
  n_perm = 1000,
  n_cores = round((parallel::detectCores() - 6) / 2)
) {
  type <- match.arg(type)

  if (type == "logistic") {
    data_i[, factor] <-
      data_i[, factor] |>
      (\(x) factor(x = x, levels = comparison))() |>
      as.numeric() -
      1
  }

  n_prot <- length(protein_names)
  cat(
    "Running logistic regression for",
    n_prot,
    "proteins on",
    n_cores,
    "cores\n"
  )
  pb <- utils::txtProgressBar(min = 0, max = n_prot, style = 3)

  # 1) iterate over indices, not names
  run_one <- function(i) {
    x <- protein_names[i]
    res <- glm_model_mcore(
      factor = factor,
      covariates = covariates,
      data_i = data_i,
      level_id = level_id,
      x = x,
      empirical = empirical,
      n_perm = n_perm,
      n_cores = n_cores,
      family = family
    )
    utils::setTxtProgressBar(pb, i)
    res
  }

  fits <- lapply(seq_along(protein_names), run_one)
  close(pb)
  names(fits) <- protein_names

  coef_mat <- dplyr::bind_rows(lapply(fits, `[[`, "coefficient"))

  log2fc <- coef_mat[, "Estimate"] / log(2)

  if (!empirical) {
    pval <- coef_mat[, grep(pattern = "Pr", x = colnames(coef_mat))[1]]
  } else {
    pval <- sapply(
      X = fits,
      FUN = function(x) x[["empirical_p_value"]],
      USE.NAMES = TRUE
    )
  }

  degs <- data.frame(
    GeneID = stringr::str_split_i(protein_names, pattern = oid_suffix, i = 1),
    Log2FC = log2fc,
    Pvalue = pval,
    row.names = protein_names,
    check.names = FALSE
  )
  colnames(degs) <- c("GeneID", "Log2FC", "Pvalue")

  return(degs)
}

#' Gaussian regression across many features
#'
#' @param protein_names Character vector of feature/column names to test.
#' @param data_i Data frame containing `factor`, `covariates`, and all of
#'   `protein_names`.
#' @param factor Character scalar, name of the (numeric) response column.
#' @param comparison Unused; kept for interface parity with [simpleWilcox()].
#' @param covariates Character vector of covariate column names, or `NA`.
#'
#' @return A data frame with columns `GeneID`, `Log2FC` (slope estimate) and
#'   `Pvalue`.
#' @export
GaussianReg <- function(protein_names, data_i, factor, comparison, covariates) {
  # turn the factor column into numerical values
  data_i[, factor] <-
    data_i[, factor] |>
    as.numeric()

  degs <- data.frame(
    GeneID = stringr::str_split(
      protein_names,
      pattern = "-",
      simplify = TRUE
    )[, 1],
    Log2FC = sapply(
      protein_names,
      FUN = function(x) {
        glm_model(factor, covariates, data_i, x, family = "gaussian")["Estimate"]
      }
    ),
    Pvalue = sapply(
      protein_names,
      FUN = function(x) {
        glm_model(factor, covariates, data_i, x, family = "gaussian")["Pr(>|t|)"]
      }
    )
  )
  return(degs)
}

#' Partial correlation via rank-based regression for one feature
#'
#' Uses [Rfit::rfit()] to fit a rank-based regression and returns the
#' coefficient row for `x`.
#'
#' @param factor Character scalar, name of the response column.
#' @param covariates Character vector of covariate column names, or `NA`.
#' @param data_i Data frame containing `factor`, `covariates`, and `x`.
#' @param x Character scalar, name of the predictor column.
#'
#' @return A named numeric vector with the coefficient row for `x`.
#' @export
rg_model <- function(factor, covariates, data_i, x) {
  if (length(na.omit(covariates)) == 0) {
    formula <- stats::as.formula(paste0("`", factor, "`~`", x, "`"))
    dat <-
      data_i[, c(factor, x)] |>
      na.omit()
  } else {
    formula <- stats::as.formula(
      paste0("`", factor, "`~`", x, "`+`", paste0(covariates, collapse = "`+`"), "`")
    )
    dat <-
      data_i[, c(factor, covariates, x)] |>
      na.omit()
  }

  model <- Rfit::rfit(formula = formula, data = dat)
  coefs <- summary(model)$coefficients
  i <- which(rownames(coefs) == paste0("`", x, "`"))
  return(coefs[i, ])
}

#' Rank-based regression across many features
#'
#' Iterates [rg_model()] over every entry in `protein_names`. All covariates
#' must already be numeric.
#'
#' @param protein_names Character vector of feature/column names to test.
#' @param data_i Data frame containing `factor`, `covariates`, and all of
#'   `protein_names`.
#' @param factor Character scalar, name of the (numeric) response column.
#' @param covariates Character vector of covariate column names, or `NA`;
#'   must already be numeric.
#'
#' @return A data frame with columns `GeneID`, `Log2FC` (rank estimate) and
#'   `Pvalue`.
#' @export
RankReg <- function(protein_names, data_i, factor, covariates) {
  # turn the factor column into numerical values
  data_i[, factor] <-
    data_i[, factor] |>
    as.numeric()

  # turn data_i into data.frame
  data_i <-
    data.frame(data_i, check.names = FALSE)

  degs <- data.frame(
    GeneID = stringr::str_split(
      protein_names,
      pattern = "-",
      simplify = TRUE
    )[, 1],
    Log2FC = sapply(
      protein_names,
      FUN = function(x) {
        rg_model(factor, covariates, data_i, x)["Estimate"]
      }
    ),
    Pvalue = sapply(
      protein_names,
      FUN = function(x) {
        rg_model(factor, covariates, data_i, x)["p.value"]
      }
    )
  )
  return(degs)
}
