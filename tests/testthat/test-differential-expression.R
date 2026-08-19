test_that("simpleWilcox returns one row per feature with expected columns", {
  set.seed(1)
  data_i <- data.frame(
    Group = rep(c("A", "B"), each = 10),
    `GENE1-OID` = c(rnorm(10, 0, 1), rnorm(10, 2, 1)),
    `GENE2-OID` = rnorm(20),
    check.names = FALSE
  )
  protein_names <- c("GENE1-OID", "GENE2-OID")

  degs <- simpleWilcox(protein_names, data_i, "Group", c("A", "B"))

  expect_equal(nrow(degs), 2)
  expect_equal(sort(degs$GeneID), c("GENE1", "GENE2"))
  expect_true(all(degs$Pvalue >= 0 & degs$Pvalue <= 1))
})

test_that("glm_model returns the coefficient row for the predictor of interest", {
  set.seed(2)
  n <- 60
  data_i <- data.frame(
    Outcome = rbinom(n, 1, 0.5),
    Predictor = rnorm(n)
  )

  res <- glm_model(
    factor = "Outcome",
    covariates = NA,
    data_i = data_i,
    x = "Predictor",
    family = stats::binomial()
  )

  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% names(res)))
})

test_that("glm_model falls back gracefully when the predictor is unidentifiable", {
  data_i <- data.frame(
    Outcome = rep(c(0, 1), 10),
    Predictor = rep(1, 20) # constant column -> collinear / dropped
  )

  res <- suppressWarnings(
    glm_model(
      factor = "Outcome",
      covariates = NA,
      data_i = data_i,
      x = "Predictor",
      family = stats::binomial()
    )
  )

  expect_type(res, "list")
  expect_equal(res$empirical_p_value, 1)
})
