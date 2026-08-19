test_that("omicFilter subsets metadata and omic tables consistently", {
  metadata <- data.frame(
    SampleID = paste0("S", 1:6),
    Group = c("A", "A", "B", "B", "C", "C"),
    stringsAsFactors = FALSE
  )
  olink <- data.frame(
    SampleID = paste0("S", 1:6),
    `GENE1-OID` = rnorm(6),
    check.names = FALSE
  )
  data <- list(metadata = metadata, olink = olink)

  out <- omicFilter(
    data = data,
    factor = "Group",
    omic_names = "olink",
    filters = c("A", "B")
  )

  expect_equal(nrow(out$metadata), 4)
  expect_equal(nrow(out$olink), 4)
  expect_true(all(out$metadata$Group %in% c("A", "B")))
  expect_true(all(out$olink$SampleID %in% out$metadata$SampleID))
})

test_that("omicFilter respects subject-level de-duplication", {
  metadata <- data.frame(
    SampleID = paste0("S", 1:4),
    SubjectID = c("P1", "P1", "P2", "P2"),
    Group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  data <- list(metadata = metadata)

  out <- omicFilter(
    data = data,
    factor = "Group",
    omic_names = character(0),
    subjectid = "SubjectID",
    filters = c("A", "B")
  )

  expect_equal(nrow(out$metadata), 2)
})

test_that("omicDeltaByVisit computes deltas from baseline", {
  metadata <- data.frame(
    SampleID = paste0("S", 1:4),
    SubjectID = c("P1", "P1", "P2", "P2"),
    Visit = c("Baseline", "Week4", "Baseline", "Week4"),
    stringsAsFactors = FALSE
  )
  omic <- data.frame(
    SampleID = paste0("S", 1:4),
    FeatureA = c(10, 15, 20, 18)
  )
  data <- list(metadata = metadata, omic = omic)

  out <- omicDeltaByVisit(
    data = data,
    omic_names = "omic",
    subjectid = "SubjectID",
    visitid = "Visit",
    baseline_level = "Baseline",
    sampleid = "SampleID"
  )

  result <- out$omic
  expect_true("delta-Week4-FeatureA" %in% names(result))

  p1_delta <- result$`delta-Week4-FeatureA`[result$SubjectID == "P1"]
  p2_delta <- result$`delta-Week4-FeatureA`[result$SubjectID == "P2"]
  expect_equal(p1_delta[!is.na(p1_delta)], 5)
  expect_equal(p2_delta[!is.na(p2_delta)], -2)
})
