test_that("repsim_jaccard computes set similarity with string repertoire ids", {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")

  ann_tbl <- tibble::tibble(
    !!receptor_col := c("r1", "r3", "r2", "r3", "r1", "r2", "r3"),
    !!repertoire_col := c("R1", "R1", "R2", "R2", "R3", "R3", "R3"),
    !!count_col := rep(1, 7),
    cdr3_aa = c("r1", "r3", "r2", "r3", "r1", "r2", "r3")
  )

  rep_tbl <- tibble::tibble(
    !!repertoire_col := c("R1", "R2", "R3"),
    Group = c("R1", "R2", "R3")
  )

  ann_tbl <- duckplyr::as_duckdb_tibble(ann_tbl)
  rep_tbl <- duckplyr::as_duckdb_tibble(rep_tbl)

  idata <- immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = ann_tbl,
    repertoires = rep_tbl
  )

  out <- repsim_jaccard(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)

  expected <- matrix(
    c(1, 1 / 3, 2 / 3,
      1 / 3, 1, 2 / 3,
      2 / 3, 2 / 3, 1),
    nrow = 3,
    byrow = TRUE
  )
  dimnames(expected) <- list(c("R1", "R2", "R3"), c("R1", "R2", "R3"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_jaccard deduplicates repeated receptor rows per repertoire", {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")

  ann_tbl <- tibble::tibble(
    !!receptor_col := c("r1", "r1", "r2", "r1", "r2", "r2"),
    !!repertoire_col := c("R1", "R1", "R1", "R2", "R2", "R2"),
    !!count_col := rep(1, 6),
    cdr3_aa = c("r1", "r1", "r2", "r1", "r2", "r2")
  )

  rep_tbl <- tibble::tibble(
    !!repertoire_col := c("R1", "R2"),
    Group = c("R1", "R2")
  )

  ann_tbl <- duckplyr::as_duckdb_tibble(ann_tbl)
  rep_tbl <- duckplyr::as_duckdb_tibble(rep_tbl)

  idata <- immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = ann_tbl,
    repertoires = rep_tbl
  )

  out <- repsim_jaccard(idata, autojoin = FALSE)

  expected <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  dimnames(expected) <- list(c("R1", "R2"), c("R1", "R2"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_jaccard returns a plottable matrix with unit diagonal", {
  skip_if_not_installed("ggplot2")

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))

  out <- repsim_jaccard(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)
  expect_equal(unname(diag(out)), rep(1, nrow(out)), tolerance = 1e-12)

  p <- vis(out)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})
