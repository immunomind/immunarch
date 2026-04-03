test_that("repsim_bray computes Bray-Curtis from proportions", {

  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")

  ann_tbl <- tibble::tibble(
    !!receptor_col := c("r1", "r2", "r1", "r2"),
    !!repertoire_col := c("R1", "R2", "R3", "R3"),
    !!count_col := c(10, 10, 5, 5),
    !!prop_col := c(1, 1, 0.5, 0.5),
    cdr3_aa = c("r1", "r2", "r1", "r2")
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

  out <- repsim_bray(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(unname(diag(out)), c(0, 0, 0), tolerance = 1e-12)

  expected <- matrix(
    c(0, 1, 0.5,
      1, 0, 0.5,
      0.5, 0.5, 0),
    nrow = 3,
    byrow = TRUE
  )
  dimnames(expected) <- list(c("R1", "R2", "R3"), c("R1", "R2", "R3"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_bray returns a symmetric [0,1] matrix and vis() works", {
  skip_if_not_installed("ggplot2")

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))

  out <- repsim_bray(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)
  expect_equal(unname(diag(out)), rep(0, nrow(out)), tolerance = 1e-12)

  upper <- out[upper.tri(out)]
  expect_true(all(upper >= 0 & upper <= 1))

  p <- vis(out)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})
