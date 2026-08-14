test_that("repsim_jaccard computes set similarity with string repertoire ids", {
  out <- repsim_jaccard(make_repsim_set_idata(), autojoin = FALSE)

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
  out <- repsim_jaccard(
    make_repsim_set_idata(duplicated = TRUE),
    autojoin = FALSE
  )

  expected <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  dimnames(expected) <- list(c("R1", "R2"), c("R1", "R2"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_jaccard returns a plottable matrix with unit diagonal", {
  skip_if_not_installed("ggplot2")

  idata <- make_aggregated_test_idata()

  out <- repsim_jaccard(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)
  expect_equal(unname(diag(out)), rep(1, nrow(out)), tolerance = 1e-12)

  p <- vis(out)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})
