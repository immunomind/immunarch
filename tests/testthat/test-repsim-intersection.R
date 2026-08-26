test_that("repsim_intersection computes shared receptor counts with string repertoire ids", {
  out <- repsim_intersection(make_repsim_set_idata(), autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)

  expected <- matrix(
    c(2, 1, 2,
      1, 2, 2,
      2, 2, 3),
    nrow = 3,
    byrow = TRUE
  )
  dimnames(expected) <- list(c("R1", "R2", "R3"), c("R1", "R2", "R3"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_intersection treats receptor sets as unique per repertoire", {
  out <- repsim_intersection(
    make_repsim_set_idata(duplicated = TRUE),
    autojoin = FALSE
  )

  expected <- matrix(c(2, 2, 2, 2), nrow = 2, byrow = TRUE)
  dimnames(expected) <- list(c("R1", "R2"), c("R1", "R2"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_intersection returns a plottable matrix", {
  skip_if_not_installed("ggplot2")

  idata <- make_aggregated_test_idata()

  out <- repsim_intersection(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)
  expect_true(all(diag(out) >= 0))

  p <- vis(out)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})
