test_that("airr_diversity_chao1 handles doubletons without singletons", {
  idata <- make_track_idata(
    tibble::tibble(sample = "R1", receptor = c("r1", "r2"), abundance = c(2L, 3L)),
    features = "receptor"
  )

  out <- airr_diversity_chao1(idata, autojoin = FALSE)

  expect_equal(out$Estimator, 2)
  expect_equal(out$SD, 0.359225030743)
  expect_equal(out$Conf.95.lo, 2)
  expect_equal(out$Conf.95.hi, 2.979890269729)
})

test_that("airr_diversity_chao1 handles one singleton and no doubletons", {
  idata <- make_track_idata(
    tibble::tibble(sample = "R1", receptor = c("r1", "r2", "r3"), abundance = c(1L, 3L, 4L)),
    features = "receptor"
  )

  out <- airr_diversity_chao1(idata, autojoin = FALSE)

  expect_equal(out$Estimator, 3)
  expect_equal(out$SD, 0.504228020)
  expect_equal(out$Conf.95.lo, 3)
  expect_equal(out$Conf.95.hi, 4.66642862)
})

test_that("airr_diversity_chao1 handles a repertoire with one observation", {
  idata <- make_track_idata(
    tibble::tibble(sample = "R1", receptor = "r1", abundance = 1L),
    features = "receptor"
  )

  out <- airr_diversity_chao1(idata, autojoin = FALSE)

  expect_equal(out$Estimator, 1)
  expect_equal(out$SD, 0.311783378)
  expect_equal(out$Conf.95.lo, 1)
  expect_equal(out$Conf.95.hi, 2.54869766)
})
