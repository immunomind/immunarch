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
