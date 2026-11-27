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

test_that("airr_diversity_chao1 handles a single singleton without doubletons", {
  idata <- make_track_idata(
    tibble::tibble(sample = "R1", receptor = c("r1", "r2", "r3"), abundance = c(1L, 3L, 4L)),
    features = "receptor"
  )

  out <- airr_diversity_chao1(idata, autojoin = FALSE)

  counts <- c(1, 3, 4)
  n <- sum(counts)
  d_obs <- length(counts)
  var_obs <- sum(exp(-counts) - exp(-2 * counts)) - (sum(counts * exp(-counts)))^2 / n
  p_unseen <- sum(exp(-counts)) / d_obs
  z <- stats::qnorm(0.975)

  expect_equal(out$Estimator, d_obs)
  expect_equal(out$SD, sqrt(var_obs))
  expect_equal(out$Conf.95.lo, max(d_obs, d_obs / (1 - p_unseen) - z * sqrt(var_obs) / (1 - p_unseen)))
  expect_equal(out$Conf.95.hi, d_obs / (1 - p_unseen) + z * sqrt(var_obs) / (1 - p_unseen))
})

test_that("airr_diversity_chao1 handles a single-clonotype singleton repertoire", {
  idata <- make_track_idata(
    tibble::tibble(sample = "R1", receptor = "r1", abundance = 1L),
    features = "receptor"
  )

  out <- airr_diversity_chao1(idata, autojoin = FALSE)

  expect_equal(out$Estimator, 1)
  expect_true(is.finite(out$SD))
  expect_equal(out$Conf.95.lo, 1)
  expect_true(is.finite(out$Conf.95.hi) && out$Conf.95.hi > 1)
})
