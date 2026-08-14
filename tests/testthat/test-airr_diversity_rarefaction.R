make_rarefaction_idata <- function(counts, repertoire_id = "R1") {
  receptor_ids <- paste0("receptor_", seq_along(counts))

  make_test_repertoire_idata(
    receptors = receptor_ids,
    repertoires = rep(repertoire_id, length(counts)),
    counts = counts,
    proportions = counts / sum(counts),
    all_repertoires = repertoire_id
  )
}


run_v1_rarefaction <- function(counts, ...) {
  airr_diversity_rarefaction(
    make_rarefaction_idata(counts),
    ...,
    verbose = FALSE,
    autojoin = FALSE
  )
}


test_that("airr_diversity_rarefaction returns normalized interpolation curves", {

  idata <- make_aggregated_test_idata()
  rep_col <- immundata::imd_schema("repertoire")

  res <- airr_diversity_rarefaction(
    idata,
    step = 100,
    extrapolation = 0,
    nboot = 0,
    norm = TRUE,
    verbose = FALSE,
    autojoin = FALSE
  )

  expect_true(all(c(rep_col, "size", "q_low", "mean", "q_high", "type") %in% names(res)))
  expect_gt(nrow(res), 0)
  expect_true(all(res$type == "interpolation"))
  expect_lte(max(res$size, na.rm = TRUE), 1 + 1e-8)
  expect_lte(max(res$mean, na.rm = TRUE), 1 + 1e-8)
})


test_that("airr_diversity_rarefaction supports extrapolation and vis()", {

  idata <- make_aggregated_test_idata()
  rep_sym <- immundata::imd_schema_sym("repertoire")
  rec_sym <- immundata::imd_schema_sym("receptor")
  cnt_sym <- immundata::imd_schema_sym("count")

  max_n <- idata$annotations |>
    dplyr::select(!!rep_sym, !!rec_sym, !!cnt_sym) |>
    dplyr::distinct(!!rec_sym, !!rep_sym, .keep_all = TRUE) |>
    dplyr::summarise(n = sum(!!cnt_sym), .by = !!rep_sym) |>
    dplyr::collect() |>
    dplyr::pull(.data$n) |>
    max()

  step <- max(1, floor(max_n / 10))
  set.seed(2026)
  res <- airr_diversity_rarefaction(
    idata,
    step = step,
    extrapolation = max_n + step,
    nboot = 2,
    norm = FALSE,
    verbose = FALSE,
    autojoin = FALSE
  )

  expect_true(any(res$type == "extrapolation"))
  expect_true(all(is.finite(res$q_low)))
  expect_true(all(is.finite(res$q_high)))

  p <- vis(res)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("Chao1 and rarefaction are consistent at interpolation/extrapolation boundaries", {

  idata <- make_aggregated_test_idata()
  rep_col <- immundata::imd_schema("repertoire")
  rep_sym <- immundata::imd_schema_sym("repertoire")
  rec_sym <- immundata::imd_schema_sym("receptor")
  cnt_sym <- immundata::imd_schema_sym("count")

  rep_counts <- idata$annotations |>
    dplyr::select(!!rep_sym, !!rec_sym, !!cnt_sym) |>
    dplyr::distinct(!!rec_sym, !!rep_sym, .keep_all = TRUE)

  observed <- rep_counts |>
    dplyr::summarise(s_obs = dplyr::n(), .by = !!rep_sym) |>
    dplyr::collect()

  max_n <- rep_counts |>
    dplyr::summarise(n = sum(!!cnt_sym), .by = !!rep_sym) |>
    dplyr::collect() |>
    dplyr::pull(.data$n) |>
    max()

  step <- max(1, floor(max_n / 20))

  raref <- airr_diversity_rarefaction(
    idata,
    step = step,
    extrapolation = max_n + 5 * step,
    nboot = 0,
    norm = FALSE,
    verbose = FALSE,
    autojoin = FALSE
  )

  chao <- airr_diversity_chao1(idata, autojoin = FALSE) |>
    dplyr::select(all_of(c(rep_col, "Estimator")))

  interp_end <- raref |>
    dplyr::filter(.data$type == "interpolation") |>
    dplyr::summarise(s_interp_end = .data$mean[which.max(.data$size)], .by = all_of(rep_col))

  interp_cmp <- observed |>
    dplyr::left_join(interp_end, by = rep_col) |>
    dplyr::arrange(.data[[rep_col]])

  # At full observed sampling depth, rarefaction interpolation equals observed richness.
  expect_equal(interp_cmp$s_interp_end, interp_cmp$s_obs, tolerance = 1e-8)

  if (any(raref$type == "extrapolation")) {
    ex_cmp <- raref |>
      dplyr::filter(.data$type == "extrapolation") |>
      dplyr::left_join(observed, by = rep_col) |>
      dplyr::left_join(chao, by = rep_col)

    # Finite extrapolation is bounded by observed richness and Chao1 asymptote.
    expect_true(all(ex_cmp$mean >= ex_cmp$s_obs - 1e-8))
    expect_true(all(ex_cmp$mean <= ex_cmp$Estimator + 1e-8))
  }
})


test_that("rarefaction interpolation matches exact hypergeometric oracles", {

  adversarial <- run_v1_rarefaction(
    c(9, 1),
    step = 1,
    quantile = c(.025, .975),
    extrapolation = 0,
    nboot = 0,
    norm = FALSE
  )
  expect_equal(adversarial$mean[adversarial$size == 1], 1)
  expect_equal(adversarial$mean, c(1, seq(1.2, 2, by = .1)), tolerance = 1e-12)

  counts <- c(4, 3, 2, 1)
  expected <- c(
    1, 1.77777777777778, 2.375, 2.82857142857143, 3.17063492063492,
    3.42857142857143, 3.625, 3.77777777777778, 3.9, 4
  )
  curve <- run_v1_rarefaction(
    counts,
    step = 1,
    quantile = c(.025, .975),
    extrapolation = 0,
    nboot = 0,
    norm = FALSE
  )
  expect_equal(curve$mean, expected, tolerance = 1e-12)
  expect_equal(curve$mean[1], 1)
  expect_equal(tail(curve$mean, 1), length(counts))
  expect_true(all(diff(curve$mean) >= 0))
  expect_true(all(curve$mean <= pmin(curve$size, length(counts))))

  singleton_curve <- run_v1_rarefaction(
    rep(1, 4),
    step = 1,
    quantile = c(.025, .975),
    extrapolation = 0,
    nboot = 0,
    norm = FALSE
  )
  expect_equal(singleton_curve$mean, 1:4)
})


test_that("rarefaction interpolation matches vegan", {

  skip_if_not_installed("vegan")
  counts <- c(9, 4, 3, 2, 1, 1)
  sizes <- seq_len(sum(counts))
  curve <- run_v1_rarefaction(
    counts,
    step = 1,
    quantile = c(.025, .975),
    extrapolation = 0,
    nboot = 0,
    norm = FALSE
  )

  vegan_estimate <- as.numeric(vegan::rarefy(counts, sample = sizes))
  expect_equal(curve$mean, vegan_estimate, tolerance = 1e-12)
})


test_that("rarefaction uses finite-sample Chao extrapolation to twice N", {

  counts <- c(3, 2, 1, 1)
  n <- sum(counts)
  f0_hat <- (n - 1) / n * 2^2 / (2 * 1)
  A <- n * f0_hat / (n * f0_hat + 2)
  expected_at_2n <- length(counts) + f0_hat * (1 - A^n)

  curve <- run_v1_rarefaction(
    counts,
    step = 2,
    quantile = c(.025, .975),
    nboot = 0,
    norm = FALSE
  )

  expect_equal(max(curve$size), 2 * n)
  expect_equal(curve$mean[curve$size == 2 * n], expected_at_2n, tolerance = 1e-12)
  expect_equal(curve$type[curve$size == n], "interpolation")
  expect_true(all(diff(curve$mean) >= 0))
})


test_that("bootstrap bounds honor both requested probabilities", {

  counts <- c(4, 3, 2, 1, 1, 1)
  set.seed(1201)
  central <- run_v1_rarefaction(
    counts,
    step = 2,
    quantile = c(.025, .975),
    nboot = 100,
    norm = FALSE
  )
  set.seed(1201)
  changed_lower <- run_v1_rarefaction(
    counts,
    step = 2,
    quantile = c(.1, .975),
    nboot = 100,
    norm = FALSE
  )

  expect_equal(central$q_high, changed_lower$q_high)
  expect_true(any(abs(central$q_low - changed_lower$q_low) > 1e-8))
  expect_true(all(central$q_low <= central$q_high))
  expect_true(all(central$q_low >= 1))
  expect_true(all(central$q_high <= central$size))
  expect_true(any(
    central$q_low[central$type == "extrapolation"] <
      central$q_high[central$type == "extrapolation"]
  ))

  without_ci <- run_v1_rarefaction(
    counts,
    step = 2,
    quantile = c(.025, .975),
    extrapolation = 0,
    nboot = 0,
    norm = FALSE
  )
  expect_true(all(is.na(without_ci$q_low)))
  expect_true(all(is.na(without_ci$q_high)))
})


test_that("v1 bootstrap handles a large estimated unseen-species pool", {

  counts <- rep(1, 1000)
  set.seed(42)
  curve <- run_v1_rarefaction(
    counts,
    step = 1000,
    quantile = c(.025, .975),
    extrapolation = 0,
    nboot = 2,
    norm = FALSE
  )

  expect_equal(curve$size, c(1, 1000))
  expect_true(all(is.finite(curve$q_low)))
  expect_true(all(is.finite(curve$q_high)))
})


test_that("rarefaction validates scientific count and CI inputs", {

  expect_error(
    run_v1_rarefaction(
      c(2.5, 1), step = 1, quantile = c(.025, .975),
      extrapolation = 0, nboot = 0, norm = FALSE
    ),
    "integer clonotype counts"
  )
  expect_error(
    run_v1_rarefaction(
      c(2, 0), step = 1, quantile = c(.025, .975),
      extrapolation = 0, nboot = 0, norm = FALSE
    ),
    "finite, positive"
  )
  expect_error(
    run_v1_rarefaction(
      c(2, 1), step = 1, quantile = c(0, .975),
      extrapolation = 0, nboot = 0, norm = FALSE
    ),
    "strictly between"
  )
  expect_error(
    run_v1_rarefaction(
      c(2, 1), step = 1, quantile = c(.025, .975),
      extrapolation = 0, nboot = 1, norm = FALSE
    ),
    "0 or an integer"
  )
  expect_error(
    run_v1_rarefaction(
      c(.Machine$integer.max, 1), step = 1, quantile = c(.025, .975),
      extrapolation = 0, nboot = 2, norm = FALSE
    ),
    "Bootstrap rarefaction supports"
  )
})
