test_that("airr_diversity_rarefaction returns normalized interpolation curves", {

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))
  rep_col <- immundata::imd_schema("repertoire")

  res <- airr_diversity_rarefaction(
    idata,
    step = 100,
    extrapolation = 0,
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

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))
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
  res <- airr_diversity_rarefaction(
    idata,
    step = step,
    extrapolation = max_n + step,
    norm = FALSE,
    verbose = FALSE,
    autojoin = FALSE
  )

  expect_true(any(res$type == "extrapolation"))

  p <- vis(res)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("Chao1 and rarefaction are consistent at interpolation/extrapolation boundaries", {

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))
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
