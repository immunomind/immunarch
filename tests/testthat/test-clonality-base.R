test_that("clonality SQL handles generated ranks and quoted bin labels", {
  idata <- make_test_repertoire_idata(
    receptors = c("r1", "r2"),
    repertoires = c("R1", "R1"),
    proportions = c(0.8, 0.2)
  )
  receptor_col <- immundata::imd_schema("receptor")

  rank_bins <- base_clonality_rank(idata, c(1L, 2L)) |>
    dplyr::collect() |>
    dplyr::arrange(!!rlang::sym(receptor_col))
  expect_equal(rank_bins$clonal_rank_bin, c(1, 2))

  prop_bins <- base_clonality_prop(
    idata,
    stats::setNames(0.5, "O'Brien")
  ) |>
    dplyr::collect() |>
    dplyr::arrange(!!rlang::sym(receptor_col))
  expect_equal(prop_bins$clonal_prop_bin, c("O'Brien", "Ultra-rare"))
})
