test_that("vis() for airr_desc_chains plots chain counts by default", {
  skip_if_not_installed("ggplot2")

  idata <- make_aggregated_test_idata("Therapy")
  res_chains <- idata |> airr_desc_chains()
  p <- vis(res_chains)

  expect_s3_class(p, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p)), NA)
  expect_equal(p$labels$title, "No. chains per sample")
  expect_equal(
    sort(ggplot2::layer_data(p)$y),
    sort(res_chains$n_chains)
  )
})
