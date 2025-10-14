test_that("vis() for airr_stats_lengths builds plots", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  idata <- get_test_immundata() |> agg_repertoires("Therapy")
  res_lengths <- idata |> airr_stats_lengths()
  p1 <- vis(res_lengths, fill = "Therapy")
  expect_s3_class(p1, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p1)), NA)

  p2 <- vis(res_lengths, facet = "Therapy")
  expect_s3_class(p2, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p2)), NA)
  expect_true(inherits(p2$facet, "FacetWrap"))

  p3 <- vis(res_lengths, fill = "Therapy", facet = "imd_repertoire_id")
  expect_s3_class(p3, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p3)), NA)
  expect_true(inherits(p3$facet, "FacetWrap"))

  p4 <- vis(res_lengths, fill = "seq_len", facet = "imd_repertoire_id")
  expect_s3_class(p4, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p4)), NA)
  expect_true(inherits(p4$facet, "FacetWrap"))

  if (!all(c("Therapy", "imd_repertoire_id") %in% names(res_lengths))) skip("Required columns missing")
  p5 <- vis(res_lengths, fill = "seq_len", facet = c("Therapy", "imd_repertoire_id"))
  expect_s3_class(p5, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p5)), NA)
  expect_true(inherits(p5$facet, "FacetGrid"))
})
