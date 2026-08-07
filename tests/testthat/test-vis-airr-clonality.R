test_that("vis() plots clonality lines on proportion or count scales", {
  skip_if_not_installed("ggplot2")

  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  proportion_col <- immundata::imd_schema("proportion")

  line_data <- tibble::tibble(
    !!repertoire_col := rep(c("R1", "R2"), each = 3),
    index = rep(1:3, times = 2),
    !!count_col := c(100, 50, 10, 20, 10, 2),
    !!proportion_col := c(0.5, 0.25, 0.05, 0.5, 0.25, 0.05),
    Group = rep(c("A", "B"), each = 3)
  ) |>
    im_as_result("airr_clonality", "line")

  proportion_plot <- vis(line_data)
  count_plot <- vis(line_data, yval = count_col)
  faceted_plot <- vis(line_data, facet = "Group", log = FALSE)

  expect_s3_class(proportion_plot, "ggplot")
  expect_s3_class(count_plot, "ggplot")
  expect_error(ggplot2::ggplot_build(proportion_plot), NA)
  expect_error(ggplot2::ggplot_build(count_plot), NA)
  expect_equal(proportion_plot$labels$y, "Receptor proportion")
  expect_equal(count_plot$labels$y, "Receptor count")
  expect_equal(
    ggplot2::layer_data(proportion_plot)$y,
    log10(line_data[[proportion_col]])
  )
  expect_equal(
    ggplot2::layer_data(count_plot)$y,
    log10(line_data[[count_col]])
  )
  expect_true(inherits(faceted_plot$facet, "FacetWrap"))
  expect_error(vis(line_data, yval = "occupied_prop"))
})


test_that("vis() makes stacked occupied-space plots for clonality bins", {
  skip_if_not_installed("ggplot2")

  repertoire_col <- immundata::imd_schema("repertoire")

  rank_data <- tibble::tibble(
    !!repertoire_col := rep(c("R1", "R2"), each = 3),
    clonal_rank_bin = rep(c(10, 100, NA), times = 2),
    occupied_prop = c(0.4, 0.3, 0.3, 0.2, 0.5, 0.3)
  ) |>
    im_as_result("airr_clonality", "rank")

  prop_data <- tibble::tibble(
    !!repertoire_col := rep(c("R1", "R2"), each = 2),
    clonal_prop_bin = rep(c("Large", "Rare"), times = 2),
    occupied_prop = c(0.7, 0.3, 0.2, 0.8)
  ) |>
    im_as_result("airr_clonality", "prop")

  rank_plot <- vis(rank_data)
  prop_plot <- vis(prop_data)

  expect_s3_class(rank_plot, "ggplot")
  expect_s3_class(prop_plot, "ggplot")
  expect_error(ggplot2::ggplot_build(rank_plot), NA)
  expect_error(ggplot2::ggplot_build(prop_plot), NA)
  expect_equal(rank_plot$labels$title, "Clonal space by receptor rank")
  expect_equal(prop_plot$labels$title, "Clonal space by receptor proportion")
  expect_equal(rank_plot$labels$fill, "Rank bin")
  expect_equal(prop_plot$labels$fill, "Proportion bin")

  rank_built <- ggplot2::ggplot_build(rank_plot)$data[[1]]
  prop_built <- ggplot2::ggplot_build(prop_plot)$data[[1]]
  expect_equal(as.vector(tapply(rank_built$y, rank_built$x, max)), c(0.7, 0.7))
  expect_equal(as.vector(tapply(prop_built$y, prop_built$x, max)), c(1, 1))
})
