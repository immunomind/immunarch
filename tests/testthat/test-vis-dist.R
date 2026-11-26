make_vis_dist_result <- function() {
  dist <- rep(c(1, 2, 3, 1, 2, 1), 2)
  edge_data <- tibble::tibble(
    subject_id = rep(c("P1", "P2"), each = 6L),
    imd_receptor_id_1 = c(1L, 1L, 1L, 2L, 2L, 3L, 5L, 5L, 5L, 6L, 6L, 7L),
    imd_receptor_id_2 = c(2L, 3L, 4L, 3L, 4L, 4L, 6L, 7L, 8L, 7L, 8L, 8L),
    seq_len = 4L,
    dist = dist,
    norm_dist = dist / 4,
    sim = 1 - norm_dist
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")
  im_as_result(edge_data, family = "dist", name = "hamm")
}


test_that("distance visualization defaults to normalized nearest neighbors", {
  p <- vis(make_vis_dist_result())

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 8L)
  expect_true(all(p$data$.idist_value == 1 / 4))
  expect_equal(names(p$facet$params$facets), "subject_id")
  expect_equal(p$labels$x, "norm_dist")
  expect_equal(p$labels$title, "Nearest-neighbor distance distribution")
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("distance visualization can plot all distances and select a column", {
  p <- vis(
    make_vis_dist_result(),
    mode = "all",
    xval = "sim",
    facet = NULL
  )

  expect_equal(nrow(p$data), 12L)
  expect_equal(
    sort(p$data$.idist_value),
    c(rep(1 / 4, 2), rep(1 / 2, 4), rep(3 / 4, 6))
  )
  expect_s3_class(p$facet, "FacetNull")
  expect_equal(p$labels$x, "sim")
  expect_equal(p$labels$title, "Pairwise distance distribution")
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("nearest-neighbor visualization respects similarity direction", {
  p <- vis(make_vis_dist_result(), xval = "sim")

  expect_equal(nrow(p$data), 8L)
  expect_equal(p$data$.idist_value, rep(3 / 4, 8L))
})


test_that("distance visualization pools data when no subject column exists", {
  edge_data <- tibble::tibble(
    imd_receptor_id_1 = c(1L, 1L, 2L),
    imd_receptor_id_2 = c(2L, 3L, 3L),
    seq_len = 3L,
    dist = c(1, 2, 1),
    norm_dist = c(1 / 3, 2 / 3, 1 / 3),
    sim = c(2 / 3, 1 / 3, 2 / 3)
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")
  edges <- im_as_result(edge_data, family = "dist", name = "hamm")
  p <- vis(edges)

  expect_s3_class(p$facet, "FacetNull")
  expect_equal(nrow(p$data), 3L)
})


test_that("distance visualization validates its common edge interface", {
  expect_error(
    vis(make_vis_dist_result(), xval = "missing"),
    "should be one of"
  )
  expect_error(
    vis(make_vis_dist_result(), facet = "missing"),
    "missing"
  )
  expect_error(
    vis(make_vis_dist_result(), binwidth = 0),
    "greater than 0"
  )
})


test_that("distance visualization is registered for the family", {
  expect_true(
    is.function(utils::getS3method(
      "vis",
      im_result_class("dist"),
      optional = TRUE
    ))
  )
})
