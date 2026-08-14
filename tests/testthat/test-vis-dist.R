make_vis_dist_result <- function() {
  dist_hamm(
    make_grouped_distance_test_idata(),
    by = "subject_id",
    autojoin = FALSE
  )
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
  receptor_col <- immundata::imd_schema("receptor")
  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(features = "cdr3_aa"),
    annotations = tibble::tibble(
      !!receptor_col := 1:3,
      cdr3_aa = c("AAA", "AAT", "ATT")
    ) |>
      duckplyr::as_duckdb_tibble(prudence = "stingy")
  )

  p <- dist_hamm(idata, autojoin = FALSE) |> vis()

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
