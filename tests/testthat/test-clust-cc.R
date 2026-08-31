# Connected-component clustering tests.


testthat::test_that("clust_cc retains transitive components and singletons", {
  source <- tibble::tibble(
    scope = 1L,
    node = 1:4,
    from = c(1L, 2L, NA_integer_, NA_integer_),
    to = c(2L, 3L, NA_integer_, NA_integer_),
    imd_receptor_id = 1:4
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")

  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(features = "node"),
    annotations = source
  )

  edges <- source |>
    dplyr::filter(!is.na(.data$from)) |>
    dplyr::select(dplyr::all_of(c("scope", "from", "to")))

  lazy <- clust_cc(
    idata = idata,
    edges = edges,
    from = "from",
    to = "to",
    by = "scope",
    cluster_col = "component",
    size_col = "component_size"
  )

  testthat::expect_s3_class(lazy, "immunarch_res_clust_cc")

  invisible(gc())

  actual <- lazy |>
    dplyr::arrange(imd_receptor_id) |>
    dplyr::collect()

  testthat::expect_equal(actual$component, c(1L, 1L, 1L, 4L))
  testthat::expect_equal(actual$component_size, c(3L, 3L, 3L, 1L))
})


testthat::test_that("clust_cc materialization survives an RNG reset", {
  source <- tibble::tibble(
    scope = 1L,
    from = c(1L, 2L, NA_integer_),
    to = c(2L, 3L, NA_integer_),
    imd_receptor_id = 1:3
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")

  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(features = "from"),
    annotations = source
  )
  edges <- source |>
    dplyr::filter(!is.na(.data$from)) |>
    dplyr::select(dplyr::all_of(c("scope", "from", "to")))

  withr::local_preserve_seed()
  set.seed(431)
  first <- clust_cc(
    idata,
    edges,
    from = "from",
    to = "to",
    by = "scope"
  )
  invisible(gc())

  set.seed(431)
  second <- clust_cc(
    idata,
    edges,
    from = "from",
    to = "to",
    by = "scope"
  )

  testthat::expect_equal(
    dplyr::collect(dplyr::arrange(second, imd_receptor_id)),
    dplyr::collect(dplyr::arrange(first, imd_receptor_id))
  )
})
