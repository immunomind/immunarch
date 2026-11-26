test_that("method wrappers compile cleanly with and without dots", {
  for (fn in list(vdj_public, vdj_track)) {
    expect_output(
      invisible(compiler::cmpfun(fn, options = list(suppressAll = FALSE))),
      NA
    )
  }
})

test_that("immunarch_methods lists methods for a registered family", {
  expect_true("airr_diversity" %in% immunarch_methods())
  expect_true("chao1" %in% immunarch_methods("airr_diversity"))
  expect_error(immunarch_methods("missing_family"), "No such family name: missing_family")
  expect_error(immunarch_methods(""), "No such family name")
})

test_that("make_immunarch_features adds a method only when supplied", {
  input <- tibble::tibble(term = c("richness", "evenness"), amount = c(2, 1)) |>
    duckplyr::as_duckdb_tibble()

  out <- make_immunarch_features(input, "diversity", "term", "amount")
  expect_equal(dplyr::collect(out), tibble::tibble(
    feature = c("richness", "evenness"), value = c(2, 1), method = "diversity"
  ))

  expect_false("method" %in% names(make_immunarch_features(input)))
  existing <- dplyr::mutate(input, method = "existing")
  expect_equal(dplyr::collect(make_immunarch_features(existing)), dplyr::collect(existing))
})
