make_autojoin_test_method <- function() {
  repertoire_col <- immundata::imd_schema("repertoire")

  core <- function(idata) {
    idata$annotations |>
      dplyr::select(dplyr::all_of(repertoire_col)) |>
      dplyr::distinct()
  }

  im_method(core, "test", "autojoin")
}


make_autojoin_test_idata <- function() {
  annotations <- tibble::tibble(
    cdr3_aa = "AAA",
    group = "case"
  )

  path <- tempfile(fileext = ".tsv")
  readr::write_tsv(annotations, path)

  immundata::read_repertoires(
    path = path,
    schema = "cdr3_aa",
    repertoire_schema = "group",
    output_folder = tempfile("immundata-test-"),
    verbose = FALSE
  )
}


test_that("autojoin uses the canonical immunarch option and defaults to TRUE", {
  old_options <- options()[intersect(
    c("immunarch.autojoin", "immundata.autojoin"),
    names(options())
  )]
  on.exit({
    options(immunarch.autojoin = NULL, immundata.autojoin = NULL)
    options(old_options)
  }, add = TRUE)

  options(immunarch.autojoin = NULL, immundata.autojoin = NULL)
  .onLoad(NULL, "immunarch")

  expect_identical(getOption(IMMUNARCH_AUTOJOIN_OPTION), TRUE)

  method <- make_autojoin_test_method()
  idata <- make_autojoin_test_idata()

  expect_true("group" %in% names(method(idata)))

  options(immunarch.autojoin = FALSE, immundata.autojoin = TRUE)
  expect_false("group" %in% names(method(idata)))

  options(immunarch.autojoin = TRUE, immundata.autojoin = FALSE)
  expect_true("group" %in% names(method(idata)))
})


test_that("package loading preserves an explicit autojoin option", {
  old_option <- getOption(IMMUNARCH_AUTOJOIN_OPTION)
  on.exit(options(immunarch.autojoin = old_option), add = TRUE)

  options(immunarch.autojoin = FALSE)
  .onLoad(NULL, "immunarch")

  expect_identical(getOption(IMMUNARCH_AUTOJOIN_OPTION), FALSE)
})
