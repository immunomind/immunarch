test_that("method wrappers compile cleanly with and without dots", {
  for (fn in list(vdj_public, vdj_track)) {
    expect_output(
      invisible(compiler::cmpfun(fn, options = list(suppressAll = FALSE))),
      NA
    )
  }
})
