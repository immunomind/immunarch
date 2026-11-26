# Keep test-wide state consistent across test_check(), test_local(), and
# test_file(), then restore it after the complete test run.
withr::local_envvar(
  c(
    DUCKPLYR_FALLBACK_INFO = "FALSE",
    DUCKPLYR_FALLBACK_COLLECT = "0",
    DUCKPLYR_FALLBACK_AUTOUPLOAD = "0",
    DUCKPLYR_FALLBACK_VERBOSE = "FALSE"
  ),
  .local_envir = testthat::teardown_env()
)

withr::local_options(
  list(
    immundata.verbose = FALSE,
    rlib_message_verbosity = "quiet"
  ),
  .local_envir = testthat::teardown_env()
)
