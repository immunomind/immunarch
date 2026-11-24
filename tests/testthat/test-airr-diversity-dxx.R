make_dxx_idata <- function() {
  make_test_repertoire_idata(
    receptors = paste0("r", seq_len(7)),
    repertoires = c(rep("R1", 3), rep("R2", 4)),
    proportions = c(0.60, 0.25, 0.15, 0.40, 0.30, 0.20, 0.10)
  )
}


test_that("airr_diversity_dxx computes one coverage threshold", {
  repertoire_col <- immundata::imd_schema("repertoire")

  out <- airr_diversity_dxx(
    make_dxx_idata(),
    perc = 50,
    autojoin = FALSE
  )

  expect_equal(
    tibble::as_tibble(out),
    tibble::tibble(
      !!repertoire_col := c("R1", "R2"),
      perc = c(50, 50),
      dxx = c(1, 2)
    )
  )
})


test_that("airr_diversity_dxx computes several coverage thresholds", {
  repertoire_col <- immundata::imd_schema("repertoire")

  out <- airr_diversity_dxx(
    make_dxx_idata(),
    perc = c(20, 50, 80),
    autojoin = FALSE
  )

  expect_equal(
    tibble::as_tibble(out),
    tibble::tibble(
      !!repertoire_col := rep(c("R1", "R2"), each = 3),
      perc = rep(c(20, 50, 80), times = 2),
      dxx = c(1, 1, 2, 1, 2, 3)
    )
  )
})


test_that("airr_diversity_dxx keeps full coverage when cumulative proportions round down", {
  counts <- c(17160L, rep(1L, 440L))
  idata <- make_track_idata(
    tibble::tibble(
      sample = "R1",
      receptor = paste0("r", seq_along(counts)),
      abundance = counts
    ),
    features = "receptor"
  )

  out <- airr_diversity_dxx(idata, perc = c(50, 100), autojoin = FALSE)

  expect_equal(out$perc, c(50, 100))
  expect_equal(out$dxx, c(1, length(counts)))
})
