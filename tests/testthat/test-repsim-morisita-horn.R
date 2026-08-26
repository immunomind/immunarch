test_that("repsim_morisita_horn matches hard-coded vegan similarities", {
  receptor_counts <- tibble::tribble(
    ~receptor, ~repertoire, ~count,
    "r1", "A", 10,
    "r2", "A", 4,
    "r3", "A", 2,
    "r4", "A", 1,
    "r5", "A", 1,
    "r14", "A", 7,
    "r1", "B", 3,
    "r2", "B", 1,
    "r3", "B", 2,
    "r4", "B", 1,
    "r6", "B", 8,
    "r7", "B", 1,
    "r1", "C", 1,
    "r2", "C", 2,
    "r8", "C", 6,
    "r9", "C", 2,
    "r10", "C", 1,
    "r11", "D", 5,
    "r12", "D", 2,
    "r13", "D", 1
  )

  observed <- repsim_morisita_horn(
    make_repsim_idata(receptor_counts, normalize = TRUE),
    autojoin = FALSE
  )

  # Generated once with vegan 2.7-5:
  # 1 - as.matrix(vegan::vegdist(repertoire_by_receptor, method = "horn"))
  expected <- matrix(
    c(
      1, 0.33270772905647505, 0.20234571139506130, 0,
      0.33270772905647505, 1, 0.08241758241758246, 0,
      0.20234571139506130, 0.08241758241758246, 1, 0,
      0, 0, 0, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
  )

  expect_equal(unclass(observed), expected, tolerance = 1e-12)
})


test_that("repsim_morisita_horn is invariant to repertoire scaling", {
  receptor_counts <- tibble::tribble(
    ~receptor, ~repertoire, ~count,
    "r1", "A", 3,
    "r2", "A", 1,
    "r1", "B", 5,
    "r2", "B", 5
  )
  scaled_counts <- receptor_counts |>
    dplyr::mutate(
      count = .data$count * dplyr::if_else(.data$repertoire == "A", 4, 10)
    )

  observed <- repsim_morisita_horn(
    make_repsim_idata(
      dplyr::mutate(receptor_counts, proportion = as.double(.data$count))
    ),
    autojoin = FALSE
  )
  observed_scaled <- repsim_morisita_horn(
    make_repsim_idata(
      dplyr::mutate(scaled_counts, proportion = as.double(.data$count))
    ),
    autojoin = FALSE
  )

  expect_equal(observed, observed_scaled, tolerance = 1e-12)
  expect_equal(observed["A", "B"], 8 / 9, tolerance = 1e-12)
})


test_that("repsim_morisita_horn handles empty repertoires", {
  receptor_counts <- tibble::tribble(
    ~receptor, ~repertoire, ~count,
    "r1", "A", 1,
    "r2", "B", 1
  )
  observed <- repsim_morisita_horn(
    make_repsim_idata(
      receptor_counts,
      all_repertoires = c("A", "B", "empty"),
      normalize = TRUE
    ),
    autojoin = FALSE
  )

  expect_equal(observed["A", "B"], 0)
  expect_true(is.na(observed["A", "empty"]))
  expect_true(is.na(observed["B", "empty"]))
  expect_equal(unname(diag(observed)), c(1, 1, 1))
})


test_that("repsim_morisita_horn returns a plottable similarity matrix", {
  skip_if_not_installed("ggplot2")

  idata <- make_aggregated_test_idata()
  observed <- repsim_morisita_horn(idata, autojoin = FALSE)

  expect_true(is.matrix(observed))
  expect_equal(observed, t(observed), tolerance = 1e-12)
  expect_equal(unname(diag(observed)), rep(1, nrow(observed)), tolerance = 1e-12)
  expect_true(all(observed[upper.tri(observed)] >= 0 & observed[upper.tri(observed)] <= 1))

  plot <- vis(observed)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot))
})
