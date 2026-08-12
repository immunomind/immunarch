make_chao_jaccard_idata <- function(
  receptor_counts,
  repertoire_ids = unique(receptor_counts$repertoire)
) {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")

  ann_tbl <- receptor_counts |>
    dplyr::transmute(
      !!receptor_col := .data$receptor,
      !!repertoire_col := .data$repertoire,
      !!count_col := .data$count,
      cdr3_aa = .data$receptor
    )
  rep_tbl <- tibble::tibble(
    !!repertoire_col := repertoire_ids,
    Group = repertoire_ids
  )

  idata <- immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = duckplyr::as_duckdb_tibble(ann_tbl, prudence = "stingy"),
    repertoires = duckplyr::as_duckdb_tibble(rep_tbl, prudence = "stingy")
  )
  idata$schema_repertoire <- "Group"
  idata
}


test_that("repsim_chao_jaccard matches hard-coded vegan similarities", {
  receptor_counts <- tibble::tribble(
    ~receptor, ~repertoire, ~count,
    "r1", "A", 6,
    "r1", "A", 4,
    "r2", "A", 4,
    "r3", "A", 2,
    "r4", "A", 1,
    "r5", "A", 1,
    "r14", "A", 7,
    "r11", "A", 0,
    "r1", "B", 3,
    "r2", "B", 1,
    "r3", "B", 2,
    "r4", "B", 1,
    "r6", "B", 3,
    "r6", "B", 5,
    "r7", "B", 1,
    "r8", "B", 0,
    "r1", "C", 1,
    "r2", "C", 2,
    "r8", "C", 6,
    "r9", "C", 2,
    "r10", "C", 1,
    "r11", "D", 5,
    "r12", "D", 2,
    "r13", "D", 1
  )

  observed <- repsim_chao_jaccard(
    make_chao_jaccard_idata(receptor_counts),
    autojoin = FALSE
  )

  # Generated once with vegan 2.7-5:
  # 1 - as.matrix(vegan::vegdist(repertoire_by_receptor, method = "chao"))
  expected <- matrix(
    c(
      1, 0.43634297395619703, 0.23013415892672862, 0,
      0.43634297395619703, 1, 0.19903019616486661, 0,
      0.23013415892672862, 0.19903019616486661, 1, 0,
      0, 0, 0, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
  )

  expect_equal(unclass(observed), expected, tolerance = 1e-12)
})


test_that("repsim_chao_jaccard handles pairs without shared doubletons", {
  receptor_counts <- tibble::tribble(
    ~receptor, ~repertoire, ~count,
    "r1", "E", 5,
    "r2", "E", 1,
    "r3", "E", 3,
    "r1", "F", 1,
    "r2", "F", 1,
    "r4", "F", 4,
    "r1", "G", 2,
    "r2", "G", 3,
    "r4", "G", 1
  )

  observed <- repsim_chao_jaccard(
    make_chao_jaccard_idata(receptor_counts),
    autojoin = FALSE
  )

  # Hard-coded output from vegan::vegdist(..., method = "chao"), converted
  # from dissimilarity to similarity.
  expected <- matrix(
    c(
      1, 0.40740740740740744, 0.66666666666666674,
      0.40740740740740744, 1, 1,
      0.66666666666666674, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("E", "F", "G"), c("E", "F", "G"))
  )

  expect_equal(unclass(observed), expected, tolerance = 1e-12)
})


test_that("repsim_chao_jaccard returns a plottable similarity matrix", {
  skip_if_not_installed("ggplot2")

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))
  observed <- repsim_chao_jaccard(idata, autojoin = FALSE)

  expect_true(is.matrix(observed))
  expect_equal(observed, t(observed), tolerance = 1e-12)
  expect_equal(unname(diag(observed)), rep(1, nrow(observed)), tolerance = 1e-12)
  upper <- observed[upper.tri(observed)]
  expect_true(all(upper >= 0 & upper <= 1))

  plot <- vis(observed)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot))
})
