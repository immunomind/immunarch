test_that("repsim_bray computes Bray-Curtis from proportions", {

  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")

  ann_tbl <- tibble::tibble(
    !!receptor_col := c("r1", "r2", "r1", "r2"),
    !!repertoire_col := c("R1", "R2", "R3", "R3"),
    !!count_col := c(10, 10, 5, 5),
    !!prop_col := c(1, 1, 0.5, 0.5),
    cdr3_aa = c("r1", "r2", "r1", "r2")
  )

  rep_tbl <- tibble::tibble(
    !!repertoire_col := c("R1", "R2", "R3"),
    Group = c("R1", "R2", "R3")
  )
  ann_tbl <- duckplyr::as_duckdb_tibble(ann_tbl)
  rep_tbl <- duckplyr::as_duckdb_tibble(rep_tbl)

  idata <- immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = ann_tbl,
    repertoires = rep_tbl
  )

  out <- repsim_bray(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(unname(diag(out)), c(0, 0, 0), tolerance = 1e-12)

  expected <- matrix(
    c(0, 1, 0.5,
      1, 0, 0.5,
      0.5, 0.5, 0),
    nrow = 3,
    byrow = TRUE
  )
  dimnames(expected) <- list(c("R1", "R2", "R3"), c("R1", "R2", "R3"))

  expect_equal(out, expected, tolerance = 1e-12)
})


test_that("repsim_bray transforms abundances", {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")

  ann_tbl <- tibble::tibble(
    !!receptor_col := c("r1", "r2", "r1", "r2"),
    !!repertoire_col := c("R1", "R1", "R2", "R2"),
    !!count_col := c(9, 1, 5, 5),
    !!prop_col := c(0.9, 0.1, 0.5, 0.5),
    cdr3_aa = c("r1", "r2", "r1", "r2")
  )
  rep_tbl <- tibble::tibble(
    !!repertoire_col := c("R1", "R2"),
    Group = c("R1", "R2")
  )

  idata <- immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = duckplyr::as_duckdb_tibble(ann_tbl),
    repertoires = duckplyr::as_duckdb_tibble(rep_tbl)
  )
  idata$schema_repertoire <- "Group"

  reference_bray <- function(transform) {
    abundance <- matrix(
      c(0.9, 0.1, 0.5, 0.5),
      nrow = 2,
      dimnames = list(c("r1", "r2"), c("R1", "R2"))
    )
    if (transform == "log1p") {
      abundance <- log1p(abundance)
    }
    sum(abs(abundance[, 1] - abundance[, 2])) / sum(abundance)
  }

  default <- repsim_bray(idata, autojoin = FALSE)
  log_transformed <- repsim_bray(
    idata,
    transform = "log1p",
    autojoin = FALSE
  )

  expect_equal(default[1, 2], reference_bray("none"), tolerance = 1e-12)
  expect_equal(
    log_transformed[1, 2],
    reference_bray("log1p"),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(log_transformed, default)))
})


test_that("repsim_bray returns a symmetric [0,1] matrix and vis() works", {
  skip_if_not_installed("ggplot2")

  idata <- get_test_immundata() |> agg_repertoires(c("Response", "Therapy"))

  out <- repsim_bray(idata, autojoin = FALSE)

  expect_true(is.matrix(out))
  expect_equal(out, t(out), tolerance = 1e-12)
  expect_equal(unname(diag(out)), rep(0, nrow(out)), tolerance = 1e-12)

  upper <- out[upper.tri(out)]
  expect_true(all(upper >= 0 & upper <= 1))

  p <- vis(out)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("repsim_bray matches an independent weighted reference", {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")

  for (seed in 1:5) {
    set.seed(seed)
    repertoire_ids <- sprintf("R%02d", 1:12)
    receptor_ids <- sprintf("r%04d", 1:500)
    ann_tbl <- tidyr::crossing(
      .receptor = receptor_ids,
      .repertoire = repertoire_ids
    ) |>
      dplyr::filter(stats::runif(dplyr::n()) < 0.07) |>
      dplyr::mutate(.weight = stats::rexp(dplyr::n(), rate = 0.5)) |>
      dplyr::transmute(
        !!receptor_col := .data$.receptor,
        !!repertoire_col := .data$.repertoire,
        !!count_col := ceiling(.data$.weight),
        !!prop_col := .data$.weight,
        cdr3_aa = .data$.receptor
      )
    rep_tbl <- tibble::tibble(
      !!repertoire_col := repertoire_ids,
      Group = repertoire_ids
    )

    abundance_matrix <- matrix(
      0,
      nrow = length(receptor_ids),
      ncol = length(repertoire_ids),
      dimnames = list(receptor_ids, repertoire_ids)
    )
    abundance_matrix[cbind(
      match(ann_tbl[[receptor_col]], receptor_ids),
      match(ann_tbl[[repertoire_col]], repertoire_ids)
    )] <- ann_tbl[[prop_col]]

    expected <- matrix(
      0,
      nrow = length(repertoire_ids),
      ncol = length(repertoire_ids),
      dimnames = list(repertoire_ids, repertoire_ids)
    )
    for (rep_i in seq_along(repertoire_ids)) {
      for (rep_j in seq_along(repertoire_ids)) {
        denominator <- sum(abundance_matrix[, rep_i] + abundance_matrix[, rep_j])
        expected[rep_i, rep_j] <- if (denominator > 0) {
          sum(abs(abundance_matrix[, rep_i] - abundance_matrix[, rep_j])) / denominator
        } else {
          NA_real_
        }
      }
    }

    idata <- immundata::ImmunData$new(
      schema = "cdr3_aa",
      annotations = duckplyr::as_duckdb_tibble(ann_tbl),
      repertoires = duckplyr::as_duckdb_tibble(rep_tbl)
    )
    idata$schema_repertoire <- "Group"

    expect_equal(
      unclass(repsim_bray(idata, autojoin = FALSE)),
      expected,
      tolerance = 1e-12,
      info = paste("seed", seed)
    )
  }
})
