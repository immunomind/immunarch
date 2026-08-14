sanitize_strata_id_for_test <- function(sid) {
  sid_chr <- as.character(sid)
  sid_chr <- gsub("[^A-Za-z0-9]+", "_", sid_chr)
  sid_chr <- gsub("^_+|_+$", "", sid_chr)
  if (sid_chr == "") sid_chr <- "na"
  sid_chr
}

build_synthetic_public_idata <- function(count_matrix, prop_matrix, strata_ids) {
  stopifnot(identical(dim(count_matrix), dim(prop_matrix)))
  stopifnot(length(strata_ids) == ncol(count_matrix))

  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  strata_col <- immundata::imd_schema("strata")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")

  receptor_ids <- rownames(count_matrix)
  repertoire_ids <- colnames(count_matrix)

  if (is.null(receptor_ids)) receptor_ids <- as.character(seq_len(nrow(count_matrix)))
  if (is.null(repertoire_ids)) repertoire_ids <- as.character(seq_len(ncol(count_matrix)))

  present_pairs <- which(count_matrix > 0, arr.ind = TRUE)
  ann_tbl <- tibble::tibble(
    !!receptor_col := receptor_ids[present_pairs[, 1]],
    !!repertoire_col := repertoire_ids[present_pairs[, 2]],
    !!count_col := as.numeric(count_matrix[present_pairs]),
    !!prop_col := as.numeric(prop_matrix[present_pairs]),
    cdr3_aa = receptor_ids[present_pairs[, 1]]
  )

  rep_tbl <- tibble::tibble(
    !!repertoire_col := repertoire_ids,
    !!strata_col := strata_ids
  )

  ann_tbl <- duckplyr::as_duckdb_tibble(ann_tbl, prudence = "stingy")
  rep_tbl <- duckplyr::as_duckdb_tibble(rep_tbl, prudence = "stingy")

  immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = ann_tbl,
    repertoires = rep_tbl
  )
}

compute_public_oracle_from_matrix <- function(count_matrix, prop_matrix, strata_ids = NULL) {
  stopifnot(identical(dim(count_matrix), dim(prop_matrix)))
  if (!is.null(strata_ids)) {
    stopifnot(length(strata_ids) == ncol(count_matrix))
  }

  receptor_col <- immundata::imd_schema("receptor")
  receptor_ids <- rownames(count_matrix)
  if (is.null(receptor_ids)) receptor_ids <- as.character(seq_len(nrow(count_matrix)))

  compute_block <- function(count_sub, prop_sub, n_repertoires_sub) {
    block_rows <- vector("list", nrow(count_sub))

    for (i in seq_len(nrow(count_sub))) {
      present <- count_sub[i, ] > 0
      count_values <- as.numeric(count_sub[i, present])
      prop_values <- as.numeric(prop_sub[i, present])

      if (length(count_values) == 0) {
        block_rows[[i]] <- tibble::tibble(
          imd_public_incidence = NA_real_,
          imd_public_count_min = NA_real_,
          imd_public_count_max = NA_real_,
          imd_public_count_mean = NA_real_,
          imd_public_count_median = NA_real_,
          imd_public_prop_min = NA_real_,
          imd_public_prop_max = NA_real_,
          imd_public_prop_mean = NA_real_,
          imd_public_prop_median = NA_real_,
          imd_public_incidence_prop = NA_real_
        )
      } else {
        block_rows[[i]] <- tibble::tibble(
          imd_public_incidence = as.numeric(length(count_values)),
          imd_public_count_min = min(count_values),
          imd_public_count_max = max(count_values),
          imd_public_count_mean = mean(count_values),
          imd_public_count_median = stats::median(count_values),
          imd_public_prop_min = min(prop_values),
          imd_public_prop_max = max(prop_values),
          imd_public_prop_mean = mean(prop_values),
          imd_public_prop_median = stats::median(prop_values),
          imd_public_incidence_prop = length(count_values) / n_repertoires_sub
        )
      }
    }

    dplyr::bind_rows(block_rows)
  }

  expected <- tibble::tibble(!!receptor_col := receptor_ids)
  expected <- dplyr::bind_cols(
    expected,
    compute_block(count_matrix, prop_matrix, ncol(count_matrix))
  )

  if (!is.null(strata_ids)) {
    for (sid in unique(strata_ids)) {
      idx <- which(strata_ids == sid)
      sid_block <- compute_block(
        count_matrix[, idx, drop = FALSE],
        prop_matrix[, idx, drop = FALSE],
        length(idx)
      )
      sid_block[[receptor_col]] <- receptor_ids
      sid_block <- sid_block |>
        dplyr::filter(!is.na(.data$imd_public_incidence))

      sid_metric_cols <- setdiff(names(sid_block), receptor_col)
      sid_suffix <- sanitize_strata_id_for_test(sid)
      col_idx <- match(sid_metric_cols, colnames(sid_block))
      colnames(sid_block)[col_idx] <- paste0(sid_metric_cols, "_strata_", sid_suffix)

      expected <- expected |>
        dplyr::left_join(sid_block, by = receptor_col)
    }
  }

  expected
}

build_public_matrices_from_idata <- function(idata) {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  strata_col <- immundata::imd_schema("strata")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")

  ann <- dplyr::collect(idata$annotations)
  rep_tbl <- idata$repertoires

  repertoire_ids <- as.character(rep_tbl[[repertoire_col]])

  count_wide <- ann |>
    dplyr::select(all_of(c(receptor_col, repertoire_col, count_col))) |>
    tidyr::pivot_wider(
      names_from = all_of(repertoire_col),
      values_from = all_of(count_col),
      values_fn = max,
      values_fill = 0
    )

  prop_wide <- ann |>
    dplyr::select(all_of(c(receptor_col, repertoire_col, prop_col))) |>
    tidyr::pivot_wider(
      names_from = all_of(repertoire_col),
      values_from = all_of(prop_col),
      values_fn = max,
      values_fill = 0
    )

  receptor_ids <- count_wide[[receptor_col]]
  count_matrix <- as.matrix(count_wide[, repertoire_ids, drop = FALSE])
  prop_matrix <- as.matrix(prop_wide[, repertoire_ids, drop = FALSE])

  rownames(count_matrix) <- as.character(receptor_ids)
  rownames(prop_matrix) <- as.character(receptor_ids)
  colnames(count_matrix) <- repertoire_ids
  colnames(prop_matrix) <- repertoire_ids

  list(
    count_matrix = count_matrix,
    prop_matrix = prop_matrix,
    strata_ids = if (strata_col %in% names(rep_tbl)) rep_tbl[[strata_col]] else NULL,
    receptor_ids = receptor_ids
  )
}


make_small_public_matrices <- function() {
  receptor_ids <- c("R1", "R2", "R3")
  repertoire_ids <- c("P1", "P2", "P3")

  list(
    counts = matrix(
      c(
        10, 0, 5,
        1, 2, 3,
        0, 4, 0
      ),
      nrow = length(receptor_ids),
      byrow = TRUE,
      dimnames = list(receptor_ids, repertoire_ids)
    ),
    proportions = matrix(
      c(
        0.50, 0.00, 0.25,
        0.10, 0.20, 0.30,
        0.00, 0.40, 0.00
      ),
      nrow = length(receptor_ids),
      byrow = TRUE,
      dimnames = list(receptor_ids, repertoire_ids)
    ),
    strata = c(1L, 1L, 2L)
  )
}

test_that("annotate_public adds global publicness metrics", {

  idata <- make_aggregated_test_idata()
  out <- annotate_public(idata)
  input_ann <- dplyr::collect(idata$annotations)
  out_ann <- dplyr::collect(out$annotations)

  expect_equal(nrow(out_ann), nrow(input_ann))
  expect_true(all(names(input_ann) %in% names(out_ann)))
  expect_true(all(c(
    "imd_public_incidence",
    "imd_public_incidence_prop",
    "imd_public_count_min",
    "imd_public_count_max",
    "imd_public_count_mean",
    "imd_public_count_median",
    "imd_public_prop_min",
    "imd_public_prop_max",
    "imd_public_prop_mean",
    "imd_public_prop_median"
  ) %in% names(out_ann)))

  expect_equal(length(grep("_strata_", names(out_ann))), 0)
})

test_that("annotate_public computes global metrics correctly", {

  idata <- make_aggregated_test_idata()
  out <- annotate_public(idata)

  matrices <- build_public_matrices_from_idata(idata)
  expected <- compute_public_oracle_from_matrix(
    matrices$count_matrix,
    matrices$prop_matrix,
    strata_ids = NULL
  )

  receptor_col <- immundata::imd_schema("receptor")
  observed <- out$annotations |>
    dplyr::select(all_of(names(expected))) |>
    dplyr::distinct() |>
    dplyr::collect()

  observed[[receptor_col]] <- as.character(observed[[receptor_col]])
  expected[[receptor_col]] <- as.character(expected[[receptor_col]])

  observed <- observed |>
    dplyr::arrange(.data[[receptor_col]])

  expected <- expected |>
    dplyr::arrange(.data[[receptor_col]])

  metric_cols <- setdiff(names(expected), receptor_col)
  observed[metric_cols] <- lapply(observed[metric_cols], as.numeric)
  expected[metric_cols] <- lapply(expected[metric_cols], as.numeric)

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("annotate_public adds per-strata metrics when idata is stratified", {

  idata <- make_aggregated_test_idata()
  stratified <- immundata::agg_strata(
    idata,
    schema = c("Response", "Therapy")
  )
  out <- annotate_public(stratified)
  out_ann <- dplyr::collect(out$annotations)

  strata_col <- immundata::imd_schema("strata")
  expect_true(strata_col %in% names(stratified$repertoires))

  n_strata <- stratified$repertoires |>
    dplyr::distinct(.data[[strata_col]]) |>
    nrow()

  inc_cols_n <- length(grep("^imd_public_incidence_strata_", names(out_ann)))
  prop_mean_cols_n <- length(grep("^imd_public_prop_mean_strata_", names(out_ann)))

  expect_equal(inc_cols_n, n_strata)
  expect_equal(prop_mean_cols_n, n_strata)
  expect_equal(out$schema_receptor, stratified$schema_receptor)
  expect_equal(out$schema_repertoire, stratified$schema_repertoire)
  expect_equal(out$schema_strata, stratified$schema_strata)
  expect_equal(
    dplyr::collect(out$repertoires),
    dplyr::collect(stratified$repertoires)
  )
  expect_equal(
    dplyr::collect(out$strata),
    dplyr::collect(stratified$strata)
  )
  expect_identical(out$provenance, stratified$provenance)
})

test_that("annotate_public computes per-strata metrics correctly", {

  idata <- make_aggregated_test_idata() |>
    immundata::agg_strata(schema = "Response")
  out <- annotate_public(idata)

  matrices <- build_public_matrices_from_idata(idata)
  expected <- compute_public_oracle_from_matrix(
    matrices$count_matrix,
    matrices$prop_matrix,
    strata_ids = matrices$strata_ids
  )

  receptor_col <- immundata::imd_schema("receptor")
  observed <- out$annotations |>
    dplyr::select(all_of(names(expected))) |>
    dplyr::distinct() |>
    dplyr::collect()

  observed[[receptor_col]] <- as.character(observed[[receptor_col]])
  expected[[receptor_col]] <- as.character(expected[[receptor_col]])

  observed <- observed |>
    dplyr::arrange(.data[[receptor_col]])

  expected <- expected |>
    dplyr::arrange(.data[[receptor_col]])

  metric_cols <- setdiff(names(expected), receptor_col)
  observed[metric_cols] <- lapply(observed[metric_cols], as.numeric)
  expected[metric_cols] <- lapply(expected[metric_cols], as.numeric)

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("annotate_public errors if repertoires are not aggregated", {

  source_idata <- make_test_idata()
  idata <- immundata::ImmunData$new(
    schema = source_idata$schema_receptor,
    annotations = source_idata$annotations
  )

  expect_error(
    annotate_public(idata),
    "agg_repertoires"
  )
})

test_that("annotate_public errors when required annotation columns are missing", {

  count_col <- immundata::imd_schema("count")
  idata <- make_aggregated_test_idata()

  bad_annotations <- idata$annotations |>
    dplyr::select(-all_of(count_col))

  bad_idata <- immundata::ImmunData$new(
    schema = idata$schema_receptor,
    annotations = bad_annotations,
    repertoires = idata$repertoires
  )

  expect_error(
    annotate_public(bad_idata),
    "missing in .*idata\\$annotations"
  )
})

test_that("annotate_public errors when repertoires table is empty", {

  idata <- make_aggregated_test_idata()

  empty_repertoires <- idata$repertoires[0, , drop = FALSE]

  bad_idata <- immundata::ImmunData$new(
    schema = idata$schema_receptor,
    annotations = idata$annotations,
    repertoires = empty_repertoires
  )

  expect_error(
    annotate_public(bad_idata),
    "No repertoires found"
  )
})

test_that("annotate_public matches deterministic matrix oracle", {

  receptor_ids <- paste0("R", 1:5)
  repertoire_ids <- paste0("P", 1:4)

  count_matrix <- matrix(
    c(
      10, 0, 4, 0,
      0, 6, 0, 0,
      7, 8, 9, 10,
      0, 0, 2, 2,
      1, 1, 1, 1
    ),
    nrow = length(receptor_ids),
    byrow = TRUE,
    dimnames = list(receptor_ids, repertoire_ids)
  )

  prop_matrix <- matrix(
    c(
      0.50, 0.00, 0.20, 0.00,
      0.00, 0.60, 0.00, 0.00,
      0.70, 0.80, 0.90, 1.00,
      0.00, 0.00, 0.05, 0.03,
      0.01, 0.02, 0.03, 0.04
    ),
    nrow = length(receptor_ids),
    byrow = TRUE,
    dimnames = list(receptor_ids, repertoire_ids)
  )

  strata_ids <- c(1L, 1L, 2L, 2L)

  idata <- build_synthetic_public_idata(count_matrix, prop_matrix, strata_ids)
  out <- annotate_public(idata)

  expected <- compute_public_oracle_from_matrix(count_matrix, prop_matrix, strata_ids)
  receptor_col <- immundata::imd_schema("receptor")

  observed <- out$annotations |>
    dplyr::select(all_of(names(expected))) |>
    dplyr::distinct() |>
    dplyr::collect()

  observed[[receptor_col]] <- as.character(observed[[receptor_col]])
  expected[[receptor_col]] <- as.character(expected[[receptor_col]])

  observed <- observed |>
    dplyr::arrange(.data[[receptor_col]])

  expected <- expected |>
    dplyr::arrange(.data[[receptor_col]])

  metric_cols <- setdiff(names(expected), receptor_col)
  observed[metric_cols] <- lapply(observed[metric_cols], as.numeric)
  expected[metric_cols] <- lapply(expected[metric_cols], as.numeric)

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("annotate_public reports global-only mode and computes deterministic no-strata metrics", {

  receptor_col <- immundata::imd_schema("receptor")
  strata_col <- immundata::imd_schema("strata")
  matrices <- make_small_public_matrices()
  receptor_ids <- rownames(matrices$counts)

  idata_with_strata <- build_synthetic_public_idata(
    count_matrix = matrices$counts,
    prop_matrix = matrices$proportions,
    strata_ids = matrices$strata
  )

  rep_tbl_no_strata <- idata_with_strata$repertoires |>
    dplyr::select(-all_of(strata_col))
  rep_tbl_no_strata <- duckplyr::as_duckdb_tibble(
    rep_tbl_no_strata,
    prudence = "stingy"
  )

  idata_no_strata <- immundata::ImmunData$new(
    schema = idata_with_strata$schema_receptor,
    annotations = idata_with_strata$annotations,
    repertoires = rep_tbl_no_strata
  )

  expect_message(
    out <- annotate_public(idata_no_strata),
    "global publicness metrics only"
  )

  expected <- tibble::tibble(
    !!receptor_col := receptor_ids,
    imd_public_incidence = c(2, 3, 1),
    imd_public_incidence_prop = c(2 / 3, 1, 1 / 3),
    imd_public_count_mean = c(7.5, 2.0, 4.0),
    imd_public_prop_mean = c(0.375, 0.2, 0.4)
  )

  observed <- out$annotations |>
    dplyr::select(all_of(names(expected))) |>
    dplyr::distinct() |>
    dplyr::collect()

  observed[[receptor_col]] <- as.character(observed[[receptor_col]])
  expected[[receptor_col]] <- as.character(expected[[receptor_col]])

  observed <- observed |>
    dplyr::arrange(.data[[receptor_col]])
  expected <- expected |>
    dplyr::arrange(.data[[receptor_col]])

  metric_cols <- setdiff(names(expected), receptor_col)
  observed[metric_cols] <- lapply(observed[metric_cols], as.numeric)
  expected[metric_cols] <- lapply(expected[metric_cols], as.numeric)

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("annotate_public computes deterministic global and per-strata incidence and proportion metrics", {

  receptor_col <- immundata::imd_schema("receptor")
  matrices <- make_small_public_matrices()
  receptor_ids <- rownames(matrices$counts)

  idata <- build_synthetic_public_idata(
    count_matrix = matrices$counts,
    prop_matrix = matrices$proportions,
    strata_ids = matrices$strata
  )

  out <- annotate_public(idata)

  expected <- tibble::tibble(
    !!receptor_col := receptor_ids,
    imd_public_incidence = c(2, 3, 1),
    imd_public_incidence_prop = c(2 / 3, 1, 1 / 3),
    imd_public_count_mean = c(7.5, 2.0, 4.0),
    imd_public_prop_mean = c(0.375, 0.2, 0.4),
    imd_public_incidence_strata_1 = c(1, 2, 1),
    imd_public_incidence_prop_strata_1 = c(0.5, 1.0, 0.5),
    imd_public_count_mean_strata_1 = c(10.0, 1.5, 4.0),
    imd_public_prop_mean_strata_1 = c(0.5, 0.15, 0.4),
    imd_public_incidence_strata_2 = c(1, 1, NA_real_),
    imd_public_incidence_prop_strata_2 = c(1, 1, NA_real_),
    imd_public_count_mean_strata_2 = c(5.0, 3.0, NA_real_),
    imd_public_prop_mean_strata_2 = c(0.25, 0.30, NA_real_)
  )

  observed <- out$annotations |>
    dplyr::select(all_of(names(expected))) |>
    dplyr::distinct() |>
    dplyr::collect()

  observed[[receptor_col]] <- as.character(observed[[receptor_col]])
  expected[[receptor_col]] <- as.character(expected[[receptor_col]])

  observed <- observed |>
    dplyr::arrange(.data[[receptor_col]])
  expected <- expected |>
    dplyr::arrange(.data[[receptor_col]])

  metric_cols <- setdiff(names(expected), receptor_col)
  observed[metric_cols] <- lapply(observed[metric_cols], as.numeric)
  expected[metric_cols] <- lapply(expected[metric_cols], as.numeric)

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("annotate_public replaces existing publicness annotations", {

  receptor_col <- immundata::imd_schema("receptor")
  matrices <- make_small_public_matrices()

  idata <- build_synthetic_public_idata(
    count_matrix = matrices$counts,
    prop_matrix = matrices$proportions,
    strata_ids = matrices$strata
  )

  once <- annotate_public(idata)
  twice <- annotate_public(once)

  once_ann <- dplyr::collect(once$annotations)
  twice_ann <- dplyr::collect(twice$annotations)
  once_public_cols <- grep("^imd_public_", names(once_ann), value = TRUE)
  twice_public_cols <- grep("^imd_public_", names(twice_ann), value = TRUE)

  expect_setequal(twice_public_cols, once_public_cols)
  expect_false(any(grepl("\\.[xy]$", names(twice_ann))))

  once_metrics <- once_ann |>
    dplyr::select(all_of(c(receptor_col, once_public_cols))) |>
    dplyr::distinct() |>
    dplyr::arrange(.data[[receptor_col]])

  twice_metrics <- twice_ann |>
    dplyr::select(all_of(c(receptor_col, twice_public_cols))) |>
    dplyr::distinct() |>
    dplyr::arrange(.data[[receptor_col]])

  expect_equal(twice_metrics, once_metrics)
})
