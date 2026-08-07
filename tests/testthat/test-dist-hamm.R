make_dist_hamm_idata <- function() {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")

  annotations <- tibble::tibble(
    !!receptor_col := c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
    !!repertoire_col := c(1L, 1L, 1L, 1L, 1L, 1L, 2L),
    cdr3_aa = c("AAA", "AAT", "ATT", "CCCC", "AAA", "AAA", "AAA"),
    v_call = c("V1", "V1", "V1", "V1", "V2", "V1", "V1"),
    j_call = c("J1", "J1", "J1", "J1", "J1", "J1", "J1")
  ) |>
    duckplyr::as_duckdb_tibble()

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = c("cdr3_aa", "v_call", "j_call")
    ),
    annotations = annotations
  )
}


make_paired_dist_hamm_idata <- function() {
  input_data <- tibble::tribble(
    ~barcode, ~locus, ~cdr3_aa, ~v_call, ~umis,
    "cell1", "TRA", "CAVR",  "TRAV1", 10L,
    "cell1", "TRB", "CASSA", "TRBV1", 20L,
    "cell2", "TRA", "CAVG",  "TRAV2", 15L,
    "cell2", "TRB", "CASSB", "TRBV2", 25L
  )

  schema <- immundata::make_receptor_schema(
    features = c("cdr3_aa", "v_call"),
    chains = c("TRA", "TRB")
  )

  annotations <- immundata::agg_receptors(
    dataset = input_data,
    schema = schema,
    barcode_col = "barcode",
    locus_col = "locus",
    umi_col = "umis"
  )

  immundata::ImmunData$new(
    schema = schema,
    annotations = annotations
  )
}


test_that("dist_hamm returns lazy receptor-level upper-triangle edges", {
  idata <- make_dist_hamm_idata()

  lazy <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    autojoin = FALSE
  )

  expect_s3_class(lazy, "duckplyr_df")

  out <- lazy |>
    dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
    dplyr::collect()

  expect_named(
    out,
    c(
      "v_call", "j_call", "imd_receptor_id_1", "imd_receptor_id_2",
      "seq_len", "dist", "norm_dist", "sim"
    ),
    ignore.order = FALSE
  )
  expect_equal(
    out$imd_receptor_id_1,
    c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 6L)
  )
  expect_equal(
    out$imd_receptor_id_2,
    c(2L, 3L, 6L, 7L, 3L, 6L, 7L, 6L, 7L, 7L)
  )
  expect_equal(out$dist, c(1, 2, 0, 0, 1, 1, 1, 2, 2, 0))
  expect_equal(out$norm_dist, out$dist / 3)
  expect_equal(out$sim, 1 - out$norm_dist)
  expect_true(all(out$seq_len == 3L))
})


test_that("dist_hamm computes paired-chain distances separately by locus", {
  idata <- make_paired_dist_hamm_idata()

  out <- dist_hamm(idata, autojoin = FALSE) |>
    dplyr::arrange(.data$locus) |>
    dplyr::collect()

  expect_named(
    out,
    c(
      "locus", "imd_receptor_id_1", "imd_receptor_id_2",
      "seq_len", "dist", "norm_dist", "sim"
    ),
    ignore.order = FALSE
  )
  expect_equal(out$locus, c("TRA", "TRB"))
  expect_equal(out$imd_receptor_id_1, c(1L, 1L))
  expect_equal(out$imd_receptor_id_2, c(2L, 2L))
  expect_equal(out$seq_len, c(4L, 5L))
  expect_equal(out$dist, c(1, 1))
  expect_equal(out$norm_dist, c(0.25, 0.2))

  bounded <- dist_hamm(
    idata,
    by = "locus",
    max_dist = 0.21,
    autojoin = FALSE
  ) |>
    dplyr::collect()

  expect_equal(bounded$locus, "TRB")
  expect_equal(bounded$dist, 1)
  expect_equal(bounded$norm_dist, 0.2)
})


test_that("dist_hamm always returns normalized Hamming distance", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    autojoin = FALSE
  ) |>
    dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(
    out$norm_dist,
    c(1 / 3, 2 / 3, 0, 0, 1 / 3, 1 / 3, 1 / 3, 2 / 3, 2 / 3, 0)
  )
})


test_that("dist_hamm bounded modes use exact candidate generation", {
  idata <- make_dist_hamm_idata()

  raw <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1,
    autojoin = FALSE
  ) |>
    dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
    dplyr::collect()

  normalized <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1 / 3,
    autojoin = FALSE
  ) |>
    dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
    dplyr::collect()

  similar <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    min_sim = 2 / 3,
    autojoin = FALSE
  ) |>
    dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(
    raw$imd_receptor_id_1,
    c(1L, 1L, 1L, 2L, 2L, 2L, 6L)
  )
  expect_equal(
    raw$imd_receptor_id_2,
    c(2L, 6L, 7L, 3L, 6L, 7L, 7L)
  )
  expect_equal(raw$dist, c(1, 0, 0, 1, 1, 1, 0))
  expect_equal(normalized$imd_receptor_id_1, raw$imd_receptor_id_1)
  expect_equal(normalized$imd_receptor_id_2, raw$imd_receptor_id_2)
  expect_equal(normalized$dist, raw$dist)
  expect_equal(normalized$norm_dist, raw$dist / 3)
  expect_equal(similar, normalized)
})


test_that("dist_hamm similarity bounds include identity", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(idata, min_sim = 1, autojoin = FALSE) |>
    dplyr::collect()

  expect_true(all(out$dist == 0))
  expect_true(all(out$sim == 1))
})


test_that("dist_hamm grouping is optional and length matching is automatic", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(
    idata,
    by = NULL,
    max_dist = 1e-6,
    autojoin = FALSE
  ) |>
    dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(out$imd_receptor_id_1, c(1L, 1L, 1L, 5L, 5L, 6L))
  expect_equal(out$imd_receptor_id_2, c(5L, 6L, 7L, 6L, 7L, 7L))
  expect_true(all(out$dist == 0))
  expect_false(any(
    out$imd_receptor_id_1 == 4L | out$imd_receptor_id_2 == 4L
  ))
})


test_that("dist_hamm uses annotation columns in by", {
  idata <- make_dist_hamm_idata()
  repertoire_col <- immundata::imd_schema("repertoire")

  global <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1e-6,
    autojoin = FALSE
  ) |>
    dplyr::collect()

  within_repertoire <- dist_hamm(
    idata,
    by = c(repertoire_col, "v_call", "j_call"),
    max_dist = 1e-6,
    autojoin = FALSE
  ) |>
    dplyr::collect()

  expect_true(any(
    global$imd_receptor_id_1 == 1L & global$imd_receptor_id_2 == 7L
  ))
  expect_equal(within_repertoire$imd_receptor_id_1, 1L)
  expect_equal(within_repertoire$imd_receptor_id_2, 6L)
  expect_equal(within_repertoire[[repertoire_col]], 1L)
})


test_that("dist_hamm block search exactly matches filtered direct distance", {
  idata <- make_dist_hamm_idata()

  direct_raw <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    autojoin = FALSE
  ) |>
    dplyr::collect()

  for (bound in c(1, 2, 10)) {
    expected <- direct_raw |>
      dplyr::filter(.data$dist <= .env$bound) |>
      dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2)

    actual <- dist_hamm(
      idata,
      by = c("v_call", "j_call"),
      max_dist = bound,
      autojoin = FALSE
    ) |>
      dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
      dplyr::collect()

    expect_equal(actual, expected)
  }

  for (bound in c(1e-6, 1 / 3, 2 / 3)) {
    expected <- direct_raw |>
      dplyr::filter(.data$norm_dist <= .env$bound + 1e-12) |>
      dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2)

    actual <- dist_hamm(
      idata,
      by = c("v_call", "j_call"),
      max_dist = bound,
      autojoin = FALSE
    ) |>
      dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2) |>
      dplyr::collect()

    expect_equal(actual, expected, tolerance = 1e-12)
  }
})


test_that("dist_hamm does not require repertoire aggregation", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 0
  )

  expect_s3_class(out, "duckplyr_df")
  expect_null(idata$repertoires)
  expect_false(immundata::imd_schema("repertoire") %in% colnames(out))
})


test_that("dist_hamm validates columns and bounds", {
  idata <- make_dist_hamm_idata()

  expect_error(dist_hamm(idata, seq_col = "missing"), "missing")
  expect_error(dist_hamm(idata, by = "missing"), "missing")
  expect_error(dist_hamm(idata, max_dist = 1.5), "integer")
  expect_error(dist_hamm(idata, min_sim = 1.1), "<= 1")
  expect_error(
    dist_hamm(idata, max_dist = 0.1, min_sim = 0.9),
    "only one"
  )

  idata_without_sequence_schema <- immundata::ImmunData$new(
    schema = "v_call",
    annotations = idata$annotations
  )
  expect_error(
    dist_hamm(idata_without_sequence_schema, seq_col = "cdr3_aa"),
    "receptor schema"
  )
})


test_that("dist_hamm is registered in the dist family", {
  expect_true("hamm" %in% ls(IMMUNARCH_METHOD_REGISTRY[["dist"]]))
})
