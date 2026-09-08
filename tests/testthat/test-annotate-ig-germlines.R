make_ig_germline_idata <- function(rows, chains = unique(rows$locus)) {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")

  if (!receptor_col %in% names(rows)) {
    rows[[receptor_col]] <- seq_len(nrow(rows))
  }
  if (!repertoire_col %in% names(rows)) {
    rows[[repertoire_col]] <- 1L
  }

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = "row_id",
      chains = chains
    ),
    annotations = duckplyr::as_duckdb_tibble(
      rows,
      prudence = "stingy"
    )
  )
}


collect_germlines <- function(idata, columns) {
  idata$annotations |>
    dplyr::select(dplyr::all_of(columns)) |>
    dplyr::collect() |>
    dplyr::arrange(row_id)
}


testthat::test_that("germlines mask the complete interval between V and J", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGH",
    germline_alignment = "AAAACCCCGGGG",
    v_alignment_end = 4L,
    j_alignment_start = 9L
  ))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(actual$germline_alignment_d_mask, "AAAANNNNGGGG")
})


testthat::test_that("light-chain germlines need no D columns", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGK",
    germline_alignment = "AAAACCCCGGGG",
    v_alignment_end = 4L,
    j_alignment_start = 9L
  ))

  testthat::expect_false(any(c(
    "d_call", "d_alignment_start", "d_alignment_end"
  ) %in% names(idata$annotations)))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(actual$germline_alignment_d_mask, "AAAANNNNGGGG")
})


testthat::test_that("adjacent V and J boundaries leave germlines unchanged", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGH",
    germline_alignment = "AAAAGGGG",
    v_alignment_end = 4L,
    j_alignment_start = 5L
  ))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(actual$germline_alignment_d_mask, "AAAAGGGG")
})


testthat::test_that("germlines are masked independently row by row", {
  idata <- make_ig_germline_idata(tibble::tribble(
    ~row_id, ~locus, ~germline_alignment, ~v_alignment_end, ~j_alignment_start,
    1L, "IGH", "AAAACCCCGGGG", 4L, 9L,
    2L, "IGH", "AATTCCGG", 2L, 7L,
    3L, "IGL", "GGGTTT", 3L, 4L
  ))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(
    actual$germline_alignment_d_mask,
    c("AAAANNNNGGGG", "AANNNNGG", "GGGTTT")
  )
})


testthat::test_that("missing row values produce missing masked germlines", {
  idata <- make_ig_germline_idata(tibble::tribble(
    ~row_id, ~locus, ~germline_alignment, ~v_alignment_end, ~j_alignment_start,
    1L, "IGH", NA_character_, 4L, 9L,
    2L, "IGH", "AAAACCCCGGGG", NA_integer_, 9L,
    3L, "IGH", "AAAACCCCGGGG", 4L, NA_integer_,
    4L, "IGH", "AAAACCCCGGGG", 4L, 9L
  ))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(nrow(actual), 4L)
  testthat::expect_equal(
    actual$germline_alignment_d_mask,
    c(NA_character_, NA_character_, NA_character_, "AAAANNNNGGGG")
  )
})


testthat::test_that("an existing masked-germline column is replaced", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGH",
    germline_alignment = "AAAACCCCGGGG",
    v_alignment_end = 4L,
    j_alignment_start = 9L,
    germline_alignment_d_mask = "incorrect"
  ))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(actual$germline_alignment_d_mask, "AAAANNNNGGGG")
  testthat::expect_equal(sum(names(out$annotations) == "germline_alignment_d_mask"), 1L)
})


testthat::test_that("custom germline and coordinate columns are supported", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGH",
    custom_germline = "AAAACCCCGGGG",
    custom_v_end = 4L,
    custom_j_start = 9L
  ))

  out <- annotate_ig_germlines(
    idata,
    germline_col = "custom_germline",
    v_end_col = "custom_v_end",
    j_start_col = "custom_j_start",
    d_mask_col = "custom_d_mask"
  )
  actual <- collect_germlines(out, c("row_id", "custom_d_mask"))

  testthat::expect_equal(actual$custom_d_mask, "AAAANNNNGGGG")
  testthat::expect_false("germline_alignment_d_mask" %in% names(out$annotations))
})


testthat::test_that("missing required germline columns use annotation validation", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGH",
    germline_alignment = "AAAACCCCGGGG",
    v_alignment_end = 4L
  ))

  testthat::expect_error(
    annotate_ig_germlines(idata),
    "Passed column name.*j_alignment_start.*not in the input ImmunData"
  )
})


testthat::test_that("germline masking has no clone-family dependency", {
  idata <- make_ig_germline_idata(tibble::tibble(
    row_id = 1L,
    locus = "IGH",
    germline_alignment = "AAAACCCCGGGG",
    v_alignment_end = 4L,
    j_alignment_start = 9L
  ))

  testthat::expect_false(any(c(
    "family_id", "heavy_family_id", "light_family_id"
  ) %in% names(idata$annotations)))

  out <- annotate_ig_germlines(idata)
  actual <- collect_germlines(out, c("row_id", "germline_alignment_d_mask"))

  testthat::expect_equal(actual$germline_alignment_d_mask, "AAAANNNNGGGG")
})


testthat::test_that("germline masking is independent of clone annotation order", {
  receptor_col <- immundata::imd_schema("receptor")
  idata <- make_ig_germline_idata(
    tibble::tibble(
      row_id = 1:4,
      locus = rep(c("IGH", "IGK"), 2),
      germline_alignment = c(
        "AAAACCCCGGGG", "TTTTAAAACCCC",
        "AAAACCCCGGGG", "TTTTAAAACCCC"
      ),
      v_alignment_end = 4L,
      j_alignment_start = 9L,
      junction = c("AAAA", "CCCC", "AAAT", "CCCT"),
      v_call = rep(c("IGHV1*01", "IGKV1*01"), 2),
      j_call = rep(c("IGHJ1*01", "IGKJ1*01"), 2),
      subject = "A",
      !!receptor_col := rep(1:2, each = 2)
    ),
    chains = c("IGH", "IGK")
  )

  germlines_first <- idata$clone(deep = TRUE) |>
    annotate_ig_germlines() |>
    annotate_ig_clones(heavy_threshold = 0.25, by_scope = "subject")
  clones_first <- idata$clone(deep = TRUE) |>
    annotate_ig_clones(heavy_threshold = 0.25, by_scope = "subject") |>
    annotate_ig_germlines()

  first <- collect_germlines(
    germlines_first,
    c("row_id", "germline_alignment_d_mask")
  )
  second <- collect_germlines(
    clones_first,
    c("row_id", "germline_alignment_d_mask")
  )

  testthat::expect_equal(first, second)
})


testthat::test_that("germline masking preserves ImmunData structure and schema", {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  idata <- make_ig_germline_idata(tibble::tribble(
    ~row_id, ~locus, ~germline_alignment, ~v_alignment_end,
    ~j_alignment_start, ~existing_annotation, ~imd_receptor_id,
    ~imd_repertoire_id,
    2L, "IGH", "AATTCCGG", 2L, 7L, "keep-b", 102L, 12L,
    1L, "IGH", "AAAACCCCGGGG", 4L, 9L, "keep-a", 101L, 11L
  ))
  before <- idata$annotations |>
    dplyr::select(dplyr::all_of(c(
      "row_id", receptor_col, repertoire_col, "existing_annotation"
    ))) |>
    dplyr::collect() |>
    dplyr::arrange(row_id)

  out <- annotate_ig_germlines(idata)
  after <- out$annotations |>
    dplyr::select(dplyr::all_of(c(
      "row_id", receptor_col, repertoire_col, "existing_annotation"
    ))) |>
    dplyr::collect() |>
    dplyr::arrange(row_id)

  testthat::expect_s3_class(out, "ImmunData")
  testthat::expect_s3_class(out$annotations, "duckplyr_df")
  testthat::expect_equal(nrow(after), nrow(before))
  testthat::expect_equal(after, before)
  testthat::expect_true(all(names(idata$annotations) %in% names(out$annotations)))
})
