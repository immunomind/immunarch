make_ig_germline_idata <- function(rows) {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")

  if (!receptor_col %in% names(rows)) rows[[receptor_col]] <- seq_len(nrow(rows))
  if (!repertoire_col %in% names(rows)) rows[[repertoire_col]] <- 1L

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(features = "row_id", chains = unique(rows$locus)),
    annotations = duckplyr::as_duckdb_tibble(rows, prudence = "stingy")
  )
}


collect_germlines <- function(idata) {
  idata$annotations |>
    dplyr::collect() |>
    dplyr::arrange(row_id)
}


germline_references <- function() {
  tibble::tribble(
    ~species, ~locus, ~segment, ~allele_id, ~sequence,
    "Homo_sapiens", "IGH", "V", "IGHV1*01", "AA..CC",
    "Homo_sapiens", "IGH", "V", "IGHV2*01", "TT..CC",
    "Homo_sapiens", "IGH", "J", "IGHJ1*01", "GG.TT",
    "Homo_sapiens", "IGH", "J", "IGHJ2*01", "CC.AA",
    "Mus_musculus", "IGH", "V", "IGHV1*01", "CCCCCC",
    "Mus_musculus", "IGH", "J", "IGHJ1*01", "AAAAAA"
  )
}


germline_rows <- function() {
  tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end, ~keep, ~germline_alignment_d_mask,
    1L, "z-long", "tied", "IGH", "AAAAAAAAAA", "IGHV1*01,IGHV99*01", "IGHJ1*01|IGHJ99*01", 2L, 4L, 1L, 4L, "z", "legacy-z",
    2L, "a-long", "tied", "IGH", "AAAAAAAAAA", "IGHV2*01", "IGHJ2*01", 1L, 4L, 1L, 4L, "a", "legacy-a",
    3L, "b-short", "tied", "IGH", "AAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L, "b", "legacy-b",
    4L, "c-short", "tied", "IGH", "AAAAAAAAA", "IGHV2*01", "IGHJ2*01", 1L, 4L, 1L, 4L, "c", "legacy-c",
    5L, "x-long", "fallback", "IGH", "AAAAAAAAAAAA", "IGHV1*01", "IGHJ2*01", 1L, 4L, 1L, 4L, "x", "legacy-x",
    6L, "z-short", "fallback", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L, "z", "legacy-z2",
    7L, "a-short", "fallback", "IGH", "AAAAAAAAAA", "IGHV2*01", "IGHJ1*01", 1L, 4L, 1L, 4L, "a", "legacy-a2"
  )
}


testthat::test_that("germlines use modal V/J representatives and reconstruct V + N + J", {
  rows <- germline_rows()
  out <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(),
    species = "Homo_sapiens"
  )
  actual <- collect_germlines(out)

  # Tied modal V/J calls retain either allele. The longest candidates tie, so
  # the character sequence_id selects a-long and its V2/J2 coordinates.
  testthat::expect_equal(actual$germline_alignment[actual$clone_id == "tied"],
    rep("TT..NNCC.A", 4L))
  # No longest row has the modal V1/J1 combination. The shorter fallback is
  # selected, extended at J to the clone maximum, and padded past reference end.
  testthat::expect_equal(actual$germline_alignment[actual$clone_id == "fallback"],
    rep("AA..NNGG.TTN", 3L))
})


testthat::test_that("germline reconstruction selects one species and auto-groups by clone and locus", {
  rows <- germline_rows()[1:2, ]
  refs <- germline_references()

  testthat::expect_error(
    annotate_ig_germlines(make_ig_germline_idata(rows), references = refs),
    "species"
  )
  testthat::expect_error(
    annotate_ig_germlines(make_ig_germline_idata(rows), references = refs, species = "Rattus_norvegicus"),
    "species"
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows), references = refs, species = "Homo_sapiens"
  ) |>
    collect_germlines()
  testthat::expect_equal(actual$germline_alignment, rep("TT..NNCC.A", 2L))
})


testthat::test_that("failed germline groups retain rows and preserve existing annotations", {
  rows <- germline_rows()[1:2, ]
  rows$v_call <- "IGHV404*01"
  rows$germline_alignment <- "obsolete"
  before <- rows[, c("row_id", "sequence_alignment", "v_call", "j_call", "keep", "germline_alignment_d_mask")]

  out <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(), species = "Homo_sapiens"
  )
  actual <- collect_germlines(out)

  testthat::expect_equal(nrow(actual), nrow(rows))
  testthat::expect_true(all(is.na(actual$germline_alignment)))
  testthat::expect_equal(
    actual[, names(before)],
    before[match(actual$row_id, before$row_id), ]
  )
  testthat::expect_s3_class(out$annotations, "duckplyr_df")
})


testthat::test_that("maximum alignment length outranks a smaller ID among modal candidates", {
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, "a-short", "length", "IGH", "AAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    2L, "z-long", "length", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ1*01", 2L, 4L, 1L, 4L
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(), species = "Homo_sapiens"
  ) |>
    collect_germlines()

  testthat::expect_equal(actual$germline_alignment, rep("A..NNNGG.T", 2L))
})


testthat::test_that("shorter fallback candidates compare sequence IDs as characters", {
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, 99, "fallback-tie", "IGH", "AAAAAAAAAAAA", "IGHV1*01", "IGHJ2*01", 1L, 4L, 1L, 4L,
    2L, 2, "fallback-tie", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ1*01", 2L, 4L, 1L, 4L,
    3L, 10, "fallback-tie", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(), species = "Homo_sapiens"
  ) |>
    collect_germlines()

  testthat::expect_equal(actual$germline_alignment, rep("AA..NNGG.TTN", 3L))
})


testthat::test_that("IMGT gaps are preserved from the reference sequence", {
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, "gapped", "gapped", "IGH", "AA.CCNNGGTT", "IGHV1*01", "IGHJ1*01", 1L, 5L, 1L, 4L
  )
  refs <- tibble::tribble(
    ~species, ~locus, ~segment, ~allele_id, ~sequence,
    "Homo_sapiens", "IGH", "V", "IGHV1*01", "AA.CCAA",
    "Homo_sapiens", "IGH", "J", "IGHJ1*01", "GGTTAA"
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows), references = refs
  ) |>
    collect_germlines()

  testthat::expect_equal(actual$germline_alignment, "AA.CCNNGGTT")
})


testthat::test_that("automatic clone grouping remains locus-scoped", {
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, "heavy", "same-family", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    2L, "light", "same-family", "IGK", "AAAAAAAAAA", "IGKV1*01", "IGKJ1*01", 1L, 4L, 1L, 4L
  )
  refs <- dplyr::bind_rows(
    germline_references(),
    tibble::tribble(
      ~species, ~locus, ~segment, ~allele_id, ~sequence,
      "Homo_sapiens", "IGK", "V", "IGKV1*01", "CCCCAA",
      "Homo_sapiens", "IGK", "J", "IGKJ1*01", "TTTTAA"
    )
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows), references = refs, species = "Homo_sapiens"
  ) |>
    collect_germlines()

  testthat::expect_equal(actual$germline_alignment, c("AA..NNGG.T", "CCCCNNTTTT"))
})


testthat::test_that("automatic clone grouping remains subject-scoped", {
  rows <- germline_rows()[1:2, ]
  rows$subject_id <- c("subject1", "subject2")

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(), species = "Homo_sapiens"
  ) |>
    collect_germlines()

  testthat::expect_equal(
    actual$germline_alignment,
    c("A..NNNGG.T", "TT..NNCC.A")
  )
})


testthat::test_that("modal V/J calls are counted before alignment length preference", {
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, "a-short", "votes", "IGH", "AAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    2L, "b-short", "votes", "IGH", "AAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    3L, "c-short", "votes", "IGH", "AAAAAAAAA", "IGHV1*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    4L, "z-long", "votes", "IGH", "AAAAAAAAAA", "IGHV2*01", "IGHJ2*01", 1L, 4L, 1L, 4L
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(), species = "Homo_sapiens"
  ) |>
    collect_germlines()

  testthat::expect_equal(actual$germline_alignment, rep("AA..NGG.TT", 4L))
})


testthat::test_that("first non-empty V and J alleles support comma and pipe calls", {
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, "leading-separators", "calls", "IGH", "AAAAAAAAAA", ",IGHV1*01,IGHV2*01", "|IGHJ1*01|IGHJ2*01", 1L, 4L, 1L, 4L
  )

  actual <- annotate_ig_germlines(
    make_ig_germline_idata(rows),
    references = germline_references(), species = "Homo_sapiens"
  ) |>
    collect_germlines()

  testthat::expect_equal(actual$germline_alignment, "AA..NNGG.T")
})


testthat::test_that("groups without a modal V/J co-occurrence retain NA germlines", {
  # This 2x2-like layout has unique modal V1 (four rows) and J1 (four rows),
  # but no row contains V1/J1: a constructively possible missing intersection.
  rows <- tibble::tribble(
    ~row_id, ~sequence_id, ~clone_id, ~locus, ~sequence_alignment, ~v_call, ~j_call, ~v_germline_start, ~v_germline_end, ~j_germline_start, ~j_germline_end,
    1L, "a", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ2*01", 1L, 4L, 1L, 4L,
    2L, "b", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ2*01", 1L, 4L, 1L, 4L,
    3L, "c", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ3*01", 1L, 4L, 1L, 4L,
    4L, "d", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV1*01", "IGHJ3*01", 1L, 4L, 1L, 4L,
    5L, "e", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV2*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    6L, "f", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV2*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    7L, "g", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV3*01", "IGHJ1*01", 1L, 4L, 1L, 4L,
    8L, "h", "no-intersection", "IGH", "AAAAAAAAAA", "IGHV3*01", "IGHJ1*01", 1L, 4L, 1L, 4L
  )
  refs <- dplyr::bind_rows(
    germline_references(),
    tibble::tribble(
      ~species, ~locus, ~segment, ~allele_id, ~sequence,
      "Homo_sapiens", "IGH", "V", "IGHV3*01", "GG..CC",
      "Homo_sapiens", "IGH", "J", "IGHJ3*01", "TT.CC"
    )
  )

  out <- annotate_ig_germlines(
    make_ig_germline_idata(rows), references = refs, species = "Homo_sapiens"
  )
  actual <- collect_germlines(out)

  testthat::expect_equal(nrow(actual), 8L)
  testthat::expect_true(all(is.na(actual$germline_alignment)))
})
