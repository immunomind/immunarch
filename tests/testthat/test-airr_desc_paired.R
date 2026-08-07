make_paired_desc_idata <- function() {
  cells <- tibble::tribble(
    ~sample, ~cell_id, ~locus, ~cdr3_aa, ~v_call, ~umi_count,
    "S1", "S1-1", "TRA", "AAA",    "TRAV1", 10,
    "S1", "S1-1", "TRB", "BBBB",   "TRBV1", 10,
    "S1", "S1-2", "TRA", "AAA",    "TRAV1", 10,
    "S1", "S1-2", "TRB", "BBBB",   "TRBV1", 10,
    "S1", "S1-3", "TRA", "CCCCC",  "TRAV2", 10,
    "S1", "S1-3", "TRB", "DDDDDD", "TRBV2", 10,
    "S2", "S2-1", "TRA", "AAA",    "TRAV1", 10,
    "S2", "S2-1", "TRB", "BBBB",   "TRBV1", 10,
    "S2", "S2-2", "TRA", "EEE",    "TRAV3", 10,
    "S2", "S2-2", "TRB", "FFFFF",  "TRBV3", 10,
    "S2", "S2-3", "TRA", "EEE",    "TRAV3", 10,
    "S2", "S2-3", "TRB", "FFFFF",  "TRBV3", 10
  ) |>
    duckplyr::as_duckdb_tibble()

  receptor_schema <- immundata::make_receptor_schema(
    features = c("cdr3_aa", "v_call"),
    chains = c("TRA", "TRB")
  )

  annotations <- immundata::agg_receptors(
    dataset = cells,
    schema = receptor_schema,
    barcode_col = "cell_id",
    locus_col = "locus",
    umi_col = "umi_count"
  )

  immundata::ImmunData$new(
    schema = receptor_schema,
    annotations = annotations
  ) |>
    immundata::agg_repertoires("sample")
}


test_that("paired-chain descriptive statistics return automatic locus groups", {
  idata <- make_paired_desc_idata()

  genes <- airr_desc_genes(idata, autojoin = FALSE)
  lengths <- airr_desc_lengths(idata, autojoin = FALSE)

  expect_named(
    genes,
    c("v_call", "imd_repertoire_id", "locus", "n"),
    ignore.order = FALSE
  )
  expect_named(
    lengths,
    c("imd_repertoire_id", "locus", "seq_len", "n", "prop", "pct"),
    ignore.order = FALSE
  )
  expect_setequal(unique(genes$locus), c("TRA", "TRB"))
  expect_setequal(unique(lengths$locus), c("TRA", "TRB"))

  prop_sums <- lengths |>
    dplyr::summarise(
      prop = sum(.data$prop),
      pct = sum(.data$pct),
      .by = c("imd_repertoire_id", "locus")
    )

  expect_equal(prop_sums$prop, rep(1, nrow(prop_sums)), tolerance = 1e-12)
  expect_equal(prop_sums$pct, rep(100, nrow(prop_sums)), tolerance = 1e-10)
})


test_that("paired-chain gene usage preserves both chains at both levels", {
  idata <- make_paired_desc_idata()

  receptor_genes <- airr_desc_genes(
    idata,
    level = "receptor",
    autojoin = FALSE
  ) |>
    dplyr::filter(.data$imd_repertoire_id == 1L) |>
    dplyr::arrange(.data$locus, .data$v_call)

  expect_equal(receptor_genes$locus, c("TRA", "TRA", "TRB", "TRB"))
  expect_equal(receptor_genes$v_call, c("TRAV1", "TRAV2", "TRBV1", "TRBV2"))
  expect_equal(receptor_genes$n, c(1, 1, 1, 1))

  barcode_genes <- airr_desc_genes(
    idata,
    level = "barcode",
    autojoin = FALSE
  ) |>
    dplyr::filter(.data$imd_repertoire_id == 1L) |>
    dplyr::arrange(.data$locus, .data$v_call)

  expect_equal(barcode_genes$locus, c("TRA", "TRA", "TRB", "TRB"))
  expect_equal(barcode_genes$v_call, c("TRAV1", "TRAV2", "TRBV1", "TRBV2"))
  expect_equal(barcode_genes$n, c(2, 1, 2, 1))

  pooled_genes <- airr_desc_genes(idata, by = NULL, autojoin = FALSE)
  expect_false("locus" %in% colnames(pooled_genes))
  expect_setequal(
    pooled_genes |>
      dplyr::filter(.data$imd_repertoire_id == 1L) |>
      dplyr::pull(.data$v_call),
    c("TRAV1", "TRAV2", "TRBV1", "TRBV2")
  )
})


test_that("paired-chain length usage preserves annotation-row counting", {
  idata <- make_paired_desc_idata()

  lengths <- airr_desc_lengths(idata, autojoin = FALSE) |>
    dplyr::filter(.data$imd_repertoire_id == 1L) |>
    dplyr::arrange(.data$locus, .data$seq_len)

  expect_equal(lengths$locus, c("TRA", "TRA", "TRB", "TRB"))
  expect_equal(lengths$seq_len, c(3, 5, 4, 6))
  expect_equal(lengths$n, c(2, 1, 2, 1))
  expect_equal(
    lengths$prop,
    c(2 / 3, 1 / 3, 2 / 3, 1 / 3),
    tolerance = 1e-12
  )

  pooled_lengths <- airr_desc_lengths(idata, by = NULL, autojoin = FALSE)
  expect_false("locus" %in% colnames(pooled_lengths))
})


test_that("descriptive grouping returns requested annotation columns", {
  idata <- make_paired_desc_idata()

  genes <- airr_desc_genes(
    idata,
    by = c("locus", "sample"),
    autojoin = FALSE
  )
  lengths <- airr_desc_lengths(
    idata,
    by = c("locus", "sample"),
    autojoin = FALSE
  )

  expect_true(all(c("locus", "sample") %in% colnames(genes)))
  expect_true(all(c("locus", "sample") %in% colnames(lengths)))

  expect_error(
    airr_desc_genes(idata, by = "missing_group"),
    "missing_group"
  )
  expect_error(
    airr_desc_lengths(idata, by = "missing_group"),
    "missing_group"
  )
})


test_that("descriptive outputs preserve names and repertoire mapping on Parquet", {
  idata <- make_paired_desc_idata()
  snapshot_dir <- tempfile("paired-desc-")
  on.exit(unlink(snapshot_dir, recursive = TRUE), add = TRUE)

  immundata::write_immundata(idata, snapshot_dir)
  parquet_idata <- immundata::read_immundata(snapshot_dir, verbose = FALSE)

  genes <- airr_desc_genes(parquet_idata)
  lengths <- airr_desc_lengths(parquet_idata)

  expect_named(
    genes,
    c("v_call", "imd_repertoire_id", "locus", "n", "sample"),
    ignore.order = FALSE
  )
  expect_named(
    lengths,
    c(
      "imd_repertoire_id", "locus", "seq_len", "n", "prop", "pct",
      "sample"
    ),
    ignore.order = FALSE
  )
  expect_false(any(grepl("\\.[xy]$", colnames(genes))))
  expect_false(any(grepl("\\.[xy]$", colnames(lengths))))

  expected_rep_map <- parquet_idata$repertoires |>
    dplyr::select(dplyr::all_of(c("imd_repertoire_id", "sample")))

  expect_equal(
    genes |>
      dplyr::distinct(.data$imd_repertoire_id, .data$sample) |>
      dplyr::arrange(.data$imd_repertoire_id) |>
      as.data.frame(),
    expected_rep_map |>
      dplyr::arrange(.data$imd_repertoire_id) |>
      as.data.frame()
  )
  expect_equal(
    lengths |>
      dplyr::distinct(.data$imd_repertoire_id, .data$sample) |>
      dplyr::arrange(.data$imd_repertoire_id) |>
      as.data.frame(),
    expected_rep_map |>
      dplyr::arrange(.data$imd_repertoire_id) |>
      as.data.frame()
  )
})
