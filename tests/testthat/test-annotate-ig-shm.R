shm_regions <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")


make_shm_rows <- function(sequence, germline = sequence, cell_id = seq_along(sequence)) {
  fixture_id <- cell_id
  rows <- tibble::tibble(
    sample = "S1",
    cell_id = paste0("cell", fixture_id),
    umi = 1L,
    sequence_alignment = sequence,
    germline_alignment = germline,
    germline_alignment_d_mask = germline,
    v_sequence_start = 1L,
    v_alignment_start = 1L,
    existing_annotation = paste0("keep", fixture_id)
  )
  for (i in seq_along(shm_regions)) {
    rows[[paste0(shm_regions[[i]], "_start")]] <- 3L * i - 2L
    rows[[paste0(shm_regions[[i]], "_end")]] <- 3L * i
  }
  rows
}


only_fwr1 <- function(rows, start = 1L, end = nchar(rows$sequence_alignment)) {
  for (region in setdiff(shm_regions, "fwr1")) {
    rows[[paste0(region, "_start")]] <- NA_integer_
    rows[[paste0(region, "_end")]] <- NA_integer_
  }
  rows$fwr1_start <- start
  rows$fwr1_end <- end
  rows
}


make_shm_idata <- function(rows, sequence_col = "sequence_alignment") {
  input <- tempfile(fileext = ".tsv")
  readr::write_tsv(rows, input, na = "")
  immundata::read_repertoires(
    path = input,
    schema = immundata::make_receptor_schema(features = sequence_col),
    barcode_col = "cell_id",
    umi_col = "umi",
    repertoire_schema = "sample",
    rename_columns = NULL,
    preprocess = NULL,
    postprocess = NULL,
    output_folder = tempfile("ig-shm-"),
    verbose = FALSE
  )
}


collect_shm <- function(idata) {
  idata$annotations |>
    dplyr::collect() |>
    dplyr::arrange(cell_id)
}


shm_output_cols <- c(
  "mu_count", "mu_freq",
  paste0(
    c("mu_count", "mu_count", "mu_freq", "mu_freq"),
    rep(c("", "_cdr", "_fwr"), each = 4L),
    rep(c("_r", "_s", "_r", "_s"), 3L)
  ),
  unlist(lapply(shm_regions, function(region) {
    paste0(
      c("mu_count_", "mu_count_", "mu_freq_", "mu_freq_"),
      region,
      c("_r", "_s", "_r", "_s")
    )
  }), use.names = FALSE)
)
shm_group_output_cols <- c("mu_freq_median", "mu_freq_mean", "mu_freq_n")


testthat::test_that("SHM annotation preserves ImmunData rows, metadata, schema, and laziness", {
  rows <- make_shm_rows(c("GCTGCTGCTGCTGCTGCTGCT", "GCCGCTGCTGCTGCTGCTGCT"))
  idata <- make_shm_idata(rows)
  before_schema <- idata$schema_receptor

  out <- annotate_ig_shm(idata)
  actual <- collect_shm(out)

  testthat::expect_s3_class(out, "ImmunData")
  testthat::expect_s3_class(out$annotations, "duckplyr_df")
  testthat::expect_equal(nrow(actual), 2L)
  testthat::expect_equal(actual$existing_annotation, c("keep1", "keep2"))
  testthat::expect_equal(out$schema_receptor, before_schema)
  testthat::expect_true(all(shm_output_cols %in% names(out$annotations)))
})


testthat::test_that("SHM has no unrelated identifier dependency", {
  rows <- make_shm_rows("GCTGCTGCTGCTGCTGCTGCT")

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count, 0L)
  testthat::expect_false(any(grepl("_(heavy|light|all)$", names(actual))))
})


testthat::test_that("automatic clone grouping is locus scoped", {
  withr::local_options(rlib_message_verbosity = "verbose")
  germline <- paste(rep("GCT", 7L), collapse = "")
  mutation_counts <- c(0L, 1L, 2L, 7L)
  sequences <- vapply(mutation_counts, function(n) {
    paste0(c(rep("ACT", n), rep("GCT", 7L - n)), collapse = "")
  }, character(1))
  rows <- make_shm_rows(sequences, rep(germline, 4L))
  rows$clone_id <- "family1"
  rows$locus <- c("IGH", "IGH", "IGK", "IGK")
  idata <- make_shm_idata(rows)

  testthat::expect_message(
    out <- annotate_ig_shm(idata),
    "resolved to.*clone_id.*locus"
  )
  actual <- collect_shm(out)

  testthat::expect_equal(actual$mu_freq_median, c(1 / 42, 1 / 42, 3 / 14, 3 / 14))
  testthat::expect_equal(actual$mu_freq_mean, c(1 / 42, 1 / 42, 3 / 14, 3 / 14))
  testthat::expect_equal(actual$mu_freq_n, c(2, 2, 2, 2))

  testthat::expect_no_message(
    pooled <- annotate_ig_shm(idata, by = "clone_id")
  )
  pooled <- collect_shm(pooled)
  testthat::expect_equal(pooled$mu_freq_median, rep(1 / 14, 4L))
  testthat::expect_equal(pooled$mu_freq_mean, rep(5 / 42, 4L))
  testthat::expect_equal(pooled$mu_freq_n, rep(4, 4L))
})


testthat::test_that("automatic external clone grouping uses subject, clone, and locus", {
  withr::local_options(rlib_message_verbosity = "verbose")
  germline <- paste(rep("GCT", 7L), collapse = "")
  sequences <- c(
    germline,
    paste0("ACT", substr(germline, 4L, nchar(germline))),
    paste0("ACTACT", substr(germline, 7L, nchar(germline)))
  )
  rows <- make_shm_rows(sequences, rep(germline, 3L))
  rows$subject_id <- c("subject1", "subject1", "subject2")
  rows$clone_id <- "clone1"
  rows$locus <- "IGH"
  idata <- make_shm_idata(rows)

  testthat::expect_message(
    out <- annotate_ig_shm(idata),
    "resolved to.*subject_id.*clone_id.*locus"
  )
  actual <- collect_shm(out)

  testthat::expect_equal(actual$mu_freq_median, c(1 / 42, 1 / 42, 2 / 21))
  testthat::expect_equal(actual$mu_freq_mean, c(1 / 42, 1 / 42, 2 / 21))
  testthat::expect_equal(actual$mu_freq_n, c(2, 2, 1))
})


testthat::test_that("NULL disables grouping and unresolved auto skips it", {
  withr::local_options(rlib_message_verbosity = "verbose")
  rows <- make_shm_rows("GCTGCTGCTGCTGCTGCTGCT")
  rows$clone_id <- "family1"
  idata <- make_shm_idata(rows)

  testthat::expect_no_message(out <- annotate_ig_shm(idata, by = NULL))
  testthat::expect_false(any(shm_group_output_cols %in% names(out$annotations)))

  rows$clone_id <- NULL
  idata <- make_shm_idata(rows)
  testthat::expect_message(
    out <- annotate_ig_shm(idata),
    "found no supported grouping columns"
  )
  testthat::expect_false(any(shm_group_output_cols %in% names(out$annotations)))
})


testthat::test_that("missing grouping keys and all-missing groups do not pool", {
  rows <- make_shm_rows(rep("GCTGCTGCTGCTGCTGCTGCT", 3L))
  rows$clone_id <- c(NA_character_, NA_character_, "family1")
  rows$locus <- "IGH"
  for (region in shm_regions) {
    rows[[paste0(region, "_start")]][[3]] <- NA_integer_
    rows[[paste0(region, "_end")]][[3]] <- NA_integer_
  }

  actual <- suppressMessages(annotate_ig_shm(make_shm_idata(rows))) |>
    collect_shm()

  testthat::expect_true(all(is.na(actual$mu_freq_median)))
  testthat::expect_true(all(is.na(actual$mu_freq_mean)))
  testthat::expect_equal(actual$mu_freq_n, c(0, 0, 0))
})


testthat::test_that("custom sequence and germline columns are supported", {
  rows <- make_shm_rows(
    "GCCGCTGCTGCTGCTGCTGCT",
    "GCTGCTGCTGCTGCTGCTGCT"
  ) |>
    dplyr::rename(observed = sequence_alignment, reference = germline_alignment_d_mask)
  idata <- make_shm_idata(rows, sequence_col = "observed")

  actual <- annotate_ig_shm(
    idata,
    sequence_col = "observed",
    germline_col = "reference"
  ) |>
    collect_shm()

  testthat::expect_equal(actual$mu_count, 1L)
  testthat::expect_equal(actual$mu_count_fwr1_s, 1L)
})


testthat::test_that("SHM defaults to the reconstructed germline alignment", {
  rows <- make_shm_rows(
    "GCCGCTGCTGCTGCTGCTGCT",
    "GCTGCTGCTGCTGCTGCTGCT"
  )
  rows$germline_alignment <- rows$germline_alignment_d_mask
  rows$germline_alignment_d_mask <- rows$sequence_alignment

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count, 1L)
})


testthat::test_that("missing canonical coordinates use standard annotation validation", {
  rows <- make_shm_rows("GCTGCTGCTGCTGCTGCTGCT") |>
    dplyr::select(-cdr2_end)
  idata <- make_shm_idata(rows)

  testthat::expect_error(
    annotate_ig_shm(idata),
    "Passed column name.*cdr2_end.*not in the input ImmunData"
  )
})


testthat::test_that("rerunning replaces SHM columns without duplication", {
  rows <- make_shm_rows("GCCGCTGCTGCTGCTGCTGCT")
  rows$clone_id <- "family1"
  rows$locus <- "IGH"
  out <- suppressMessages(
    make_shm_idata(rows) |>
      annotate_ig_shm() |>
      annotate_ig_shm()
  )

  output_cols <- c(shm_output_cols, shm_group_output_cols)
  testthat::expect_true(all(vapply(output_cols, function(x) {
    sum(names(out$annotations) == x) == 1L
  }, logical(1))))
})


testthat::test_that("large diverse file-backed input remains DuckDB-backed", {
  n_rows <- 10000L
  germline <- paste(rep("GCT", 120L), collapse = "")
  bases <- c("A", "C", "G", "T")
  encode_suffix <- function(index) {
    digits <- integer(7L)
    for (position in seq_along(digits)) {
      digits[[position]] <- index %% 4L
      index <- index %/% 4L
    }
    paste0(bases[digits + 1L], collapse = "")
  }
  suffixes <- vapply(0:(n_rows - 1L), encode_suffix, character(1))
  sequences <- paste0(substr(germline, 1L, 353L), suffixes)
  rows <- make_shm_rows(
    sequences,
    rep(germline, n_rows),
    cell_id = seq_len(n_rows)
  )

  out <- make_shm_idata(rows) |>
    annotate_ig_shm()

  testthat::expect_s3_class(out$annotations, "duckplyr_df")
  summary <- out$annotations |>
    dplyr::summarise(
      n = dplyr::n(),
      annotated = sum(!is.na(mu_count))
    ) |>
    dplyr::collect()
  testthat::expect_equal(summary$n, n_rows)
  testthat::expect_equal(summary$annotated, n_rows)
})


testthat::test_that("overall mutation eligibility follows the alignment contract", {
  cases <- tibble::tribble(
    ~cell_id, ~sequence, ~germline, ~count, ~denominator,
    1L, "GCT", "GCT", 0L, 3L,
    2L, "GCC", "GCT", 1L, 3L, # GCT -> GCC, synonymous
    3L, "ACT", "GCT", 1L, 3L, # GCT -> ACT, replacement
    4L, "ANN", "GCT", 0L, 1L,
    5L, "NCT", "GCT", 0L, 2L,
    6L, "GCT", "GNT", 0L, 2L,
    7L, "G.T", "GCT", 0L, 2L,
    8L, "GCT", "G-T", 0L, 2L,
    9L, "NNN", "NNN", 0L, 0L
  )
  rows <- make_shm_rows(cases$sequence, cases$germline, cases$cell_id) |>
    only_fwr1(start = 1L, end = 3L)

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count, cases$count)
  expected_freq <- ifelse(cases$denominator == 0L, NA_real_, cases$count / cases$denominator)
  testthat::expect_equal(actual$mu_freq, expected_freq)
  testthat::expect_false(any(is.nan(actual$mu_freq)))
})


testthat::test_that("missing alignments make every SHM output missing", {
  rows <- make_shm_rows(
    c(NA_character_, "GCTGCTGCTGCTGCTGCTGCT"),
    c("GCTGCTGCTGCTGCTGCTGCT", NA_character_)
  )
  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_true(all(is.na(actual[shm_output_cols])))
})


testthat::test_that("mutations are classified independently against the germline codon", {
  # GCT -> GCC: Ala -> Ala (silent)
  # GCT -> ACT: Ala -> Thr (replacement)
  # TGG -> TAG: Trp -> stop (unclassified)
  # AAA -> GGA: A1G gives GAA (Glu), A2G gives AGA (Arg): two replacements
  cases <- tibble::tribble(
    ~cell_id, ~sequence, ~germline, ~r, ~s, ~total,
    1L, "GCC", "GCT", 0L, 1L, 1L,
    2L, "ACT", "GCT", 1L, 0L, 1L,
    3L, "TAG", "TGG", 0L, 0L, 0L,
    4L, "GGA", "AAA", 2L, 0L, 2L,
    5L, "GNC", "GCT", 0L, 0L, 0L,
    6L, "GC-", "GCT", 0L, 0L, 0L,
    7L, "G-CC", "G-CT", 0L, 0L, 0L
  )
  rows <- make_shm_rows(cases$sequence, cases$germline, cases$cell_id) |>
    only_fwr1(start = 1L, end = 3L)
  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count_fwr1_r, cases$r)
  testthat::expect_equal(actual$mu_count_fwr1_s, cases$s)
  testthat::expect_equal(actual$mu_count, cases$total)
  testthat::expect_equal(actual$mu_freq[actual$cell_id == "cell3"], 0)
})


testthat::test_that("query coordinates map through alignment starts, N, dots, and gaps", {
  rows <- make_shm_rows(
    c("ACT", "ACT", "TTTACT", "TTTACT", "NCTACT", "-GCTACT", ".GCTACT", "GCTACT"),
    c("GCT", "GCT", "TTTGCT", "TTTGCT", "GCTGCT", "-GCTGCT", ".GCTGCT", "-CTGCT")
  ) |>
    only_fwr1()
  rows$cell_id <- paste0("cell", 1:8)
  rows$v_sequence_start <- c(1L, 4L, 1L, 4L, 1L, 1L, 1L, 1L)
  rows$v_alignment_start <- c(1L, 1L, 4L, 4L, 1L, 1L, 1L, 1L)
  rows$fwr1_start <- c(1L, 4L, 1L, 4L, 4L, 4L, 4L, 4L)
  rows$fwr1_end <- c(3L, 6L, 3L, 6L, 6L, 6L, 6L, 6L)

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count, rep(1L, 8L))
  testthat::expect_equal(actual$mu_freq, rep(1 / 3, 8L))
})


testthat::test_that("AIRR region ends are inclusive and adjacent regions remain distinct", {
  rows <- make_shm_rows(
    "GCCACTGCTGCTGCTGCTGCT",
    "GCTGCTGCTGCTGCTGCTGCT"
  )
  rows$fwr1_start <- 1L
  rows$fwr1_end <- 3L
  rows$cdr1_start <- 4L
  rows$cdr1_end <- 6L

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count_fwr1_s, 1L)
  testthat::expect_equal(actual$mu_count_cdr1_r, 1L)
  testthat::expect_equal(actual$mu_count, 2L)
})


testthat::test_that("each of seven regions reports replacement and silent mutations", {
  germline <- paste(rep("GCT", 7L), collapse = "")
  replacement <- paste(rep("ACT", 7L), collapse = "")
  silent <- paste(rep("GCC", 7L), collapse = "")
  actual <- make_shm_idata(make_shm_rows(c(replacement, silent), c(germline, germline))) |>
    annotate_ig_shm() |>
    collect_shm()

  for (region in shm_regions) {
    testthat::expect_equal(actual[[paste0("mu_count_", region, "_r")]], c(1L, 0L))
    testthat::expect_equal(actual[[paste0("mu_count_", region, "_s")]], c(0L, 1L))
    testthat::expect_equal(actual[[paste0("mu_freq_", region, "_r")]], c(1 / 3, 0))
    testthat::expect_equal(actual[[paste0("mu_freq_", region, "_s")]], c(0, 1 / 3))
  }
  testthat::expect_equal(actual$mu_count, c(7L, 7L))
  regional_count_cols <- unlist(lapply(shm_regions, function(region) {
    paste0("mu_count_", region, c("_r", "_s"))
  }), use.names = FALSE)
  regional_sum <- rowSums(actual[regional_count_cols])
  testthat::expect_equal(actual$mu_count, regional_sum)
})


testthat::test_that("missing region boundaries are local and overall uses available regions", {
  rows <- make_shm_rows(c(
    "ACTACTACTACTACTACTACT",
    "ACTACTACTACTACTACTACT"
  ), germline = rep(paste(rep("GCT", 7L), collapse = ""), 2L))
  rows$fwr4_start[[1]] <- NA_integer_
  rows$fwr4_end[[1]] <- NA_integer_
  rows$cdr2_start[[2]] <- NA_integer_

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_true(all(is.na(actual[1, grep("fwr4", names(actual))])))
  testthat::expect_equal(actual$mu_count[[1]], 6L)
  testthat::expect_true(all(is.na(actual[2, grep("^mu_.*cdr2", names(actual))])))
  testthat::expect_equal(actual$mu_count_cdr1_r[[2]], 1L)
  testthat::expect_equal(actual$mu_count_fwr3_r[[2]], 1L)
  testthat::expect_equal(actual$mu_count[[2]], 6L)
})


testthat::test_that("combined CDR metrics use available regions", {
  germline <- paste(rep("GCT", 7L), collapse = "")
  sequence <- paste(rep("ACT", 7L), collapse = "")
  rows <- make_shm_rows(
    rep(sequence, 2L),
    rep(germline, 2L)
  )
  rows$cdr3_start[[1]] <- NA_integer_
  rows$cdr3_end[[1]] <- NA_integer_
  for (region in paste0("cdr", 1:3)) {
    rows[[paste0(region, "_start")]][[2]] <- NA_integer_
    rows[[paste0(region, "_end")]][[2]] <- NA_integer_
  }

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm(by = NULL) |>
    collect_shm()

  cdr3_cols <- paste0(
    c("mu_count_", "mu_count_", "mu_freq_", "mu_freq_"),
    "cdr3",
    c("_r", "_s", "_r", "_s")
  )
  combined_cdr_cols <- c(
    "mu_count_cdr_r", "mu_count_cdr_s",
    "mu_freq_cdr_r", "mu_freq_cdr_s"
  )
  testthat::expect_true(all(is.na(actual[1, cdr3_cols])))
  testthat::expect_equal(
    unlist(actual[1, combined_cdr_cols], use.names = FALSE),
    c(2, 0, 1 / 3, 0)
  )
  testthat::expect_true(all(is.na(actual[2, combined_cdr_cols])))
})


testthat::test_that("all SHM outputs are missing when every region is unavailable", {
  rows <- make_shm_rows("GCTGCTGCTGCTGCTGCTGCT")
  for (region in shm_regions) {
    rows[[paste0(region, "_start")]] <- NA_integer_
    rows[[paste0(region, "_end")]] <- NA_integer_
  }

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_true(all(is.na(actual[shm_output_cols])))
})


testthat::test_that("D-masked CDR3 positions are excluded from its denominator", {
  rows <- make_shm_rows("GCCACT", "GCTNNN")
  for (region in setdiff(shm_regions, "cdr3")) {
    rows[[paste0(region, "_start")]] <- NA_integer_
    rows[[paste0(region, "_end")]] <- NA_integer_
  }
  rows$cdr3_start <- 1L
  rows$cdr3_end <- 6L

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  testthat::expect_equal(actual$mu_count_cdr3_r, 0L)
  testthat::expect_equal(actual$mu_count_cdr3_s, 1L)
  testthat::expect_equal(actual$mu_freq_cdr3_r, 0)
  testthat::expect_equal(actual$mu_freq_cdr3_s, 1 / 3)
})


testthat::test_that("locus labels do not alter per-row SHM", {
  germline <- paste(rep("GCT", 7L), collapse = "")
  sequence <- paste(c("ACT", "GCC", rep("GCT", 5L)), collapse = "")
  rows <- make_shm_rows(
    rep(sequence, 3L),
    rep(germline, 3L),
    cell_id = 1:3
  )
  rows$locus <- c("IGH", "IGK", "arbitrary")
  actual_idata <- make_shm_idata(rows) |>
    annotate_ig_shm()
  actual <- collect_shm(actual_idata)

  testthat::expect_equal(actual$locus, c("IGH", "IGK", "arbitrary"))
  testthat::expect_equal(actual$mu_count, rep(2L, 3L))
  testthat::expect_true(all(vapply(shm_output_cols, function(column) {
    length(unique(actual[[column]])) == 1L
  }, logical(1))))
})


testthat::test_that("full IMGT fixture matches values checked with SHazaM", {
  # Independently checked against shazam::observedMutations() 1.3.2 using a
  # custom seven-level RegionDefinition over the real IMGT VDJ boundaries for
  # a 360-nt alignment with junction_length 45. Expected values are hard-coded
  # so SHazaM is not a test or runtime dependency.
  germline <- paste(rep("GCT", 120L), collapse = "")
  sequence <- strsplit(germline, "", fixed = TRUE)[[1]]
  sequence[c(1L, 115L, 196L, 197L, 352L)] <- "A"
  sequence[c(81L, 168L, 315L)] <- "C"
  sequence <- paste(sequence, collapse = "")

  rows <- make_shm_rows(sequence, germline)
  starts <- c(1L, 79L, 115L, 166L, 196L, 313L, 352L)
  ends <- c(78L, 114L, 165L, 195L, 312L, 351L, 360L)
  for (i in seq_along(shm_regions)) {
    rows[[paste0(shm_regions[[i]], "_start")]] <- starts[[i]]
    rows[[paste0(shm_regions[[i]], "_end")]] <- ends[[i]]
  }

  actual <- make_shm_idata(rows) |>
    annotate_ig_shm() |>
    collect_shm()

  expected_r <- c(fwr1 = 1L, cdr1 = 0L, fwr2 = 1L, cdr2 = 0L,
    fwr3 = 2L, cdr3 = 0L, fwr4 = 1L)
  expected_s <- c(fwr1 = 0L, cdr1 = 1L, fwr2 = 0L, cdr2 = 1L,
    fwr3 = 0L, cdr3 = 1L, fwr4 = 0L)
  region_lengths <- stats::setNames(ends - starts + 1L, shm_regions)

  testthat::expect_equal(actual$mu_count, 8L)
  testthat::expect_equal(actual$mu_freq, 8 / 360)
  testthat::expect_equal(actual$mu_count_r, 5L)
  testthat::expect_equal(actual$mu_count_s, 3L)
  testthat::expect_equal(actual$mu_freq_r, 5 / 360)
  testthat::expect_equal(actual$mu_freq_s, 3 / 360)
  testthat::expect_equal(actual$mu_count_cdr_r, 0L)
  testthat::expect_equal(actual$mu_count_cdr_s, 3L)
  testthat::expect_equal(actual$mu_freq_cdr_r, 0)
  testthat::expect_equal(actual$mu_freq_cdr_s, 3 / 105)
  testthat::expect_equal(actual$mu_count_fwr_r, 5L)
  testthat::expect_equal(actual$mu_count_fwr_s, 0L)
  testthat::expect_equal(actual$mu_freq_fwr_r, 5 / 255)
  testthat::expect_equal(actual$mu_freq_fwr_s, 0)
  for (region in shm_regions) {
    testthat::expect_equal(
      actual[[paste0("mu_count_", region, "_r")]], expected_r[[region]]
    )
    testthat::expect_equal(
      actual[[paste0("mu_count_", region, "_s")]], expected_s[[region]]
    )
    testthat::expect_equal(
      actual[[paste0("mu_freq_", region, "_r")]],
      expected_r[[region]] / region_lengths[[region]]
    )
    testthat::expect_equal(
      actual[[paste0("mu_freq_", region, "_s")]],
      expected_s[[region]] / region_lengths[[region]]
    )
  }
})
