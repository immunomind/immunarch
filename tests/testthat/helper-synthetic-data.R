# Import tracking fixtures through the public reader, including bulk weights
# and paired chains. Keep backing files until the test run has finished.
make_track_idata <- function(rows = NULL, paired = FALSE,
                             repertoires = "sample", features = c("junction", "v_call")) {
  if (is.null(rows)) {
    rows <- tibble::tibble(
      sample = c(rep("S1", 4), rep("S2", 4), rep("S3", 3), rep("S4", 2), rep("S5", 2)),
      junction = c("AAA", "CCC", "GGG", "AAA", "AAA", "CCC", "TTT", "GGG",
                   "AAA", "GGG", "ACG", "CCC", "ACG", "GGA", "ACG"),
      v_call = c("V1", "V2", "V3", "V2", "V1", "V2", "V4", "V3",
                 "V1", "V3", "V5", "V2", "V5", "V6", "V5"),
      abundance = c(8L, 2L, 1L, 1L, 1L, 9L, 4L, 1L, 2L, 3L, 5L, 6L, 1L, 7L, 1L),
      trait = c("other", "hit", rep("other", 13))
    )
    rows$group <- ifelse(rows$sample %in% c("S1", "S3", "S5"), "Case", "Control")
    rows$row_id <- seq_len(nrow(rows))
  }
  input <- tempfile(fileext = ".tsv")
  output <- tempfile("track-test-")
  withr::defer(unlink(c(input, output), recursive = TRUE),
               envir = testthat::teardown_env())
  readr::write_tsv(rows, input)
  immundata::read_repertoires(
    input,
    schema = immundata::make_receptor_schema(
      features = features, chains = if (paired) c("IGH", "IGL") else NULL
    ),
    count_col = if (paired) NULL else "abundance",
    barcode_col = if (paired) "cell" else NULL,
    umi_col = if (paired) "umi" else NULL,
    locus_col = if (paired) "locus" else NULL,
    repertoire_schema = repertoires, output_folder = output,
    preprocess = NULL, postprocess = NULL, rename_columns = NULL,
    verbose = FALSE
  )
}


make_public_idata <- function() {
  rows <- tibble::tibble(
    sample = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S4"),
    junction = c("AAA", "CCC", "ACG", "AAA", "CCC", "AAA", "GGG", "TTT"),
    v_call = "V1", abundance = c(5L, 1L, 1L, 10L, 4L, 20L, 8L, 2L),
    group = c("Case", "Case", "Case", "Case", "Case", "Control", "Control", "Control"),
    row_id = seq_len(8L)
  )
  make_track_idata(rows, repertoires = c("sample", "group")) |>
    immundata::agg_strata("group")
}


make_test_repertoire_idata <- function(
  receptors,
  repertoires,
  counts = NULL,
  proportions = NULL,
  all_repertoires = unique(repertoires)
) {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  proportion_col <- immundata::imd_schema("proportion")

  annotations <- tibble::tibble(
    !!receptor_col := receptors,
    !!repertoire_col := repertoires,
    cdr3_aa = receptors
  )
  if (!is.null(counts)) {
    annotations[[count_col]] <- counts
  }
  if (!is.null(proportions)) {
    annotations[[proportion_col]] <- proportions
  }

  repertoire_table <- tibble::tibble(
    !!repertoire_col := all_repertoires,
    Group = all_repertoires
  )

  idata <- immundata::ImmunData$new(
    schema = "cdr3_aa",
    annotations = duckplyr::as_duckdb_tibble(
      annotations,
      prudence = "stingy"
    ),
    repertoires = duckplyr::as_duckdb_tibble(
      repertoire_table,
      prudence = "stingy"
    )
  )
  idata$schema_repertoire <- "Group"
  idata
}


make_repsim_idata <- function(
  receptor_counts,
  all_repertoires = unique(receptor_counts$repertoire),
  normalize = FALSE
) {
  proportions <- NULL
  if ("proportion" %in% names(receptor_counts)) {
    proportions <- receptor_counts$proportion
  } else if (normalize) {
    proportions <- receptor_counts |>
      dplyr::mutate(
        proportion = .data$count / sum(.data$count),
        .by = "repertoire"
      ) |>
      dplyr::pull(.data$proportion)
  }

  make_test_repertoire_idata(
    receptors = receptor_counts$receptor,
    repertoires = receptor_counts$repertoire,
    counts = receptor_counts$count,
    proportions = proportions,
    all_repertoires = all_repertoires
  )
}


make_repsim_set_idata <- function(duplicated = FALSE) {
  receptor_counts <- if (duplicated) {
    tibble::tribble(
      ~receptor, ~repertoire, ~count,
      "r1", "R1", 1,
      "r1", "R1", 1,
      "r2", "R1", 1,
      "r1", "R2", 1,
      "r2", "R2", 1,
      "r2", "R2", 1
    )
  } else {
    tibble::tribble(
      ~receptor, ~repertoire, ~count,
      "r1", "R1", 1,
      "r3", "R1", 1,
      "r2", "R2", 1,
      "r3", "R2", 1,
      "r1", "R3", 1,
      "r2", "R3", 1,
      "r3", "R3", 1
    )
  }

  make_repsim_idata(receptor_counts)
}


# Import the package example once, but isolate each test from R6 mutations.
.test_idata_cache <- new.env(parent = emptyenv())


make_test_idata <- function() {
  if (is.null(.test_idata_cache$base)) {
    .test_idata_cache$base <- get_test_idata()
  }

  .test_idata_cache$base$clone(deep = TRUE)
}


make_aggregated_test_idata <- function(
  schema = c("Response", "Therapy")
) {
  cache_key <- paste(schema, collapse = "\r")
  if (is.null(.test_idata_cache[[cache_key]])) {
    .test_idata_cache[[cache_key]] <- make_test_idata() |>
      agg_repertoires(schema)
  }

  .test_idata_cache[[cache_key]]$clone(deep = TRUE)
}


make_grouped_distance_test_idata <- function() {
  receptor_col <- immundata::imd_schema("receptor")

  annotations <- tibble::tibble(
    !!receptor_col := seq_len(8L),
    subject_id = rep(c("P1", "P2"), each = 4L),
    cdr3_aa = c(
      "AAAA", "AAAT", "AATT", "ATTT",
      "CCCC", "CCCT", "CCTT", "CTTT"
    ),
    v_call = "V1",
    j_call = "J1"
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = c("cdr3_aa", "v_call", "j_call")
    ),
    annotations = annotations
  )
}
