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
    .test_idata_cache$base <- get_test_immundata()
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
