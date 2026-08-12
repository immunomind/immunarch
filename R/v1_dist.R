#' @title Immune receptor distances
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to calculate **sequence distances between immune receptors**.
#' These methods help you explore receptor similarity, select thresholds for
#' Immunoglobulin clone assignment, and construct sparse receptor-similarity graphs.
#'
#' ## Available functions
#'
#' The following methods are available.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams dist_hamm
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @section Visualisation:
#' Distance results can be passed directly to [vis()].
#'
#' ## 1) Hamming distances (`dist_hamm`)
#'
#' `vis()` plots the normalized distance to the nearest non-identical neighbor
#' for each receptor by default. This is the distribution used for manual IG
#' clone-threshold selection. When a subject column can be inferred using the
#' same priority as `sample_by = "<auto>"`, one panel is drawn per subject;
#' otherwise, a single pooled distribution is drawn.
#'
#' Use `mode = "all"` to plot every pairwise distance instead of nearest neighbors.
#' Pass one of the following values to `xval` to change what will be plotted:
#' `"norm_dist"` (the default), `"dist"` (non-normalised distnaces), or `"sim"` (similarity).
#' Set `facet = NULL` to have one plot per all input samples, or name a result column to facet by it.
#' `binwidth` and `title` customize the histogram.
#'
#' @name dist
#' @concept Receptor distance
NULL


#' @keywords internal
infer_dist_group_col <- function(columns) {
  candidates <- c(
    "subject_id",
    "donor",
    "donor_id",
    "patient_id",
    "patient",
    "subject",
    immundata::imd_schema("repertoire")
  )
  index <- match(TRUE, candidates %in% columns)

  if (is.na(index)) NULL else candidates[[index]]
}


#' @keywords internal
dist_hamm_impl <- function(
  idata,
  seq_col = "cdr3_aa",
  by = NULL,
  max_dist = NULL,
  min_sim = NULL,
  sample_n = 0L,
  sample_by = "<auto>"
) {
  checkmate::assert_string(seq_col)
  checkmate::assert_number(max_dist, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(
    min_sim,
    lower = 0, upper = 1, finite = TRUE, null.ok = TRUE
  )
  checkmate::assert_count(sample_n)

  if (!is.null(max_dist) && !is.null(min_sim)) {
    cli::cli_abort("Provide only one of {.arg max_dist} and {.arg min_sim}.")
  }

  if (!is.null(max_dist) && max_dist >= 1 && max_dist != floor(max_dist)) {
    cli::cli_abort(
      "{.arg max_dist} must be an integer when it is at least 1."
    )
  }

  if (is.null(by)) {
    by <- character()
  } else {
    checkmate::assert_character(
      by,
      min.len = 1L,
      any.missing = FALSE,
      unique = TRUE
    )
  }

  receptor_col <- immundata::imd_schema("receptor")
  locus_col <- immundata::imd_schema("locus")
  annotation_cols <- colnames(idata$annotations)
  paired <- length(idata$schema_receptor$chains) == 2L

  resolved_sample_by <- NULL
  if (sample_n > 0L) {
    if (!is.null(sample_by)) {
      checkmate::assert_string(sample_by)
    }

    if (identical(sample_by, "<auto>")) {
      resolved_sample_by <- infer_dist_group_col(annotation_cols)

      if (is.null(resolved_sample_by)) {
        cli::cli_abort(c(
          "Could not infer a sampling group from {.field idata$annotations}.",
          "i" = "Set {.arg sample_by} to a column name, or to {.code NULL} for global sampling."
        ))
      }
    } else {
      resolved_sample_by <- sample_by
    }

    if (!is.null(resolved_sample_by) &&
      !resolved_sample_by %in% annotation_cols) {
      cli::cli_abort(
        "Column {.field {resolved_sample_by}} selected by {.arg sample_by} is missing from {.field idata$annotations}."
      )
    }

    if (!is.null(resolved_sample_by) && !resolved_sample_by %in% by) {
      cli::cli_abort(c(
        "Sampling is grouped by {.field {resolved_sample_by}}, but that column is absent from {.arg by}.",
        "i" = "Add {.field {resolved_sample_by}} to {.arg by} to prevent distances between different sampling groups."
      ))
    }
  }

  if (!seq_col %in% idata$schema_receptor$features) {
    cli::cli_abort(
      "{.arg seq_col} ({.field {seq_col}}) must be part of the receptor schema."
    )
  }

  required_cols <- unique(c(
    receptor_col,
    seq_col,
    by,
    if (paired) locus_col
  ))
  missing_cols <- setdiff(required_cols, annotation_cols)

  if (length(missing_cols)) {
    cli::cli_abort(
      "Column(s) [{.field {missing_cols}}] are missing from {.field idata$annotations}."
    )
  }

  duplicated_roles <- intersect(by, c(receptor_col, seq_col))
  if (length(duplicated_roles)) {
    cli::cli_abort(
      "{.arg by} must not contain the receptor ID or sequence column: {.field {duplicated_roles}}."
    )
  }

  by <- unique(c(by, if (paired) locus_col))
  source <- idata$annotations |>
    dplyr::select(dplyr::all_of(c(receptor_col, by, seq_col))) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()

  con <- dbplyr::remote_con(source)

  quote_ident <- function(x) {
    as.character(dbplyr::sql_escape_ident(con, x))
  }
  q_receptor <- quote_ident(receptor_col)
  q_sequence <- quote_ident(seq_col)
  q_by <- vapply(by, quote_ident, character(1))
  source_sql <- as.character(dbplyr::sql_render(source))

  group_select <- if (length(by)) {
    paste0(
      ", ",
      paste0(
        q_by,
        " AS \".idist_group_",
        seq_along(q_by),
        "\"",
        collapse = ", "
      )
    )
  } else {
    ""
  }

  node_group_fields <- if (length(by)) {
    paste0(
      "a.\".idist_group_",
      seq_along(by),
      "\" AS ",
      q_by,
      collapse = ", "
    )
  } else {
    ""
  }

  group_join <- if (length(by)) {
    paste0(
      " AND ",
      paste0(
        "a.\".idist_group_",
        seq_along(by),
        "\" IS NOT DISTINCT FROM b.\".idist_group_",
        seq_along(by),
        "\"",
        collapse = " AND "
      )
    )
  } else {
    ""
  }

  base_nodes_name <- if (sample_n > 0L) "all_nodes" else "nodes"
  nodes_cte <- paste0(
    base_nodes_name, " AS MATERIALIZED (",
    " SELECT ",
    q_receptor, " AS \".idist_receptor\", ",
    q_sequence, " AS \".idist_sequence\", ",
    "length(", q_sequence, ")::BIGINT AS \".idist_length\"",
    group_select,
    " FROM (", source_sql, ") source",
    " WHERE ", q_receptor, " IS NOT NULL",
    " AND ", q_sequence, " IS NOT NULL",
    " AND length(", q_sequence, ") > 0",
    ")"
  )

  if (sample_n > 0L) {
    if (is.null(resolved_sample_by)) {
      sample_group_select <- ""
      sample_partition <- ""
      sample_group_join <- ""
    } else {
      sample_group_index <- match(resolved_sample_by, by)
      sample_group_field <- paste0(
        "\".idist_group_", sample_group_index, "\""
      )
      sample_group_select <- paste0(
        ", ", sample_group_field, " AS \".idist_sample_group\""
      )
      sample_partition <- "PARTITION BY \".idist_sample_group\" "
      sample_group_join <- paste0(
        " AND n.", sample_group_field,
        " IS NOT DISTINCT FROM s.\".idist_sample_group\""
      )
    }

    nodes_cte <- paste0(
      nodes_cte,
      ", sampled_receptors AS MATERIALIZED (",
      " SELECT \".idist_receptor\"",
      if (!is.null(resolved_sample_by)) ", \".idist_sample_group\"" else "",
      " FROM (",
      " SELECT *, row_number() OVER (",
      sample_partition,
      "ORDER BY random()) AS \".idist_sample_rank\"",
      " FROM (",
      " SELECT DISTINCT \".idist_receptor\"",
      sample_group_select,
      " FROM all_nodes",
      ") sampling_groups",
      ") ranked_receptors",
      " WHERE \".idist_sample_rank\" <= ", as.character(sample_n),
      ")",
      ", nodes AS MATERIALIZED (",
      " SELECT n.*",
      " FROM all_nodes n",
      " JOIN sampled_receptors s",
      " ON n.\".idist_receptor\" = s.\".idist_receptor\"",
      sample_group_join,
      ")"
    )
  }

  pair_join <- paste0(
    "a.\".idist_length\" = b.\".idist_length\"",
    group_join,
    " AND a.\".idist_receptor\" < b.\".idist_receptor\""
  )

  if (is.null(min_sim) && (is.null(max_dist) || max_dist == 0)) {
    pair_query <- paste0(
      "SELECT ",
      if (length(by)) paste0(node_group_fields, ", ") else "",
      "a.\".idist_receptor\" AS imd_receptor_id_1",
      ", b.\".idist_receptor\" AS imd_receptor_id_2",
      ", a.\".idist_length\" AS seq_len",
      ", hamming(",
      "a.\".idist_sequence\", b.\".idist_sequence\"",
      ") AS mismatches",
      " FROM nodes a",
      " JOIN nodes b ON ", pair_join
    )

    query <- paste0(
      "WITH ", nodes_cte, ", distances AS (",
      pair_query,
      ")",
      " SELECT ",
      if (length(by)) paste0(paste(q_by, collapse = ", "), ", ") else "",
      "imd_receptor_id_1, imd_receptor_id_2, seq_len",
      ", mismatches AS dist",
      ", 1.0 * mismatches / seq_len AS norm_dist",
      ", 1.0 - 1.0 * mismatches / seq_len AS sim",
      " FROM distances"
    )
  } else {
    norm_bound <- if (!is.null(min_sim)) 1 - min_sim else max_dist
    max_mismatches_sql <- if (!is.null(min_sim) || max_dist < 1) {
      paste0(
        "least(\".idist_length\", floor(",
        as.character(dbplyr::escape(norm_bound, con = con)),
        " * \".idist_length\" + 1e-12)::BIGINT)"
      )
    } else {
      paste0(
        "least(\".idist_length\", ",
        as.integer(max_dist),
        "::BIGINT)"
      )
    }

    block_group_fields <- if (length(by)) {
      paste0(
        "a.\".idist_group_",
        seq_along(by),
        "\"",
        collapse = ", "
      )
    } else {
      ""
    }

    candidate_group_fields <- if (length(by)) {
      paste0(
        "\".idist_group_",
        seq_along(by),
        "\" AS ",
        q_by,
        collapse = ", "
      )
    } else {
      ""
    }

    block_group_join <- if (length(by)) {
      paste0(
        " AND ",
        paste0(
          "a.\".idist_group_",
          seq_along(by),
          "\" IS NOT DISTINCT FROM b.\".idist_group_",
          seq_along(by),
          "\"",
          collapse = " AND "
        )
      )
    } else {
      ""
    }

    query <- paste0(
      "WITH ", nodes_cte,
      ", bounded_nodes AS MATERIALIZED (",
      " SELECT *, ", max_mismatches_sql, " AS max_mismatches",
      " FROM nodes",
      ")",
      ", block_index AS MATERIALIZED (",
      " SELECT n.*, block",
      ", substr(",
      "n.\".idist_sequence\", ",
      "floor(block * n.\".idist_length\" / (n.max_mismatches + 1))::BIGINT + 1, ",
      "(floor((block + 1) * n.\".idist_length\" / (n.max_mismatches + 1)) - ",
      "floor(block * n.\".idist_length\" / (n.max_mismatches + 1)))::BIGINT",
      ") AS token",
      " FROM bounded_nodes n",
      " CROSS JOIN LATERAL range(0, n.max_mismatches + 1) blocks(block)",
      ")",
      ", candidates AS (",
      " SELECT DISTINCT ",
      if (length(by)) paste0(block_group_fields, ", ") else "",
      "a.\".idist_receptor\" AS imd_receptor_id_1",
      ", b.\".idist_receptor\" AS imd_receptor_id_2",
      ", a.\".idist_sequence\" AS sequence_from",
      ", b.\".idist_sequence\" AS sequence_to",
      ", a.\".idist_length\" AS seq_len",
      ", a.max_mismatches",
      " FROM block_index a",
      " JOIN block_index b ON ",
      "a.\".idist_length\" = b.\".idist_length\"",
      block_group_join,
      " AND a.max_mismatches = b.max_mismatches",
      " AND a.block = b.block",
      " AND a.token = b.token",
      " AND a.\".idist_receptor\" < b.\".idist_receptor\"",
      ")",
      ", distances AS (",
      " SELECT *, hamming(sequence_from, sequence_to) AS mismatches",
      " FROM candidates",
      ")",
      " SELECT ",
      if (length(by)) paste0(candidate_group_fields, ", ") else "",
      "imd_receptor_id_1, imd_receptor_id_2, seq_len",
      ", mismatches AS dist",
      ", 1.0 * mismatches / seq_len AS norm_dist",
      ", 1.0 - 1.0 * mismatches / seq_len AS sim",
      " FROM distances",
      " WHERE mismatches <= max_mismatches"
    )
  }

  dplyr::tbl(con, dbplyr::sql(query)) |>
    duckplyr::as_duckdb_tibble()
}


#' @description
#' **1) Hamming distances (`dist_hamm`).** Compute an upper-triangular sparse
#' edge table of Hamming distances between receptor IDs represented in
#' `idata$annotations`.
#' Comparisons are restricted to equal sequence lengths and equal values of
#' `by`. Please mind that the data is not grouped by repertoires automatically.
#' For paired-chain receptor schemas, the canonical locus column is added to `by` automatically,
#' so distances are computed separately per locus. Add it explicitly if you have a
#' non-AIRR-C standard locus column name.
#'
#' @param seq_col Name of a sequence feature in the receptor schema.
#' @param by Optional columns from `idata$annotations` that must be equal within
#'   a receptor pair. For B-cell clone assignment, a typical choice is
#'   `c("v_call", "j_call")`; add `imd_repertoire_id` when comparisons must stay
#'   within repertoires, or `subject_id` if your data has such a dedicated column for
#'   subject-specific metadata. Locus grouping is also always applied for paired-chain receptor schemas.
#' @param max_dist Maximum distance. Values in `(0, 1)` bound normalized
#'   Hamming distance, integer values which are `>=1` bound raw
#'   Hamming distance. `NULL` (the default) and the legacy value `0` disable
#'   bounding. Mutually exclusive with `min_sim`.
#' @param min_sim Minimum normalized Hamming similarity in `[0, 1]`. Mutually
#'   exclusive with `max_dist`.
#' @param sample_n Maximum number of distinct receptor IDs to sample before
#'   distance calculation within each `sample_by` group. `0` (the default)
#'   disables sampling.
#' @param sample_by Select a column to sample by. The default,
#'   `"<auto>"`, uses the first available column among `subject_id`, `donor`,
#'   `donor_id`, `patient_id`, `patient`, `subject`, and `imd_repertoire_id`.
#'   The resolved column must also be supplied through `by`. Set to `NULL` to
#'   sample from the whole dataset without accounting for groups.
#'   For obvious reasons, this argument is completely ignored when `sample_n = 0`.
#'
#' @return
#'
#' ## 1) Hamming distances (`dist_hamm`)
#' A lazy duckplyr table with:
#' * grouping columns supplied through `by`
#' * `locus` for paired-chain receptor schemas
#' * `imd_receptor_id_1`, `imd_receptor_id_2` -- canonical receptor identifiers,
#'   with `imd_receptor_id_1 < imd_receptor_id_2`
#' * `seq_len`
#' * `dist` -- raw Hamming distance
#' * `norm_dist` -- `dist / seq_len`
#' * `sim` -- normalized Hamming similarity, `1 - norm_dist`
#'
#' Call [dplyr::collect()] only when the edge table is small enough for memory.
#'
#' @examples
#' idata <- get_test_idata()
#'
#' # Build the computing query
#' res <- dist_hamm(idata)
#' # Check how results would look like
#' res
#' # Run the computation query and upload the results to the memory
#' res <- res |> collect()
#'
#' # Bound Hamming distance.
#' # Non-normalised -- no receptors with more than 1 mismatch in the output
#' res <- dist_hamm(idata, max_dist = 1) |> collect()
#' res
#'
#' # Normalised -- no receptors with more than .2 normalised distance in the output
#' res <- dist_hamm(idata, max_dist = .2) |> collect()
#' res
#'
#' # Similarity, inverse to normalised distance -- no receptors with
#' # less than 0.9 similarity in the output
#' # The results are the same as for the previous normalised run
#' res <- dist_hamm(idata, min_sim = .8) |> collect()
#' res
#'
#' # Plot normalized nearest-non-identical-neighbor distances.
#' vis(dist_hamm(idata))
#' # Plot all raw pairwise distances instead.
#' vis(dist_hamm(idata), mode = "all", xval = "dist")
#'
#' \dontrun{
#' #
#' # Typical use case for BCR data
#' #
#'
#' # BCR-style V/J/length blocking with at most three mismatches, inside patients
#' clone_edges <- dist_hamm(
#'   idata_bcr,
#'   seq_col = "junction",
#'   by = c("patient", "v_call", "j_call"),
#'   max_dist = 3
#' )
#' }
#'
#' @rdname dist
#' @concept Receptor distance
#' @export
dist_hamm <- register_immunarch_method(
  core = dist_hamm_impl,
  family = "dist",
  name = "hamm",
  required_cols = c("seq_col", "by"),
  need_repertoires = FALSE
)
