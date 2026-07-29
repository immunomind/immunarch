#' @title Receptor distances
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions that constructs sparse receptor-pair distance tables
#' from the receptors in one `ImmunData` object.
#'
#' ## Available functions
#'
#' * `dist_hamm()` computes Hamming distance between equal-length receptor
#'   sequences.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams dist_hamm
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @name dist
#' @concept Receptor distance
NULL


#' @keywords internal
dist_hamm_impl <- function(
    idata,
    seq_col = "cdr3_aa",
    by = NULL,
    max_dist = 0) {
  checkmate::assert_string(seq_col)
  checkmate::assert_number(max_dist, lower = 0, finite = TRUE)

  if (max_dist >= 1 && max_dist != floor(max_dist)) {
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

  nodes_cte <- paste0(
    "nodes AS MATERIALIZED (",
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

  pair_join <- paste0(
    "a.\".idist_length\" = b.\".idist_length\"",
    group_join,
    " AND a.\".idist_receptor\" < b.\".idist_receptor\""
  )

  if (max_dist == 0) {
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
      " FROM distances"
    )
  } else {
    max_mismatches_sql <- if (max_dist < 1) {
      paste0(
        "least(\".idist_length\", floor(",
        as.character(dbplyr::escape(max_dist, con = con)),
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
      " FROM distances",
      " WHERE mismatches <= max_mismatches"
    )
  }

  dplyr::tbl(con, dbplyr::sql(query)) |>
    duckplyr::as_duckdb_tibble()
}


#' @description `dist_hamm()` computes an upper-triangular sparse edge table of
#' Hamming distances between receptor IDs represented in `idata$annotations`.
#' Comparisons are restricted to equal sequence lengths and equal values of
#' `by`; repertoires are not an implicit grouping boundary. For paired-chain
#' receptor schemas, the canonical locus column is added to `by` automatically,
#' so distances are computed separately per locus.
#'
#' With `max_dist = 0`, all eligible receptor pairs are returned. With a
#' positive bound, exact `(k + 1)`-block candidate generation is used before
#' Hamming verification; this avoids constructing the complete Cartesian
#' product when the bound is selective.
#'
#' @param seq_col Name of a sequence feature in the receptor schema.
#' @param by Optional columns from `idata$annotations` that must be equal within
#'   a receptor pair. For B-cell clone assignment, a typical choice is
#'   `c("v_call", "j_call")`; add `imd_repertoire_id` when comparisons must stay
#'   within repertoires. Sequence length grouping is always applied. Locus
#'   grouping is also always applied for paired-chain receptor schemas.
#' @param max_dist Distance bound. The default `0` disables bounding. Values in
#'   `(0, 1)` bound normalized Hamming distance; integer values greater than or
#'   equal to `1` bound raw Hamming distance.
#'
#' @return A lazy duckplyr table with:
#' * grouping columns supplied through `by`
#' * `locus` for paired-chain receptor schemas
#' * `imd_receptor_id_1`, `imd_receptor_id_2` -- canonical receptor identifiers,
#'   with `imd_receptor_id_1 < imd_receptor_id_2`
#' * `seq_len`
#' * `dist` -- raw Hamming distance
#' * `norm_dist` -- `dist / seq_len`
#'
#' Call [dplyr::collect()] only when the edge table is small enough for memory.
#'
#' @examples
#' \dontrun{
#' # Direct Hamming distance between receptors.
#' edges <- dist_hamm(immdata, seq_col = "junction")
#'
#' # BCR-style V/J/length blocking with at most four mismatches.
#' clone_edges <- dist_hamm(
#'   immdata,
#'   seq_col = "junction",
#'   by = c("v_call", "j_call"),
#'   max_dist = 4
#' )
#'
#' # Bound normalized Hamming distance.
#' clone_edges_norm <- dist_hamm(
#'   immdata,
#'   seq_col = "junction",
#'   by = c("v_call", "j_call"),
#'   max_dist = 0.1
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
