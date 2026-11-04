#' @title Immune receptor distances
#'
#' @description
#'
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
#' @param idata An [immundata::ImmunData] object.
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
dist_hamm_impl <- function(
  idata,
  seq_col = "cdr3_aa",
  by = NULL,
  max_dist = NULL,
  min_sim = NULL,
  sample_n = 0L,
  sample_by = "<auto>",
  mode = c("pairs", "graph")
) {
  mode <- tryCatch(
    match.arg(mode),
    error = function(cnd) {
      cli::cli_abort(
        "{.arg mode} must be one of {.val pairs} or {.val graph}.",
        parent = cnd
      )
    }
  )
  checkmate::assert_string(seq_col)
  table_bound <- is.data.frame(max_dist)
  if (!table_bound) {
    checkmate::assert_number(
      max_dist,
      lower = 0, finite = TRUE, null.ok = TRUE
    )
  }
  checkmate::assert_number(
    min_sim,
    lower = 0, upper = 1, finite = TRUE, null.ok = TRUE
  )
  checkmate::assert_count(sample_n)

  if (!is.null(max_dist) && !is.null(min_sim)) {
    cli::cli_abort("Provide only one of {.arg max_dist} and {.arg min_sim}.")
  }

  if (!table_bound &&
    !is.null(max_dist) && max_dist >= 1 && max_dist != floor(max_dist)) {
    cli::cli_abort(
      "{.arg max_dist} must be an integer when it is at least 1."
    )
  }

  checkmate::assert_character(
    by,
    min.len = 1L,
    any.missing = FALSE,
    unique = TRUE,
    null.ok = TRUE
  )

  receptor_col <- immundata::imd_schema("receptor")
  locus_col <- immundata::imd_schema("locus")
  annotation_cols <- colnames(idata$annotations)
  paired <- length(idata$schema_receptor$chains) == 2L

  resolved_sample_by <- NULL
  if (sample_n > 0L) {
    checkmate::assert_string(sample_by, null.ok = TRUE)

    if (identical(sample_by, "<auto>")) {
      candidates <- c(
        "subject_id",
        "donor",
        "donor_id",
        "patient_id",
        "patient",
        "subject",
        immundata::imd_schema("repertoire")
      )
      index <- match(TRUE, candidates %in% annotation_cols)
      resolved_sample_by <- if (is.na(index)) NULL else candidates[[index]]

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
  reserved_by <- if (mode == "graph") {
    c(
      ".idist_record_type", "from", "to", "node_id",
      "seq_len", "dist", "norm_dist", "sim"
    )
  } else {
    c(
      "imd_receptor_id_1", "imd_receptor_id_2",
      "seq_len", "dist", "norm_dist", "sim"
    )
  }
  conflicting_by <- intersect(by, reserved_by)
  if (length(conflicting_by)) {
    cli::cli_abort(
      "Grouping column(s) conflict with reserved {.arg mode} output columns: {.field {conflicting_by}}."
    )
  }

  threshold_spec <- NULL
  if (table_bound) {
    thresholds <- max_dist
    checkmate::assert_data_frame(
      thresholds,
      min.rows = 1L,
      min.cols = 2L
    )

    if (anyDuplicated(names(thresholds)) ||
      !"threshold" %in% names(thresholds)) {
      cli::cli_abort(
        "A {.arg max_dist} table must have unique columns, including {.field threshold}."
      )
    }

    keys <- setdiff(names(thresholds), "threshold")
    missing_by <- setdiff(keys, by)
    if (length(missing_by)) {
      cli::cli_abort(
        "Threshold key column(s) must be included in {.arg by}: {.field {missing_by}}."
      )
    }

    checkmate::assert_numeric(
      thresholds$threshold,
      finite = TRUE,
      any.missing = FALSE
    )
    if (any(thresholds$threshold <= 0 | thresholds$threshold >= 1)) {
      cli::cli_abort("Table thresholds must be strictly between 0 and 1.")
    }

    if (anyNA(thresholds[keys])) {
      cli::cli_abort("Threshold key columns must not contain missing values.")
    }

    if (anyDuplicated(thresholds[keys])) {
      cli::cli_abort("Threshold key combinations must be unique.")
    }

    input_keys <- idata$annotations |>
      dplyr::select(dplyr::all_of(keys)) |>
      dplyr::distinct() |>
      dplyr::collect()
    threshold_keys <- dplyr::select(thresholds, dplyr::all_of(keys))

    if (nrow(dplyr::anti_join(threshold_keys, input_keys, by = keys))) {
      cli::cli_abort(
        "The threshold table contains a key combination absent from the input."
      )
    }

    if (nrow(dplyr::anti_join(input_keys, threshold_keys, by = keys))) {
      cli::cli_abort(
        "The threshold table must cover every input key combination."
      )
    }

    threshold_spec <- list(
      keys = keys,
      table = dplyr::rename(
        thresholds,
        .idist_threshold = "threshold"
      )
    )
  }

  source <- idata$annotations |>
    dplyr::select(dplyr::all_of(c(receptor_col, by, seq_col))) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()
  con <- dbplyr::remote_con(source)

  if (!is.null(threshold_spec)) {
    source <- source |>
      dplyr::inner_join(
        dbplyr::copy_inline(con, threshold_spec$table),
        by = threshold_spec$keys
      )
  }

  empty_sql <- DBI::SQL("")
  group_aliases <- paste0(".idist_group_", seq_along(by))
  qualified_columns <- function(columns, table = NULL, output_names = NULL) {
    table_prefix <- if (is.null(table)) {
      empty_sql
    } else {
      DBI::SQL(paste0(table, "."))
    }
    fields <- if (is.null(output_names)) {
      glue::glue_sql("{table_prefix}{`columns`}", .con = con)
    } else {
      glue::glue_sql(
        "{table_prefix}{`columns`} AS {`output_names`}",
        .con = con
      )
    }
    glue::glue_sql_collapse(fields, sep = ", ")
  }

  if (length(by)) {
    fields <- qualified_columns(by, output_names = group_aliases)
    group_select <- glue::glue_sql(", {fields}", .con = con)

    group_predicates <- glue::glue_sql(
      "a.{`group_aliases`} IS NOT DISTINCT FROM b.{`group_aliases`}",
      .con = con
    )
    predicates <- glue::glue_sql_collapse(group_predicates, sep = " AND ")
    group_join <- glue::glue_sql(" AND {predicates}", .con = con)
  } else {
    group_select <- empty_sql
    group_join <- empty_sql
  }

  threshold_select <- if (!is.null(threshold_spec)) {
    DBI::SQL(", \".idist_threshold\"")
  } else {
    empty_sql
  }

  base_receptors_name <- if (sample_n > 0L) "all_receptors" else "receptors"
  receptors_cte <- glue::glue_sql(
    "
    {`base_receptors_name`} AS MATERIALIZED (
      SELECT
        {`receptor_col`} AS \".idist_receptor\",
        {`seq_col`} AS \".idist_sequence\",
        length({`seq_col`})::BIGINT AS \".idist_length\"
        {group_select}{threshold_select}
      FROM ({source}) source
      WHERE {`receptor_col`} IS NOT NULL
        AND {`seq_col`} IS NOT NULL
        AND length({`seq_col`}) > 0
    )",
    source = DBI::SQL(dbplyr::sql_render(source)),
    .con = con
  )

  if (sample_n > 0L) {
    if (is.null(resolved_sample_by)) {
      sample_group_select <- empty_sql
      sample_partition <- empty_sql
      sample_group_join <- empty_sql
      sampled_group_field <- empty_sql
    } else {
      sample_group_index <- match(resolved_sample_by, by)
      sample_group_field <- group_aliases[[sample_group_index]]
      sample_group_select <- glue::glue_sql(
        ", {`sample_group_field`} AS \".idist_sample_group\"",
        .con = con
      )
      sample_partition <- DBI::SQL(
        "PARTITION BY \".idist_sample_group\" "
      )
      sample_group_join <- glue::glue_sql(
        " AND n.{`sample_group_field`} IS NOT DISTINCT FROM s.\".idist_sample_group\"",
        .con = con
      )
      sampled_group_field <- DBI::SQL(", \".idist_sample_group\"")
    }

    receptors_cte <- glue::glue_sql(
      "
      {receptors_cte}, sampling_draws AS MATERIALIZED (
        SELECT *, random() AS \".idist_sample_draw\"
        FROM (
          SELECT DISTINCT \".idist_receptor\"{sample_group_select}
          FROM all_receptors
        ) sampling_groups
      ), sampled_receptors AS MATERIALIZED (
        SELECT \".idist_receptor\"{sampled_group_field}
        FROM (
          SELECT *, row_number() OVER (
            {sample_partition}ORDER BY \".idist_sample_draw\"
          ) AS \".idist_sample_rank\"
          FROM sampling_draws
        ) ranked_receptors
        WHERE \".idist_sample_rank\" <= {sample_n}
      ), receptors AS MATERIALIZED (
        SELECT n.*
        FROM all_receptors n
        JOIN sampled_receptors s
          ON n.\".idist_receptor\" = s.\".idist_receptor\"{sample_group_join}
      )",
      .con = con
    )
  }

  if (length(by)) {
    signature_groups <- qualified_columns(group_aliases)
    signature_group_columns <- glue::glue_sql(
      ", {signature_groups}",
      .con = con
    )
    signature_order <- glue::glue_sql(
      "{signature_groups}, \".idist_sequence\"",
      .con = con
    )
    graph_groups <- qualified_columns(by)
    graph_output_columns <- glue::glue_sql(
      "{graph_groups}, ",
      .con = con
    )
    qualified_membership_groups <- qualified_columns(group_aliases, "r")
    membership_group_columns <- glue::glue_sql(
      ", {qualified_membership_groups}",
      .con = con
    )
    membership_output_groups <- qualified_columns(group_aliases, "m", by)
    membership_output_prefix <- glue::glue_sql(
      "{membership_output_groups}, ",
      .con = con
    )
    pair_output_groups <- qualified_columns(group_aliases, "p", by)
    pair_output_prefix <- glue::glue_sql(
      "{pair_output_groups}, ",
      .con = con
    )
    qualified_distance_groups <- qualified_columns(group_aliases, "a")
    distance_group_columns <- glue::glue_sql(
      "{qualified_distance_groups}, ",
      .con = con
    )
    qualified_pair_groups <- qualified_columns(group_aliases, "d")
    distance_pair_groups <- glue::glue_sql(
      "{qualified_pair_groups}, ",
      .con = con
    )
    qualified_node_groups <- qualified_columns(group_aliases, "n")
    node_pair_groups <- glue::glue_sql(
      "{qualified_node_groups}, ",
      .con = con
    )
    edge_output_groups <- qualified_columns(group_aliases, "d", by)
    edge_output_prefix <- glue::glue_sql(
      "{edge_output_groups}, ",
      .con = con
    )
    signature_predicates <- glue::glue_sql_collapse(
      glue::glue_sql(
        "r.{`group_aliases`} IS NOT DISTINCT FROM n.{`group_aliases`}",
        .con = con
      ),
      sep = " AND "
    )
    signature_group_join <- glue::glue_sql(
      " AND {signature_predicates}",
      .con = con
    )
  } else {
    signature_group_columns <- empty_sql
    signature_order <- DBI::SQL("\".idist_sequence\"")
    graph_output_columns <- empty_sql
    membership_group_columns <- empty_sql
    membership_output_prefix <- empty_sql
    pair_output_prefix <- empty_sql
    distance_group_columns <- empty_sql
    distance_pair_groups <- empty_sql
    node_pair_groups <- empty_sql
    edge_output_prefix <- empty_sql
    signature_group_join <- empty_sql
  }

  contraction_ctes <- glue::glue_sql(
    "
    {receptors_cte}, sequence_nodes AS MATERIALIZED (
      SELECT row_number() OVER (
        ORDER BY {signature_order}
      )::BIGINT AS node_id, signatures.*
      FROM (
        SELECT DISTINCT
          \".idist_sequence\", \".idist_length\"
          {signature_group_columns}{threshold_select}
        FROM receptors
      ) signatures
    ), memberships AS MATERIALIZED (
      SELECT DISTINCT
        r.\".idist_receptor\"{membership_group_columns}, n.node_id
      FROM receptors r
      JOIN sequence_nodes n
        ON r.\".idist_sequence\" IS NOT DISTINCT FROM n.\".idist_sequence\"
        {signature_group_join}
    )",
    .con = con
  )

  node_pair_join <- glue::glue_sql(
    "
    a.\".idist_length\" = b.\".idist_length\"{group_join}
      AND a.node_id < b.node_id",
    .con = con
  )

  unbounded <- is.null(threshold_spec) &&
    is.null(min_sim) && (is.null(max_dist) || max_dist == 0)
  if (unbounded) {
    distance_ctes <- glue::glue_sql(
      "
      graph_distances AS (
        SELECT
          {distance_group_columns}a.node_id AS \"from\",
          b.node_id AS \"to\",
          a.\".idist_length\" AS seq_len,
          hamming(
            a.\".idist_sequence\", b.\".idist_sequence\"
          ) AS mismatches
        FROM sequence_nodes a
        JOIN sequence_nodes b ON {node_pair_join}
      )",
      .con = con
    )
  } else {
    norm_bound_sql <- if (!is.null(threshold_spec)) {
      DBI::SQL("\".idist_threshold\"")
    } else {
      norm_bound <- if (!is.null(min_sim)) 1 - min_sim else max_dist
      norm_bound
    }
    normalized_bound <- !is.null(threshold_spec) ||
      !is.null(min_sim) || max_dist < 1
    max_mismatches_sql <- if (normalized_bound) {
      glue::glue_sql(
        "least(\".idist_length\", floor(
          {norm_bound_sql} * \".idist_length\" + 1e-12
        )::BIGINT)",
        .con = con
      )
    } else {
      max_dist <- as.integer(max_dist)
      glue::glue_sql(
        "least(\".idist_length\", {max_dist}::BIGINT)",
        .con = con
      )
    }

    distance_ctes <- glue::glue_sql(
      "
      bounded_nodes AS MATERIALIZED (
        SELECT *, {max_mismatches_sql} AS max_mismatches
        FROM sequence_nodes
      ), block_index AS MATERIALIZED (
        SELECT
          n.*, block,
          substr(
            n.\".idist_sequence\",
            floor(
              block * n.\".idist_length\" / (n.max_mismatches + 1)
            )::BIGINT + 1,
            (floor(
              (block + 1) * n.\".idist_length\" / (n.max_mismatches + 1)
            ) - floor(
              block * n.\".idist_length\" / (n.max_mismatches + 1)
            ))::BIGINT
          ) AS token
        FROM bounded_nodes n
        CROSS JOIN LATERAL range(
          0, n.max_mismatches + 1
        ) blocks(block)
      ), candidates AS (
        SELECT DISTINCT
          {distance_group_columns}a.node_id AS \"from\",
          b.node_id AS \"to\",
          a.\".idist_sequence\" AS sequence_from,
          b.\".idist_sequence\" AS sequence_to,
          a.\".idist_length\" AS seq_len,
          a.max_mismatches
        FROM block_index a
        JOIN block_index b
          ON a.\".idist_length\" = b.\".idist_length\"{group_join}
          AND a.max_mismatches = b.max_mismatches
          AND a.block = b.block
          AND a.token = b.token
          AND a.node_id < b.node_id
      ), measured_distances AS (
        SELECT *, hamming(sequence_from, sequence_to) AS mismatches
        FROM candidates
      ), graph_distances AS (
        SELECT *
        FROM measured_distances
        WHERE mismatches <= max_mismatches
      )",
      .con = con
    )
  }

  materialize <- function(query, pattern) {
    dplyr::tbl(con, dbplyr::sql(query)) |>
      duckplyr::as_duckdb_tibble() |>
      dplyr::compute(
        name = basename(tempfile(pattern = pattern)),
        prudence = "stingy"
      )
  }

  if (mode == "graph") {
    graph_query <- glue::glue_sql(
      "
      WITH {contraction_ctes}, {distance_ctes}, graph_rows AS (
        SELECT
          {edge_output_prefix}'edge' AS \".idist_record_type\",
          d.\"from\", d.\"to\",
          NULL AS {`receptor_col`}, NULL::BIGINT AS node_id,
          d.seq_len, d.mismatches AS dist,
          1.0 * d.mismatches / d.seq_len AS norm_dist,
          1.0 - 1.0 * d.mismatches / d.seq_len AS sim
        FROM graph_distances d

        UNION ALL

        SELECT
          {membership_output_prefix}'membership' AS \".idist_record_type\",
          NULL::BIGINT AS \"from\", NULL::BIGINT AS \"to\",
          m.\".idist_receptor\" AS {`receptor_col`}, m.node_id,
          NULL::BIGINT AS seq_len, NULL::BIGINT AS dist,
          NULL::DOUBLE AS norm_dist, NULL::DOUBLE AS sim
        FROM memberships m
      )
      SELECT
        {graph_output_columns}\".idist_record_type\", \"from\", \"to\",
        {`receptor_col`}, node_id, seq_len, dist, norm_dist, sim
      FROM graph_rows",
      .con = con
    )

    graph_rows <- materialize(graph_query, "immunarch_dist_hamm_graph_")
    out <- list(
      edges = graph_rows |>
        dplyr::filter(.idist_record_type == "edge") |>
        dplyr::select(dplyr::all_of(c(
          by, "from", "to", "seq_len", "dist", "norm_dist", "sim"
        ))),
      membership = graph_rows |>
        dplyr::filter(.idist_record_type == "membership") |>
        dplyr::select(dplyr::all_of(c(by, receptor_col, "node_id")))
    )
    class(out) <- c("immunarch_sequence_graph", class(out))
    return(out)
  }

  pairs_query <- glue::glue_sql(
    "
    WITH {contraction_ctes}, {distance_ctes}, expanded_pairs AS (
      SELECT
        {distance_pair_groups}
        least(m1.\".idist_receptor\", m2.\".idist_receptor\") AS receptor_1,
        greatest(m1.\".idist_receptor\", m2.\".idist_receptor\") AS receptor_2,
        d.seq_len, d.mismatches
      FROM graph_distances d
      JOIN memberships m1 ON d.\"from\" = m1.node_id
      JOIN memberships m2 ON d.\"to\" = m2.node_id
      WHERE m1.\".idist_receptor\" <> m2.\".idist_receptor\"

      UNION ALL

      SELECT
        {node_pair_groups}
        m1.\".idist_receptor\" AS receptor_1,
        m2.\".idist_receptor\" AS receptor_2,
        n.\".idist_length\" AS seq_len, 0::BIGINT AS mismatches
      FROM sequence_nodes n
      JOIN memberships m1 ON n.node_id = m1.node_id
      JOIN memberships m2 ON n.node_id = m2.node_id
        AND m1.\".idist_receptor\" < m2.\".idist_receptor\"
    )
    SELECT
      {pair_output_prefix}p.receptor_1 AS imd_receptor_id_1,
      p.receptor_2 AS imd_receptor_id_2, p.seq_len,
      p.mismatches AS dist,
      1.0 * p.mismatches / p.seq_len AS norm_dist,
      1.0 - 1.0 * p.mismatches / p.seq_len AS sim
    FROM expanded_pairs p",
    .con = con
  )

  materialize(pairs_query, "immunarch_dist_hamm_")
}


#' @description
#' **1) Hamming distances (`dist_hamm`).** Compute an upper-triangular sparse
#' edge table of Hamming distances between receptor IDs represented in
#' the receptor annotations in `idata`.
#' Comparisons are restricted to equal sequence lengths and equal values of
#' `by`. The data are not grouped by repertoire automatically.
#' For paired-chain receptor schemas, the canonical locus column is added to `by` automatically,
#' so distances are computed separately per locus. Add it explicitly if you have a
#' non-AIRR-C standard locus column name.
#'
#' @param seq_col A single character string naming a sequence feature in the
#'   receptor schema.
#' @param by A character vector or `NULL`. These optional annotation columns
#'   must be equal within a receptor pair. For B-cell clone assignment, a
#'   typical choice is
#'   `c("v_call", "j_call")`; add `imd_repertoire_id` when comparisons must stay
#'   within repertoires, or `subject_id` if your data has such a dedicated column for
#'   subject-specific metadata. Locus grouping is also always applied for
#'   paired-chain receptor schemas.
#' @param max_dist A single nonnegative number, a keyed threshold table, or
#'   `NULL`. Scalar values in `(0, 1)` bound normalized Hamming distance;
#'   integer values which are `>=1` bound raw Hamming distance. `NULL` (the
#'   default) and the legacy scalar value `0` disable bounding. A threshold
#'   table must contain one or more key columns plus a numeric `threshold`
#'   column strictly inside `(0, 1)`. Its key columns must be included in `by`
#'   and must exactly cover the corresponding input key combinations. Mutually
#'   exclusive with `min_sim`.
#' @param min_sim A single number in `[0, 1]` or `NULL`, giving the minimum
#'   normalized Hamming similarity. It is mutually exclusive with `max_dist`.
#' @param sample_n A nonnegative whole number giving the maximum number of
#'   distinct receptor IDs to sample before
#'   distance calculation within each `sample_by` group. `0` (the default)
#'   disables sampling.
#' @param sample_by A single character string or `NULL` selecting the column to
#'   sample by. The default,
#'   `"<auto>"`, uses the first available column among `subject_id`, `donor`,
#'   `donor_id`, `patient_id`, `patient`, `subject`, and `imd_repertoire_id`.
#'   The resolved column must also be supplied through `by`. Set to `NULL` to
#'   sample from the whole dataset without accounting for groups.
#'   This argument is ignored when `sample_n = 0`.
#' @param mode Output representation. `"pairs"` (the default) returns the
#'   receptor-pair table described below. `"graph"` returns contracted
#'   sequence-node edges and receptor-to-node membership. Grouping columns
#'   cannot use names reserved for structural columns in the selected output.
#'
#' @return
#'
#' ## 1) Hamming distances (`dist_hamm`)
#' With `mode = "pairs"`, a lazy duckplyr table with:
#' * grouping columns supplied through `by`
#' * `locus` for paired-chain receptor schemas
#' * `imd_receptor_id_1`, `imd_receptor_id_2` -- canonical receptor identifiers,
#'   with `imd_receptor_id_1 < imd_receptor_id_2`
#' * `seq_len`
#' * `dist` -- raw Hamming distance
#' * `norm_dist` -- `dist / seq_len`
#' * `sim` -- normalized Hamming similarity, `1 - norm_dist`
#'
#' With `mode = "graph"`, an `immunarch_sequence_graph` object containing two
#' lazy duckplyr tables. `edges` contains the effective grouping columns,
#' globally unique `from` and `to` node IDs, and `seq_len`, `dist`, `norm_dist`,
#' and `sim`. `membership` contains the effective grouping columns, the
#' canonical receptor ID, and `node_id`. A receptor can have multiple membership
#' rows, and membership includes nodes without edges. Nodes represent distinct
#' combinations of the full effective comparison group and a non-missing,
#' non-empty sequence; sequences themselves are not returned.
#'
#' Call [dplyr::collect()] only when a returned table is small enough for memory.
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
dist_hamm <- register_method_dist(
  core = dist_hamm_impl,
  family = "dist",
  name = "hamm",
  required_cols = c("seq_col", "by"),
  need_repertoires = FALSE
)
