#' @title Connected-component clustering for receptor graphs
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions for clustering sparse receptor graphs represented by
#' an `ImmunData` object and an `immunarch_sequence_graph`.
#'
#' @param idata An [immundata::ImmunData] object. Its receptor annotations
#'   provide the graph nodes through the canonical receptor-ID column and, when
#'   requested, repertoire metadata.
#'
#' @name clust
#' @concept Clustering
NULL


clust_sql_columns <- function(con, columns, alias = NULL, trailing = FALSE) {
  if (!length(columns)) {
    return(DBI::SQL(""))
  }
  qualifier <- if (is.null(alias)) "" else paste0(alias, ".")
  fields <- glue::glue_sql_collapse(
    glue::glue_sql("{DBI::SQL(qualifier)}{`columns`}", .con = con),
    sep = ", "
  )
  if (trailing) glue::glue_sql("{fields}, ", .con = con) else fields
}


clust_sql_join <- function(con, columns, left, right) {
  if (!length(columns)) {
    return(DBI::SQL(""))
  }
  predicates <- glue::glue_sql_collapse(
    glue::glue_sql(
      "{DBI::SQL(paste0(left, '.'))}{`columns`} IS NOT DISTINCT FROM {DBI::SQL(paste0(right, '.'))}{`columns`}",
      .con = con
    ),
    sep = " AND "
  )
  glue::glue_sql(" AND {predicates}", .con = con)
}


clust_cc_impl <- function(
  idata,
  edges,
  by = NULL,
  cluster_col = "cluster_id",
  size_col = "cluster_size"
) {
  if (!inherits(edges, "immunarch_sequence_graph")) {
    cli::cli_abort(c(
      "{.fn clust_cc} requires an {.cls immunarch_sequence_graph}.",
      "i" = "Create it with {.code dist_hamm(..., mode = \"graph\")}.",
      "x" = "Ordinary receptor edge tables are not accepted."
    ))
  }
  checkmate::assert_character(
    by,
    any.missing = FALSE,
    unique = TRUE,
    null.ok = TRUE
  )
  checkmate::assert_string(cluster_col)
  checkmate::assert_string(size_col)
  by <- if (is.null(by)) character() else by

  receptor_col <- immundata::imd_schema("receptor")
  if (receptor_col %in% by) {
    cli::cli_abort(
      "The canonical receptor-ID column must not also appear in {.arg by}."
    )
  }

  missing_scope_cols <- setdiff(c(by, receptor_col), colnames(idata$annotations))
  if (length(missing_scope_cols)) {
    cli::cli_abort(
      "Column(s) missing from receptor annotations: {.field {missing_scope_cols}}."
    )
  }
  valid_contract <- is.list(edges) &&
    all(c("edges", "membership") %in% names(edges)) &&
    all(c("from", "to") %in% colnames(edges$edges)) &&
    all(c(receptor_col, "node_id") %in% colnames(edges$membership))
  if (!valid_contract) {
    cli::cli_abort(
      "The {.cls immunarch_sequence_graph} has an invalid {.fn dist_hamm} graph contract."
    )
  }

  universe <- idata$annotations |>
    dplyr::select(dplyr::all_of(c(by, receptor_col))) |>
    dplyr::filter(!is.na(!!rlang::sym(receptor_col))) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()
  membership_scope <- intersect(by, colnames(edges$membership))
  membership <- edges$membership |>
    dplyr::select(dplyr::all_of(c(
      membership_scope, receptor_col, "node_id"
    ))) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()
  distance_edges <- edges$edges |>
    dplyr::select("from", "to") |>
    dplyr::distinct() |>
    duckplyr::as_tbl()

  con <- dbplyr::remote_con(universe)
  graph_cons <- list(
    dbplyr::remote_con(membership),
    dbplyr::remote_con(distance_edges)
  )
  if (!all(vapply(graph_cons, identical, logical(1), y = con))) {
    cli::cli_abort(
      "{.field idata} and the sequence graph must use the same DuckDB connection."
    )
  }

  prefix <- function(alias = NULL) {
    clust_sql_columns(con, by, alias, trailing = TRUE)
  }
  null_join <- function(columns, left, right) {
    clust_sql_join(con, columns, left, right)
  }

  scoped_membership_sql <- glue::glue_sql(
    "
    SELECT DISTINCT {prefix('u')}u.{`receptor_col`}, m.node_id
    FROM ({universe}) u
    JOIN ({membership}) m
      ON u.{`receptor_col`} = m.{`receptor_col`}{null_join(membership_scope, 'u', 'm')}
    WHERE m.node_id IS NOT NULL",
    universe = DBI::SQL(dbplyr::sql_render(universe)),
    membership = DBI::SQL(dbplyr::sql_render(membership)),
    .con = con
  )
  scoped_membership <- dplyr::tbl(con, dbplyr::sql(scoped_membership_sql)) |>
    dplyr::compute(
      name = basename(tempfile(pattern = "immunarch_clust_membership_"))
    )
  scoped_nodes <- scoped_membership |>
    dplyr::select(dplyr::all_of(c(by, "node_id"))) |>
    dplyr::distinct()

  scoped_edges_sql <- glue::glue_sql(
    "
    WITH scoped_nodes AS MATERIALIZED ({scoped_nodes}),
    distance_edges AS (
      SELECT DISTINCT {prefix('a')}e.\"from\", e.\"to\"
      FROM ({distance_edges}) e
      JOIN scoped_nodes a ON e.\"from\" = a.node_id
      JOIN scoped_nodes b
        ON e.\"to\" = b.node_id{null_join(by, 'a', 'b')}
      WHERE e.\"from\" IS NOT NULL AND e.\"to\" IS NOT NULL
    ), membership_anchors AS (
      SELECT m.*, min(m.node_id) OVER (
        PARTITION BY {prefix('m')}m.{`receptor_col`}
      ) AS anchor
      FROM ({scoped_membership}) m
    ), identity_edges AS (
      SELECT DISTINCT {prefix('a')}a.anchor AS \"from\", a.node_id AS \"to\"
      FROM membership_anchors a
      WHERE a.anchor <> a.node_id
    )
    SELECT * FROM distance_edges
    UNION
    SELECT * FROM identity_edges",
    scoped_nodes = DBI::SQL(dbplyr::sql_render(scoped_nodes)),
    distance_edges = DBI::SQL(dbplyr::sql_render(distance_edges)),
    scoped_membership = DBI::SQL(dbplyr::sql_render(scoped_membership)),
    .con = con
  )
  scoped_edges <- dplyr::tbl(con, dbplyr::sql(scoped_edges_sql)) |>
    duckplyr::as_duckdb_tibble()

  node_labels <- clust_cc_tables(
    nodes = scoped_nodes,
    edges = scoped_edges,
    from = "from",
    to = "to",
    by = by,
    cluster_col = ".cc_sequence_component",
    node_col = "node_id"
  ) |>
    duckplyr::as_tbl()

  output_sql <- glue::glue_sql(
    "
    WITH member_components AS MATERIALIZED (
      SELECT DISTINCT
        {prefix('m')}m.{`receptor_col`}, l.\".cc_sequence_component\"
      FROM ({scoped_membership}) m
      JOIN ({node_labels}) l
        ON m.node_id = l.node_id{null_join(by, 'm', 'l')}
    ), connected AS (
      SELECT
        {prefix('x')}x.{`receptor_col`},
        min(x.{`receptor_col`}) OVER (
          PARTITION BY {prefix('x')}x.\".cc_sequence_component\"
        ) AS {`cluster_col`},
        count(*) OVER (
          PARTITION BY {prefix('x')}x.\".cc_sequence_component\"
        ) AS {`size_col`}
      FROM member_components x
    ), singletons AS (
      SELECT
        {prefix('u')}u.{`receptor_col`},
        u.{`receptor_col`} AS {`cluster_col`},
        1::BIGINT AS {`size_col`}
      FROM ({universe}) u
      ANTI JOIN member_components x
        ON u.{`receptor_col`} = x.{`receptor_col`}{null_join(by, 'u', 'x')}
    )
    SELECT * FROM connected
    UNION ALL
    SELECT * FROM singletons",
    scoped_membership = DBI::SQL(dbplyr::sql_render(scoped_membership)),
    node_labels = DBI::SQL(dbplyr::sql_render(node_labels)),
    universe = DBI::SQL(dbplyr::sql_render(universe)),
    .con = con
  )

  dplyr::tbl(con, dbplyr::sql(output_sql)) |>
    duckplyr::as_duckdb_tibble() |>
    dplyr::compute(
      name = basename(tempfile(pattern = "immunarch_clust_cc_")),
      prudence = "stingy"
    )
}


# Minimum-label propagation with one retained label per scoped node.
clust_cc_tables <- function(
  nodes,
  edges,
  from,
  to,
  by,
  cluster_col,
  node_col
) {
  by <- if (is.null(by)) character() else by
  required_receptor_cols <- c(by, node_col)
  required_edge_cols <- c(by, from, to)
  node_source <- nodes |>
    dplyr::select(dplyr::all_of(required_receptor_cols)) |>
    dplyr::distinct()
  if (!inherits(node_source, "tbl_sql")) {
    node_source <- duckplyr::as_tbl(node_source)
  }

  edge_source <- edges |>
    dplyr::select(dplyr::all_of(required_edge_cols)) |>
    dplyr::distinct()
  if (!inherits(edge_source, "tbl_sql")) {
    edge_source <- duckplyr::as_tbl(edge_source)
  }

  node_con <- dbplyr::remote_con(node_source)
  edge_con <- dbplyr::remote_con(edge_source)

  if (!identical(node_con, edge_con)) {
    cli::cli_abort(
      "The node and edge tables must use the same DuckDB connection."
    )
  }

  con <- node_con
  group_prefix <- clust_sql_columns(con, by, trailing = TRUE)
  candidate_group_prefix <- clust_sql_columns(
    con, by, "e",
    trailing = TRUE
  )
  label_group_prefix <- clust_sql_columns(con, by, "l", trailing = TRUE)
  grouped_candidate_prefix <- clust_sql_columns(
    con, by, "c",
    trailing = TRUE
  )
  edge_group_join <- clust_sql_join(con, by, "l", "e")
  current_group_join <- clust_sql_join(con, by, "c", "current")
  cte_columns <- clust_sql_columns(
    con, c(by, ".cc_node", ".cc_label")
  )
  key_columns <- clust_sql_columns(con, c(by, ".cc_node"))

  query <- glue::glue_sql(
    "
    WITH RECURSIVE cc_nodes AS MATERIALIZED (
      SELECT DISTINCT {group_prefix}{`node_col`} AS \".cc_node\"
      FROM ({node_source}) node_source
      WHERE {`node_col`} IS NOT NULL
    ), cc_edge_source AS MATERIALIZED (
      SELECT DISTINCT
        {group_prefix}{`from`} AS \".cc_from\", {`to`} AS \".cc_to\"
      FROM ({edge_source}) edge_source
      WHERE {`from`} IS NOT NULL
        AND {`to`} IS NOT NULL
        AND {`from`} <> {`to`}
    ), cc_edges AS MATERIALIZED (
      SELECT
        {group_prefix}\".cc_from\" AS \".cc_src\", \".cc_to\" AS \".cc_dst\"
      FROM cc_edge_source
      UNION
      SELECT
        {group_prefix}\".cc_to\" AS \".cc_src\", \".cc_from\" AS \".cc_dst\"
      FROM cc_edge_source
    ), cc_labels({cte_columns}) USING KEY ({key_columns}) AS (
      SELECT
        {group_prefix}\".cc_node\", \".cc_node\" AS \".cc_label\"
      FROM cc_nodes
      UNION
      (
        SELECT
          {grouped_candidate_prefix}c.\".cc_node\",
          min(c.\".cc_label\") AS \".cc_label\"
        FROM (
          SELECT
            {candidate_group_prefix}e.\".cc_dst\" AS \".cc_node\",
            l.\".cc_label\"
          FROM cc_labels l
          JOIN cc_edges e
            ON l.\".cc_node\" = e.\".cc_src\"{edge_group_join}
        ) c
        JOIN recurring.cc_labels current
          ON c.\".cc_node\" = current.\".cc_node\"{current_group_join}
        GROUP BY {grouped_candidate_prefix}c.\".cc_node\", current.\".cc_label\"
        HAVING min(c.\".cc_label\") < current.\".cc_label\"
      )
    )
    SELECT
      {label_group_prefix}l.\".cc_node\" AS {`node_col`},
      l.\".cc_label\" AS {`cluster_col`}
    FROM cc_labels l",
    node_source = DBI::SQL(dbplyr::sql_render(node_source)),
    edge_source = DBI::SQL(dbplyr::sql_render(edge_source)),
    .con = con
  )

  dplyr::tbl(con, dbplyr::sql(query)) |>
    duckplyr::as_duckdb_tibble() |>
    dplyr::compute(
      name = basename(tempfile(pattern = "immunarch_clust_cc_")),
      prudence = "stingy"
    )
}


#' @description
#' `clust_cc()` finds receptor-level connected components from the contracted
#' sequence graph returned by [dist_hamm()] with `mode = "graph"`.
#'
#' @param edges An `immunarch_sequence_graph` returned by [dist_hamm()] with
#'   `mode = "graph"`. Ordinary receptor edge tables are not accepted.
#' @param by A character vector of column names or `NULL`. These columns come
#'   from the receptor annotations in `idata`; each distinct combination is an
#'   independent clustering scope. The scope may be broader than the comparison
#'   groups used by [dist_hamm()]. The canonical receptor-ID column cannot also
#'   be used in `by`.
#' @param cluster_col A single character string giving the output
#'   component-identifier column name. Within each `by` group, the identifier
#'   is the minimum original receptor identifier in that component. Defaults to
#'   `"cluster_id"`.
#' @param size_col A single character string giving the output component-size column
#'   name. The size is the number of distinct original receptors in the component. Defaults
#'   to `"cluster_size"`.
#'
#' @details
#'
#' Sequence-node distances are applied only where both nodes participate in the
#' same clustering scope. If one receptor belongs to several sequence nodes in
#' a scope, those nodes are joined by identity edges before components are
#' computed. Receptors without a usable sequence remain singleton components.
#'
#' @return A lazy duckplyr table with one row per distinct, non-missing receptor
#'   ID per clustering scope and
#'   the following columns:
#'
#' * the grouping columns supplied through `by`, in the same order;
#' * the canonical receptor identifier;
#' * the component identifier, named by `cluster_col`; and
#' * the number of distinct receptors in that component, named by `size_col`.
#'
#' Call [dplyr::collect()] when an in-memory tibble is needed.
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this example.
#' # Generally, you should NOT do this in your session.
#' db_exec("SET threads TO 1")
#'
#' idata <- immundata::get_test_idata()
#'
#' # Connect receptors whose CDR3 amino-acid sequences differ by at most one
#' # residue. Request the contracted graph used by clust_cc().
#' graph <- dist_hamm(idata, max_dist = 1, mode = "graph")
#' res <- clust_cc(
#'   idata = idata,
#'   edges = graph
#' )
#' res
#' res <- res |> collect()
#'
#' # Build independent graphs within each repertoire and V-gene group. The
#' # same columns must be supplied to by in both functions.
#' group_cols <- c("imd_repertoire_id", "v_call")
#' graph_by_group <- dist_hamm(
#'   idata,
#'   by = group_cols,
#'   min_sim = 0.8,
#'   mode = "graph"
#' )
#'
#' res_by_group <- clust_cc(
#'   idata = idata,
#'   edges = graph_by_group,
#'   by = group_cols
#' )
#'
#' res_by_group |>
#'   dplyr::filter(cluster_size > 1L) |>
#'   dplyr::arrange(dplyr::desc(cluster_size)) |>
#'   dplyr::collect()
#'
#' @rdname clust
#' @concept Clustering
#' @export
clust_cc <- register_method_clust(
  core = clust_cc_impl,
  family = "clust",
  name = "cc",
  need_repertoires = FALSE
)
