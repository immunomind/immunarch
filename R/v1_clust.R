# Clustering receptor graphs.


#' Connected-component clustering for receptor graphs
#'
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions for clustering sparse receptor graphs represented by
#' an `ImmunData` object and a lazy edge table.
#'
#' @param idata An [immundata::ImmunData] object. Its annotations provide the
#'   graph nodes through the canonical receptor-ID column and, when requested,
#'   repertoire metadata.
#' @inheritParams im_common_args
#'
#' @name clust
#' @concept Clustering
NULL


#' @keywords internal
clust_cc_impl <- function(
  idata,
  edges,
  from = "imd_receptor_id_1",
  to = "imd_receptor_id_2",
  by = NULL,
  cluster_col = "cluster_id",
  size_col = "cluster_size"
) {
  checkmate::assert_string(from)
  checkmate::assert_string(to)
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

  quote_ident <- function(con, x) {
    as.character(dbplyr::sql_escape_ident(con, x))
  }

  group_select <- function(quoted_by, table_alias = NULL) {
    if (!length(quoted_by)) {
      return("")
    }

    prefix <- if (is.null(table_alias)) "" else paste0(table_alias, ".")
    paste0(prefix, quoted_by, collapse = ", ")
  }

  group_join <- function(
    quoted_by,
    left_alias,
    right_alias,
    prefix = " AND "
  ) {
    if (!length(quoted_by)) {
      return("")
    }

    paste0(
      prefix,
      paste0(
        left_alias,
        ".",
        quoted_by,
        " IS NOT DISTINCT FROM ",
        right_alias,
        ".",
        quoted_by,
        collapse = " AND "
      )
    )
  }

  required_receptor_cols <- c(by, receptor_col)
  required_edge_cols <- c(by, from, to)
  missing_receptor_cols <- setdiff(
    required_receptor_cols,
    colnames(idata$annotations)
  )
  missing_edge_cols <- setdiff(required_edge_cols, colnames(edges))

  if (length(missing_receptor_cols)) {
    cli::cli_abort(
      "Column(s) missing from {.field idata$annotations}: {.field {missing_receptor_cols}}."
    )
  }

  if (length(missing_edge_cols)) {
    cli::cli_abort(
      "Column(s) missing from {.arg edges}: {.field {missing_edge_cols}}."
    )
  }

  node_source <- idata$annotations |>
    dplyr::select(dplyr::all_of(required_receptor_cols)) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()

  edge_source <- edges |>
    dplyr::select(dplyr::all_of(required_edge_cols)) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()

  node_con <- dbplyr::remote_con(node_source)
  edge_con <- dbplyr::remote_con(edge_source)

  if (!identical(node_con, edge_con)) {
    cli::cli_abort(
      "{.field idata$annotations} and {.arg edges} must use the same DuckDB connection."
    )
  }

  con <- node_con
  q_by <- vapply(by, quote_ident, character(1), con = con)
  q_node <- quote_ident(con, receptor_col)
  q_from <- quote_ident(con, from)
  q_to <- quote_ident(con, to)
  q_cluster <- quote_ident(con, cluster_col)
  q_size <- quote_ident(con, size_col)

  node_sql <- as.character(dbplyr::sql_render(node_source))
  edge_sql <- as.character(dbplyr::sql_render(edge_source))

  selected_groups <- group_select(q_by)
  group_prefix <- if (length(q_by)) paste0(selected_groups, ", ") else ""
  group_r <- group_select(q_by, "r")
  group_l <- group_select(q_by, "l")
  group_n <- group_select(q_by, "n")
  reach_group_prefix <- if (length(q_by)) paste0(group_r, ", ") else ""
  label_group_prefix <- if (length(q_by)) paste0(group_l, ", ") else ""
  node_group_prefix <- if (length(q_by)) paste0(group_n, ", ") else ""

  reach_group_join <- group_join(q_by, "r", "e")
  destination_group_join <- group_join(q_by, "e", "n")
  partition_cols <- paste0(
    if (length(q_by)) paste0(group_l, ", ") else "",
    "l.\".cc_component\""
  )

  query <- paste0(
    "WITH RECURSIVE ",
    "cc_nodes AS MATERIALIZED (",
    " SELECT DISTINCT ",
    group_prefix,
    q_node,
    " AS \".cc_node\"",
    " FROM (",
    node_sql,
    ") node_source",
    " WHERE ",
    q_node,
    " IS NOT NULL",
    "), ",
    "cc_edge_source AS MATERIALIZED (",
    " SELECT DISTINCT ",
    group_prefix,
    q_from,
    " AS \".cc_from\", ",
    q_to,
    " AS \".cc_to\"",
    " FROM (",
    edge_sql,
    ") edge_source",
    " WHERE ",
    q_from,
    " IS NOT NULL",
    " AND ",
    q_to,
    " IS NOT NULL",
    " AND ",
    q_from,
    " <> ",
    q_to,
    "), ",
    "cc_edges AS (",
    " SELECT ",
    group_prefix,
    "\".cc_from\" AS \".cc_src\", ",
    "\".cc_to\" AS \".cc_dst\"",
    " FROM cc_edge_source",
    " UNION ",
    " SELECT ",
    group_prefix,
    "\".cc_to\" AS \".cc_src\", ",
    "\".cc_from\" AS \".cc_dst\"",
    " FROM cc_edge_source",
    "), ",
    "cc_reach AS (",
    " SELECT ",
    group_prefix,
    "\".cc_node\", \".cc_node\" AS \".cc_label\"",
    " FROM cc_nodes",
    " UNION ",
    " SELECT ",
    node_group_prefix,
    "n.\".cc_node\", r.\".cc_label\"",
    " FROM cc_reach r",
    " JOIN cc_edges e",
    " ON r.\".cc_node\" = e.\".cc_src\"",
    reach_group_join,
    " JOIN cc_nodes n",
    " ON e.\".cc_dst\" = n.\".cc_node\"",
    destination_group_join,
    "), ",
    "cc_labels AS (",
    " SELECT ",
    reach_group_prefix,
    "r.\".cc_node\", min(r.\".cc_label\") AS \".cc_component\"",
    " FROM cc_reach r",
    " GROUP BY ",
    reach_group_prefix,
    "r.\".cc_node\"",
    ") ",
    "SELECT ",
    label_group_prefix,
    "l.\".cc_node\" AS ",
    q_node,
    ", l.\".cc_component\" AS ",
    q_cluster,
    ", count(*) OVER (PARTITION BY ",
    partition_cols,
    ") AS ",
    q_size,
    " FROM cc_labels l"
  )

  dplyr::tbl(con, dbplyr::sql(query)) |>
    duckplyr::as_duckdb_tibble()
}


#' @description
#' `clust_cc()` finds connected components in an undirected receptor graph. The
#' nodes come from the canonical receptor-ID column in `idata$annotations`, and
#' `edges` supplies the receptor pairs. Consequently, every receptor is
#' returned, including receptors with no edges as one-node components.
#'
#' @param edges A lazy duckplyr table with one row per edge. It must contain the
#'   columns named by `from` and `to`, plus every column in `by`, and must use
#'   the same DuckDB connection as `idata$annotations`. Endpoint values must use
#'   the same type and identifiers as the canonical receptor-ID column.
#'   Edges are treated as undirected; duplicate edges, reversed duplicates,
#'   self-loops, and rows with a missing endpoint do not change the result.
#' @param from,to Single strings naming the endpoint columns in `edges`. The
#'   defaults, `"imd_receptor_id_1"` and `"imd_receptor_id_2"`, match the output
#'   of [dist_hamm()].
#' @param by `NULL` or a character vector of column names present in both
#'   `idata$annotations` and `edges`. Each distinct combination defines an
#'   independent graph, so edges never connect nodes across groups. For example,
#'   `by = "imd_repertoire_id"` clusters each repertoire separately. The
#'   canonical receptor-ID column cannot also be used in `by`.
#' @param cluster_col A single string giving the output component-identifier
#'   column name. Within each `by` group, the identifier is the minimum node
#'   identifier in that component. Defaults to `"cluster_id"`.
#' @param size_col A single string giving the output component-size column
#'   name. The size is the number of distinct nodes in the component. Defaults
#'   to `"cluster_size"`.
#'
#' @details
#' Only receptor identifiers listed in `idata$annotations` are returned. Edges
#' whose endpoints are not present there cannot introduce new nodes. Duplicate
#' annotation rows are collapsed, and rows with a missing receptor identifier
#' are omitted. The order of output rows is not guaranteed; use
#' [dplyr::arrange()] before collecting when order matters.
#'
#' [dist_hamm()] produces the expected edge columns using the default `from`
#' and `to` names. Its distance columns (`seq_len`, `dist`, `norm_dist`, and
#' `sim`) may remain in `edges`; `clust_cc()` selects only the endpoint and `by`
#' columns it needs. Because [dist_hamm()] derives its edges from the same
#' `idata`, the node and edge inputs automatically share a DuckDB connection.
#'
#' @return A lazy duckplyr table with one row per distinct, non-missing node and
#'   the following columns:
#'
#' * the grouping columns supplied through `by`, in the same order;
#' * the canonical receptor identifier;
#' * the component identifier, named by `cluster_col`; and
#' * the number of nodes in that component, named by `size_col`.
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
#' # residue. dist_hamm() and clust_cc() use matching receptor-ID defaults.
#' edges <- dist_hamm(idata, max_dist = 1, autojoin = FALSE)
#' res <- clust_cc(
#'   idata = idata,
#'   edges = edges
#' )
#' res
#' res <- res |> collect()
#'
#' # Build independent graphs within each repertoire and V-gene group. The
#' # same columns must be supplied to by in both functions.
#' group_cols <- c("imd_repertoire_id", "v_call")
#' edges_by_group <- dist_hamm(
#'   idata,
#'   by = group_cols,
#'   min_sim = 0.8
#' )
#'
#' res_by_group <- clust_cc(
#'   idata = idata,
#'   edges = edges_by_group,
#'   by = group_cols
#' )
#'
#' res_by_group |>
#'   dplyr::filter(.data$cluster_size > 1L) |>
#'   dplyr::arrange(dplyr::desc(.data$cluster_size)) |>
#'   dplyr::collect()
#'
#' @rdname clust
#' @concept Clustering
#' @export
clust_cc <- register_immunarch_method(
  core = clust_cc_impl,
  family = "clust",
  name = "cc",
  need_repertoires = FALSE
)
