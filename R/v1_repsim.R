#' @title Immune repertoire similarity
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **public or shared receptors** between repertoire.
#'
#' ## Available functions
#'
#' Supported methods are the following.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams repsim_intersection
#' @inheritParams repsim_jaccard
#' @inheritParams repsim_bray
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this session.
#' # Change this only if you know what you're doing (e.g., multi-user machines, shared CI/servers).
#' db_exec("SET threads TO 1")
#' # Load data
#' \dontrun{
#' immdata <- get_test_idata() |> agg_repertoires("Therapy")
#' }
#'
#' @name repsim
#' @concept Repertoire similarity
NULL


#' @keywords internal
repsim_intersection_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  required_cols <- c(receptor_id_col, repertoire_id_col)

  repertoire_meta <- idata$repertoires |>
    dplyr::collect()
  repertoire_ids <- repertoire_meta[[repertoire_id_col]]

  repertoire_labels <- repertoire_meta |>
    tidyr::unite(
      ".label",
      dplyr::all_of(idata$schema_repertoire),
      sep = "|",
      na.rm = TRUE
    ) |>
    dplyr::pull(.data$.label)

  result_matrix <- matrix(
    0,
    nrow = length(repertoire_ids),
    ncol = length(repertoire_ids),
    dimnames = list(repertoire_labels, repertoire_labels)
  )
  if (length(repertoire_ids) == 0L) {
    return(result_matrix)
  }

  repertoire_ids_tbl <- idata$repertoires |>
    dplyr::select(dplyr::all_of(repertoire_id_col))

  receptor_sets <- idata$annotations |>
    dplyr::select(dplyr::all_of(required_cols)) |>
    dplyr::distinct(
      !!rlang::sym(receptor_id_col),
      !!rlang::sym(repertoire_id_col)
    )

  shared_receptors <- receptor_sets |>
    dplyr::summarise(
      .n_repertoires = dplyr::n(),
      .by = dplyr::all_of(receptor_id_col)
    ) |>
    dplyr::filter(.data$.n_repertoires > 1L) |>
    dplyr::select(dplyr::all_of(receptor_id_col))

  shared_receptor_sets <- receptor_sets |>
    dplyr::inner_join(shared_receptors, by = receptor_id_col)

  rep_x <- paste0(repertoire_id_col, ".x")
  rep_y <- paste0(repertoire_id_col, ".y")

  overlaps <- shared_receptor_sets |>
    dplyr::inner_join(
      shared_receptor_sets,
      by = receptor_id_col,
      suffix = c(".x", ".y")
    ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::summarise(
      value = dplyr::n(),
      .by = dplyr::all_of(c(rep_x, rep_y))
    )

  off_diagonal <- dplyr::cross_join(
    repertoire_ids_tbl,
    repertoire_ids_tbl,
    suffix = c(".x", ".y")
  ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::left_join(overlaps, by = c(rep_x, rep_y)) |>
    dplyr::mutate(value = dplyr::coalesce(.data$value, 0L))

  diagonal <- repertoire_ids_tbl |>
    dplyr::left_join(
      receptor_sets |>
        dplyr::summarise(
          value = dplyr::n(),
          .by = dplyr::all_of(repertoire_id_col)
        ),
      by = repertoire_id_col
    ) |>
    dplyr::transmute(
      !!rlang::sym(rep_x) := .data[[repertoire_id_col]],
      !!rlang::sym(rep_y) := .data[[repertoire_id_col]],
      value = dplyr::coalesce(.data$value, 0L)
    )

  pair_values <- dplyr::union_all(off_diagonal, diagonal) |>
    dplyr::collect()

  row_index <- match(pair_values[[rep_x]], repertoire_ids)
  col_index <- match(pair_values[[rep_y]], repertoire_ids)
  result_matrix[cbind(row_index, col_index)] <- pair_values$value
  result_matrix[cbind(col_index, row_index)] <- pair_values$value

  result_matrix
}

#' @description `repsim_intersection` - number of **shared receptors** between
#' each pair of repertoires (intersection size). Handy for quick overlap heatmaps,
#' QC of replicate similarity, or spotting donor-shared "public" clonotypes.
#'
#' @return
#'
#' ## `repsim_intersection`
#' A **symmetric numeric matrix** where rows/columns are `repertoire_id` and each
#' cell is the count of shared unique receptors. The diagonal contains per-repertoire
#' richness (total unique receptors). Row/column names are repertoire IDs.
#'
#' @examples
#' #
#' # repsim_intersection
#' #
#' \dontrun{
#' m_pub <- repsim_intersection(immdata)
#' }
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_intersection <- register_immunarch_method(repsim_intersection_impl, "repsim", "intersection")


#' @keywords internal
repsim_jaccard_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  required_cols <- c(receptor_id_col, repertoire_id_col)

  repertoire_meta <- idata$repertoires |>
    dplyr::collect()
  repertoire_ids <- repertoire_meta[[repertoire_id_col]]

  repertoire_labels <- repertoire_meta |>
    tidyr::unite(
      ".label",
      dplyr::all_of(idata$schema_repertoire),
      sep = "|",
      na.rm = TRUE
    ) |>
    dplyr::pull(.data$.label)

  result_matrix <- matrix(
    NA_real_,
    nrow = length(repertoire_ids),
    ncol = length(repertoire_ids),
    dimnames = list(repertoire_labels, repertoire_labels)
  )
  if (length(repertoire_ids) == 0L) {
    return(result_matrix)
  }

  repertoire_ids_tbl <- idata$repertoires |>
    dplyr::select(dplyr::all_of(repertoire_id_col))

  receptor_sets <- idata$annotations |>
    dplyr::select(dplyr::all_of(required_cols)) |>
    dplyr::distinct(
      !!rlang::sym(receptor_id_col),
      !!rlang::sym(repertoire_id_col)
    )

  shared_receptors <- receptor_sets |>
    dplyr::summarise(
      .n_repertoires = dplyr::n(),
      .by = dplyr::all_of(receptor_id_col)
    ) |>
    dplyr::filter(.data$.n_repertoires > 1L) |>
    dplyr::select(dplyr::all_of(receptor_id_col))

  shared_receptor_sets <- receptor_sets |>
    dplyr::inner_join(shared_receptors, by = receptor_id_col)

  rep_x <- paste0(repertoire_id_col, ".x")
  rep_y <- paste0(repertoire_id_col, ".y")

  overlaps <- shared_receptor_sets |>
    dplyr::inner_join(
      shared_receptor_sets,
      by = receptor_id_col,
      suffix = c(".x", ".y")
    ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::summarise(
      intersection = dplyr::n(),
      .by = dplyr::all_of(c(rep_x, rep_y))
    )

  repertoire_sizes <- repertoire_ids_tbl |>
    dplyr::left_join(
      receptor_sets |>
        dplyr::summarise(
          size = dplyr::n(),
          .by = dplyr::all_of(repertoire_id_col)
        ),
      by = repertoire_id_col
    ) |>
    dplyr::mutate(size = dplyr::coalesce(.data$size, 0L))

  size_x <- "size.x"
  size_y <- "size.y"
  off_diagonal <- dplyr::cross_join(
    repertoire_sizes,
    repertoire_sizes,
    suffix = c(".x", ".y")
  ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::left_join(overlaps, by = c(rep_x, rep_y)) |>
    dplyr::mutate(
      intersection = dplyr::coalesce(.data$intersection, 0L),
      union = !!rlang::sym(size_x) + !!rlang::sym(size_y) - .data$intersection,
      value = dplyr::if_else(
        .data$union > 0L,
        .data$intersection / .data$union,
        NA_real_
      )
    ) |>
    dplyr::select(dplyr::all_of(c(rep_x, rep_y, "value")))

  diagonal <- repertoire_ids_tbl |>
    dplyr::transmute(
      !!rlang::sym(rep_x) := .data[[repertoire_id_col]],
      !!rlang::sym(rep_y) := .data[[repertoire_id_col]],
      value = 1.0
    )

  pair_values <- dplyr::union_all(off_diagonal, diagonal) |>
    dplyr::collect()

  row_index <- match(pair_values[[rep_x]], repertoire_ids)
  col_index <- match(pair_values[[rep_y]], repertoire_ids)
  result_matrix[cbind(row_index, col_index)] <- pair_values$value
  result_matrix[cbind(col_index, row_index)] <- pair_values$value

  result_matrix
}


#' @description `repsim_jaccard` - **Jaccard similarity** of receptor
#' sets between repertoires (\eqn{A \cap B}{A cap B} / \eqn{A \cup B}{A cup B}). Best when comparing cohorts with
#' different sizes to get a scale-invariant overlap score.
#'
#' @inheritParams im_common_args
#'
#' @return
#'
#' ## `repsim_jaccard`
#' A **symmetric numeric matrix** where rows/columns are `repertoire_id` and each
#' cell is the Jaccard similarity in `[0, 1]`. The diagonal is `1`. Row/column
#' names are repertoire IDs.
#'
#' @examples
#' #
#' # repsim_jaccard
#' #
#' \dontrun{
#' m_jac <- repsim_jaccard(immdata)
#' }
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_jaccard <- register_immunarch_method(repsim_jaccard_impl, "repsim", "jaccard")


#' @keywords internal
repsim_bray_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  prop_col <- immundata::imd_schema("proportion")
  required_cols <- c(receptor_id_col, repertoire_id_col, prop_col)

  repertoire_meta <- idata$repertoires |>
    dplyr::collect()
  repertoire_ids <- repertoire_meta[[repertoire_id_col]]

  repertoire_labels <- repertoire_meta |>
    tidyr::unite(
      ".label",
      dplyr::all_of(idata$schema_repertoire),
      sep = "|",
      na.rm = TRUE
    ) |>
    dplyr::pull(.data$.label)

  result_matrix <- matrix(
    NA_real_,
    nrow = length(repertoire_ids),
    ncol = length(repertoire_ids),
    dimnames = list(repertoire_labels, repertoire_labels)
  )
  if (length(repertoire_ids) == 0L) {
    return(result_matrix)
  }

  repertoire_ids_tbl <- idata$repertoires |>
    dplyr::select(dplyr::all_of(repertoire_id_col))

  # Keep the reusable projection inside DuckDB; no receptor rows enter R.
  receptor_abundances <- idata$annotations |>
    dplyr::select(dplyr::all_of(required_cols)) |>
    dplyr::distinct(
      !!rlang::sym(receptor_id_col),
      !!rlang::sym(repertoire_id_col),
      .keep_all = TRUE
    ) |>
    dplyr::compute()

  rep_x <- paste0(repertoire_id_col, ".x")
  rep_y <- paste0(repertoire_id_col, ".y")
  prop_x <- paste0(prop_col, ".x")
  prop_y <- paste0(prop_col, ".y")

  overlaps <- receptor_abundances |>
    dplyr::inner_join(
      receptor_abundances,
      by = receptor_id_col,
      suffix = c(".x", ".y")
    ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::mutate(
      min_prop = dplyr::if_else(
        !!rlang::sym(prop_x) <= !!rlang::sym(prop_y),
        !!rlang::sym(prop_x),
        !!rlang::sym(prop_y)
      )
    ) |>
    dplyr::summarise(
      intersection = sum(.data$min_prop),
      .by = dplyr::all_of(c(rep_x, rep_y))
    )

  repertoire_totals <- repertoire_ids_tbl |>
    dplyr::left_join(
      receptor_abundances |>
        dplyr::summarise(
          total = sum(!!rlang::sym(prop_col)),
          .by = dplyr::all_of(repertoire_id_col)
        ),
      by = repertoire_id_col
    ) |>
    dplyr::mutate(total = dplyr::coalesce(.data$total, 0))

  total_x <- "total.x"
  total_y <- "total.y"
  off_diagonal <- dplyr::cross_join(
    repertoire_totals,
    repertoire_totals,
    suffix = c(".x", ".y")
  ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::left_join(overlaps, by = c(rep_x, rep_y)) |>
    dplyr::mutate(
      intersection = dplyr::coalesce(.data$intersection, 0),
      denominator = !!rlang::sym(total_x) + !!rlang::sym(total_y),
      value = dplyr::if_else(
        .data$denominator > 0,
        (.data$denominator - 2 * .data$intersection) / .data$denominator,
        NA_real_
      )
    ) |>
    dplyr::select(dplyr::all_of(c(rep_x, rep_y, "value")))

  diagonal <- repertoire_ids_tbl |>
    dplyr::transmute(
      !!rlang::sym(rep_x) := .data[[repertoire_id_col]],
      !!rlang::sym(rep_y) := .data[[repertoire_id_col]],
      value = 0.0
    )

  pair_values <- dplyr::union_all(off_diagonal, diagonal) |>
    dplyr::collect()

  row_index <- match(pair_values[[rep_x]], repertoire_ids)
  col_index <- match(pair_values[[rep_y]], repertoire_ids)
  result_matrix[cbind(row_index, col_index)] <- pair_values$value
  result_matrix[cbind(col_index, row_index)] <- pair_values$value

  result_matrix
}

#' @description `repsim_bray` - **Bray-Curtis dissimilarity** between repertoires:
#' \eqn{\sum_k |p_{ik} - p_{jk}| / \sum_k (p_{ik} + p_{jk})}, where \eqn{p_{ik}}
#' are receptor abundances in repertoire `i`.
#'
#' The method uses `imd_proportion` values from `idata$annotations`. The reusable
#' receptor-abundance projection is materialized inside DuckDB; only repertoire
#' metadata and the final dense pair table are collected into R.
#'
#' @return
#'
#' ## `repsim_bray`
#' A **symmetric numeric matrix** where rows/columns are repertoire labels and
#' each cell is Bray-Curtis dissimilarity in `[0, 1]` (`0` = identical, `1` = no overlap).
#' The diagonal is `0`.
#'
#' @examples
#' #
#' # repsim_bray
#' #
#' \dontrun{
#' m_bray <- repsim_bray(immdata)
#'
#' # Optional transforms can be done on proportions beforehand:
#' immdata_log <- immundata::mutate(immdata, imd_proportion = log1p(imd_proportion))
#' m_bray_log <- repsim_bray(immdata_log)
#' }
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_bray <- register_immunarch_method(repsim_bray_impl, "repsim", "bray")
