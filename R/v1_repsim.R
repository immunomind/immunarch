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
  repertoire_ids <- idata$repertoires |>
    pull({{ repertoire_id_col }}) |>
    unique() |>
    sort()

  rep_labels <- if (is.null(idata$schema_repertoire) || length(idata$schema_repertoire) == 0) {
    tibble::tibble(rep_id = repertoire_ids, label = as.character(repertoire_ids))
  } else {
    idata$repertoires |>
      dplyr::select(dplyr::all_of(c(repertoire_id_col, idata$schema_repertoire))) |>
      dplyr::distinct(.data[[repertoire_id_col]], .keep_all = TRUE) |>
      dplyr::arrange(.data[[repertoire_id_col]]) |>
      tidyr::unite(".label", dplyr::all_of(idata$schema_repertoire), sep = "|", na.rm = TRUE) |>
      dplyr::transmute(rep_id = .data[[repertoire_id_col]], label = .data$.label)
  }

  result_matrix <- matrix(NA_real_, nrow = length(repertoire_ids), ncol = length(repertoire_ids), dimnames = list(rep_labels$label, rep_labels$label))

  target_cols <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  set_tbl <- idata$annotations |>
    dplyr::select(dplyr::all_of(target_cols)) |>
    dplyr::distinct(
      !!rlang::sym(receptor_id_col),
      !!rlang::sym(repertoire_id_col),
      .keep_all = TRUE
    )
  rep_x <- paste0(immundata::imd_schema("repertoire"), ".x")
  rep_y <- paste0(immundata::imd_schema("repertoire"), ".y")
  pairs <- set_tbl |>
    dplyr::inner_join(set_tbl, by = receptor_id_col) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::summarise(n = dplyr::n(), .by = dplyr::all_of(c(rep_x, rep_y))) |>
    dplyr::collect()

  rep_sizes <- set_tbl |>
    dplyr::summarise(n = dplyr::n(), .by = dplyr::all_of(repertoire_id_col)) |>
    dplyr::collect()
  size_map <- stats::setNames(rep_sizes$n, as.character(rep_sizes[[repertoire_id_col]]))

  if (length(repertoire_ids) >= 2L) {
    for (rep_i in 1:(length(repertoire_ids) - 1)) {
      for (rep_j in (rep_i + 1):length(repertoire_ids)) {
        id_i <- repertoire_ids[[rep_i]]
        id_j <- repertoire_ids[[rep_j]]

        val <- pairs |>
          filter(!!rlang::sym(rep_x) == id_i, !!rlang::sym(rep_y) == id_j) |>
          pull(n)

        val <- if (length(val) == 0) 0 else val

        result_matrix[rep_i, rep_j] <- val
        result_matrix[rep_j, rep_i] <- val
      }
    }
  }

  diag_vals <- vapply(
    repertoire_ids,
    function(x) {
      val <- size_map[[as.character(x)]]
      if (is.null(val) || length(val) == 0 || is.na(val)) 0 else as.numeric(val)
    },
    numeric(1)
  )
  diag(result_matrix) <- diag_vals

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
repsim_intersection <- register_immunarch_method(repsim_intersection_impl, "repsim", "intersection", )


#' @keywords internal
repsim_jaccard_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  repertoire_ids <- idata$repertoires |>
    pull({{ repertoire_id_col }}) |>
    unique() |>
    sort()

  rep_labels <- if (is.null(idata$schema_repertoire) || length(idata$schema_repertoire) == 0) {
    tibble::tibble(rep_id = repertoire_ids, label = as.character(repertoire_ids))
  } else {
    idata$repertoires |>
      dplyr::select(dplyr::all_of(c(repertoire_id_col, idata$schema_repertoire))) |>
      dplyr::distinct(.data[[repertoire_id_col]], .keep_all = TRUE) |>
      dplyr::arrange(.data[[repertoire_id_col]]) |>
      tidyr::unite(".label", dplyr::all_of(idata$schema_repertoire), sep = "|", na.rm = TRUE) |>
      dplyr::transmute(rep_id = .data[[repertoire_id_col]], label = .data$.label)
  }

  result_matrix <- matrix(NA_real_, nrow = length(repertoire_ids), ncol = length(repertoire_ids), dimnames = list(rep_labels$label, rep_labels$label))

  target_cols <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  set_tbl <- idata$annotations |>
    dplyr::select(dplyr::all_of(target_cols)) |>
    dplyr::distinct(
      !!rlang::sym(receptor_id_col),
      !!rlang::sym(repertoire_id_col),
      .keep_all = TRUE
    )
  rep_x <- paste0(immundata::imd_schema("repertoire"), ".x")
  rep_y <- paste0(immundata::imd_schema("repertoire"), ".y")
  pairs <- set_tbl |>
    dplyr::inner_join(set_tbl, by = receptor_id_col) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::summarise(n = dplyr::n(), .by = dplyr::all_of(c(rep_x, rep_y))) |>
    dplyr::collect()

  rep_sizes <- set_tbl |>
    dplyr::summarise(n = dplyr::n(), .by = dplyr::all_of(repertoire_id_col)) |>
    dplyr::collect()
  size_map <- stats::setNames(rep_sizes$n, as.character(rep_sizes[[repertoire_id_col]]))

  if (length(repertoire_ids) >= 2L) {
    for (rep_i in 1:(length(repertoire_ids) - 1)) {
      for (rep_j in (rep_i + 1):length(repertoire_ids)) {
        id_i <- repertoire_ids[[rep_i]]
        id_j <- repertoire_ids[[rep_j]]

        inter_val <- pairs |>
          filter(!!rlang::sym(rep_x) == id_i, !!rlang::sym(rep_y) == id_j) |>
          pull(n)
        inter_val <- if (length(inter_val) == 0) 0 else inter_val

        size_i <- size_map[[as.character(id_i)]]
        size_i <- if (is.null(size_i) || length(size_i) == 0 || !is.finite(size_i)) 0 else size_i

        size_j <- size_map[[as.character(id_j)]]
        size_j <- if (is.null(size_j) || length(size_j) == 0 || !is.finite(size_j)) 0 else size_j

        union_val <- size_i + size_j - inter_val

        val <- ifelse(union_val > 0, inter_val / union_val, NA_real_)

        result_matrix[rep_i, rep_j] <- val
        result_matrix[rep_j, rep_i] <- val
      }
    }
  }
  diag(result_matrix) <- 1

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

  if (!(prop_col %in% colnames(idata$annotations))) {
    cli::cli_abort(
      paste0(
        "Column {.field ", prop_col, "} is required in {.field idata$annotations}. ",
        "Run {.code immundata::agg_repertoires()} to compute repertoire proportions."
      )
    )
  }

  repertoire_ids <- idata$repertoires |>
    pull({{ repertoire_id_col }}) |>
    unique() |>
    sort()

  rep_labels <- if (is.null(idata$schema_repertoire) || length(idata$schema_repertoire) == 0) {
    tibble::tibble(rep_id = repertoire_ids, label = as.character(repertoire_ids))
  } else {
    idata$repertoires |>
      dplyr::select(dplyr::all_of(c(repertoire_id_col, idata$schema_repertoire))) |>
      dplyr::distinct(.data[[repertoire_id_col]], .keep_all = TRUE) |>
      dplyr::arrange(.data[[repertoire_id_col]]) |>
      tidyr::unite(".label", dplyr::all_of(idata$schema_repertoire), sep = "|", na.rm = TRUE) |>
      dplyr::transmute(rep_id = .data[[repertoire_id_col]], label = .data$.label)
  }

  result_matrix <- matrix(NA_real_, nrow = length(repertoire_ids), ncol = length(repertoire_ids), dimnames = list(rep_labels$label, rep_labels$label))

  target_cols <- c(receptor_id_col, repertoire_id_col, prop_col)
  rep_x <- paste0(repertoire_id_col, ".x")
  rep_y <- paste0(repertoire_id_col, ".y")
  prop_x <- paste0(prop_col, ".x")
  prop_y <- paste0(prop_col, ".y")

  prop_tbl <- idata$annotations |>
    dplyr::select(dplyr::all_of(target_cols)) |>
    dplyr::distinct(
      !!rlang::sym(receptor_id_col),
      !!rlang::sym(repertoire_id_col),
      .keep_all = TRUE
    )

  pairs <- prop_tbl |>
    dplyr::inner_join(prop_tbl, by = receptor_id_col) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::mutate(
      min_prop = dplyr::if_else(
        !!rlang::sym(prop_x) <= !!rlang::sym(prop_y),
        !!rlang::sym(prop_x),
        !!rlang::sym(prop_y)
      )
    ) |>
    dplyr::summarise(intersection = sum(.data$min_prop), .by = dplyr::all_of(c(rep_x, rep_y))) |>
    dplyr::collect()

  rep_totals <- prop_tbl |>
    dplyr::summarise(sum_prop = sum(!!rlang::sym(prop_col)), .by = dplyr::all_of(repertoire_id_col)) |>
    dplyr::collect()

  if (length(repertoire_ids) >= 2L) {
    for (rep_i in 1:(length(repertoire_ids) - 1)) {
      for (rep_j in (rep_i + 1):length(repertoire_ids)) {
        id_i <- repertoire_ids[[rep_i]]
        id_j <- repertoire_ids[[rep_j]]

        inter_val <- pairs |>
          dplyr::filter(
            !!rlang::sym(rep_x) == id_i,
            !!rlang::sym(rep_y) == id_j
          ) |>
          dplyr::pull(intersection)
        inter_val <- if (length(inter_val) == 0) 0 else inter_val

        size_i <- rep_totals |>
          dplyr::filter(!!rlang::sym(repertoire_id_col) == id_i) |>
          dplyr::pull(sum_prop)
        size_i <- if (length(size_i) == 0 || !is.finite(size_i)) 0 else size_i

        size_j <- rep_totals |>
          dplyr::filter(!!rlang::sym(repertoire_id_col) == id_j) |>
          dplyr::pull(sum_prop)
        size_j <- if (length(size_j) == 0 || !is.finite(size_j)) 0 else size_j

        denom <- size_i + size_j
        val <- ifelse(denom > 0, (denom - 2 * inter_val) / denom, NA_real_)

        result_matrix[rep_i, rep_j] <- val
        result_matrix[rep_j, rep_i] <- val
      }
    }
  }

  diag(result_matrix) <- 0
  result_matrix
}


#' @description `repsim_bray` - **Bray-Curtis dissimilarity** between repertoires:
#' \eqn{\sum_k |p_{ik} - p_{jk}| / \sum_k (p_{ik} + p_{jk})}, where \eqn{p_{ik}}
#' are receptor abundances in repertoire `i`.
#'
#' The method uses `imd_proportion` values from `idata$annotations`. This avoids
#' constructing full receptor-by-repertoire matrices in RAM and only collects
#' repertoire-level summaries and pairwise overlaps.
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
