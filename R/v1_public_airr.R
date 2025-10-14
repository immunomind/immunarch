#' @title Public indices - pairwise repertoire overlap
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
#' @inheritParams airr_public_intersection
#' @inheritParams airr_public_jaccard
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this session.
#' # Change this only if you know what you're doing (e.g., multi-user machines, shared CI/servers).
#' db_exec("SET threads TO 1")
#' # Load data
#' immdata <- get_test_idata() |> agg_repertoires("Therapy")
#'
#' @name airr_public
#' @concept Public indices
NULL


#' @keywords internal
airr_public_intersection_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  repertoire_ids <- idata$repertoires |>
    pull({{ repertoire_id_col }}) |>
    unique() |>
    sort()

  rep_labels <- idata$repertoires |>
    dplyr::select(dplyr::all_of(c(repertoire_id_col, idata$schema_repertoire))) |>
    tidyr::unite(".label", dplyr::all_of(idata$schema_repertoire), sep = "|", na.rm = TRUE) |>
    dplyr::transmute(rep_id = .data[[repertoire_id_col]], label = .data$.label)

  result_matrix <- matrix(NA, nrow = length(repertoire_ids), ncol = length(repertoire_ids), dimnames = list(rep_labels$label, rep_labels$label))

  target_cols <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  rep_x <- paste0(immundata::imd_schema("repertoire"), ".x")
  rep_y <- paste0(immundata::imd_schema("repertoire"), ".y")
  pairs <- idata$annotations |>
    select(all_of(target_cols)) |>
    inner_join(idata$annotations |> select(all_of(target_cols)), by = immundata::imd_schema("receptor")) |>
    filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    summarise(n = n(), .by = all_of(c(rep_x, rep_y))) |>
    collect()

  # yes-yes, I know it is inefficient
  for (rep_i in 1:(length(repertoire_ids) - 1)) {
    for (rep_j in (rep_i + 1):length(repertoire_ids)) {
      val <- pairs |>
        filter(!!rlang::sym(rep_x) == rep_i, !!rlang::sym(rep_y) == rep_j) |>
        pull(n)

      val <- if (length(val) == 0) 0 else val

      result_matrix[rep_i, rep_j] <- val
      result_matrix[rep_j, rep_i] <- val
    }
  }

  result_matrix
}

#' @description `airr_public_intersection` - number of **shared receptors** between
#' each pair of repertoires (intersection size). Handy for quick overlap heatmaps,
#' QC of replicate similarity, or spotting donor-shared "public" clonotypes.
#'
#' @return
#'
#' ## `airr_public_intersection`
#' A **symmetric numeric matrix** where rows/columns are `repertoire_id` and each
#' cell is the count of shared unique receptors. The diagonal contains per-repertoire
#' richness (total unique receptors). Row/column names are repertoire IDs.
#'
#' @examples
#' #
#' # airr_public_intersection
#' #
#' \dontrun{
#' m_pub <- airr_public_intersection(immdata)
#' }
#'
#' @rdname airr_public
#' @concept Public indices
#' @export
airr_public_intersection <- register_immunarch_method(airr_public_intersection_impl, "airr_public", "intersection", )


#' @keywords internal
airr_public_jaccard_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  repertoire_ids <- idata$repertoires |>
    pull({{ repertoire_id_col }}) |>
    unique() |>
    sort()

  rep_labels <- idata$repertoires |>
    dplyr::select(dplyr::all_of(c(repertoire_id_col, idata$schema_repertoire))) |>
    tidyr::unite(".label", dplyr::all_of(idata$schema_repertoire), sep = "|", na.rm = TRUE) |>
    dplyr::transmute(rep_id = .data[[repertoire_id_col]], label = .data$.label)

  result_matrix <- matrix(NA, nrow = length(repertoire_ids), ncol = length(repertoire_ids), dimnames = list(rep_labels$label, rep_labels$label))

  target_cols <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  rep_x <- paste0(immundata::imd_schema("repertoire"), ".x")
  rep_y <- paste0(immundata::imd_schema("repertoire"), ".y")
  pairs <- idata$annotations |>
    select(all_of(target_cols)) |>
    inner_join(idata$annotations |> select(all_of(target_cols)), by = immundata::imd_schema("receptor")) |>
    filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    summarise(n = n(), .by = all_of(c(rep_x, rep_y))) |>
    collect()

  # yes-yes, I know it is inefficient
  for (rep_i in 1:(length(repertoire_ids) - 1)) {
    for (rep_j in (rep_i + 1):length(repertoire_ids)) {
      inter_val <- pairs |>
        filter(!!rlang::sym(rep_x) == rep_i, !!rlang::sym(rep_y) == rep_j) |>
        pull(n)
      inter_val <- if (length(inter_val) == 0) 0 else inter_val

      size_i <- idata$repertoires |>
        filter(!!rlang::sym(repertoire_id_col) == rep_i) |>
        pull(immundata::imd_schema("n_receptors"))
      size_j <- idata$repertoires |>
        filter(!!rlang::sym(repertoire_id_col) == rep_j) |>
        pull(immundata::imd_schema("n_receptors"))

      union_val <- size_i + size_j - inter_val

      val <- ifelse(union_val > 0, inter_val / union_val, NA_real_)

      result_matrix[rep_i, rep_j] <- val
      result_matrix[rep_j, rep_i] <- val
    }
  }

  result_matrix
}


#' @description `airr_public_jaccard` - **Jaccard similarity** of receptor
#' sets between repertoires (\eqn{A \cap B}{A cap B} / \eqn{A \cup B}{A cup B}). Best when comparing cohorts with
#' different sizes to get a scale-invariant overlap score.
#'
#' @inheritParams im_common_args
#'
#' @return
#'
#' ## `airr_public_jaccard`
#' A **symmetric numeric matrix** where rows/columns are `repertoire_id` and each
#' cell is the Jaccard similarity in `[0, 1]`. The diagonal is `1`. Row/column
#' names are repertoire IDs.
#'
#' @examples
#' #
#' # airr_public_jaccard
#' #
#' \dontrun{
#' m_jac <- airr_public_jaccard(immdata)
#' }
#'
#' @rdname airr_public
#' @concept Public indices
#' @export
airr_public_jaccard <- register_immunarch_method(airr_public_jaccard_impl, "airr_public", "jaccard")
