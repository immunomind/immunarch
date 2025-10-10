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
    unique()

  result_matrix <- matrix(-1, nrow = length(repertoire_ids), ncol = length(repertoire_ids))

  for (i in seq_along(repertoire_ids[-length(repertoire_ids)])) {
    for (j in seq_along(repertoire_ids[(i + 1):length(repertoire_ids)])) {
      rep_1_index <- repertoire_ids[i]
      rep_2_index <- repertoire_ids[i + j]

      val_pub <- idata$annotations |>
        filter(!!rlang::sym(repertoire_id_col) == rep_1_index, !!rlang::sym(repertoire_id_col) == rep_2_index) |>
        distinct(!!rlang::sym(repertoire_id_col)) |>
        count() |>
        pull("n")

      result_matrix[rep_1_index, rep_2_index] <- val_pub
      result_matrix[rep_2_index, rep_1_index] <- val_pub
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
#' m_pub <- airr_public_intersection(immdata)
#'
#' @rdname airr_public
#' @concept Public indices
#' @export
airr_public_intersection <- register_immunarch_method(airr_public_intersection_impl, "airr_public", "intersection")


#' @keywords internal
airr_public_jaccard_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  repertoire_ids <- idata$repertoires |>
    pull({{ repertoire_id_col }}) |>
    unique()

  result_matrix <- matrix(-1, nrow = length(repertoire_ids), ncol = length(repertoire_ids))

  for (i in seq_along(repertoire_ids[-length(repertoire_ids)])) {
    for (j in seq_along(repertoire_ids[(i + 1):length(repertoire_ids)])) {
      rep_1_index <- repertoire_ids[i]
      rep_2_index <- repertoire_ids[i + j]

      repertoire_pair <- idata$annotations |>
        filter(!!rlang::sym(repertoire_id_col) %in% c(rep_1_index, rep_2_index)) |>
        select({{ repertoire_id_col }}, {{ receptor_id_col }})

      rep_1 <- repertoire_pair |>
        filter(!!rlang::sym(repertoire_id_col) == rep_1_index) |>
        select({{ receptor_id_col }})
      rep_2 <- repertoire_pair |>
        filter(!!rlang::sym(repertoire_id_col) == rep_2_index) |>
        select({{ receptor_id_col }})

      val_inter <- intersect(rep_1, rep_2) |>
        count() |>
        pull("n")
      val_union <- union(rep_1, rep_2) |>
        count() |>
        pull("n")

      result_matrix[rep_1_index, rep_2_index] <- val_inter / val_union
      result_matrix[rep_2_index, rep_1_index] <- val_inter / val_union
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
#' m_jac <- airr_public_jaccard(immdata)
#'
#' @rdname airr_public
#' @concept Public indices
#' @export
airr_public_jaccard <- register_immunarch_method(airr_public_jaccard_impl, "airr_public", "jaccard")
