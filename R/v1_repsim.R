#' @title Immune repertoire similarity
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to compare **shared receptors** between repertoires.
#' These methods help you measure repertoire similarity, identify public
#' receptors, compare biological groups, and check whether replicate samples
#' have similar receptor composition.
#'
#' ## Available functions
#'
#' The following methods are available.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams repsim_intersection
#' @inheritParams repsim_jaccard
#' @inheritParams repsim_morisita_horn
#' @inheritParams repsim_bray
#' @inheritParams repsim_chao_jaccard
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @section Visualisation:
#' All five `repsim_*()` results can be passed directly to [vis()]. `vis()`
#' displays the result as a dot plot, with repertoire labels on both axes. Each
#' dot represents one pair of repertoires; its size and fill show the value for
#' that pair.
#'
#' ## 1) Shared-receptor counts (`repsim_intersection`)
#'
#' Larger and darker dots mean that the repertoires share more unique receptors.
#' The diagonal shows the total number of unique receptors in each repertoire.
#'
#' ## 2) Jaccard similarity (`repsim_jaccard`)
#'
#' Larger and darker dots mean greater similarity between the sets of unique
#' receptors. Values range from `0` (no shared receptors) to `1` (identical
#' receptor sets).
#'
#' ## 3) Chao-Jaccard similarity (`repsim_chao_jaccard`)
#'
#' Larger and darker dots mean greater estimated similarity after correction
#' for receptors that may be missing because of limited sampling. Values range
#' from `0` (no observed shared receptors) to `1` (maximum estimated
#' similarity).
#'
#' ## 4) Morisita-Horn similarity (`repsim_morisita_horn`)
#'
#' Larger and darker dots mean greater similarity between receptor abundance
#' profiles. Values range from `0` (no overlap) to `1` (identical abundance
#' composition).
#'
#' ## 5) Bray-Curtis dissimilarity (`repsim_bray`)
#'
#' Larger and darker dots mean greater dissimilarity between receptor abundance
#' profiles. Values range from `0` (identical abundance composition) to `1` (no
#' overlap).
#'
#' All repertoire-similarity plots accept:
#'
#' * `size_max_mm` to set the maximum dot diameter.
#' * `row_order` and `col_order` to set the repertoire order on each axis.
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this example.
#' # Generally, you should NOT do this in your session.
#' db_exec("SET threads TO 1")
#'
#' # Load example data.
#' immdata <- get_test_idata()
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

#' @description
#' **1) Shared-receptor counts (`repsim_intersection`).** Count the unique
#' receptors shared by each pair of repertoires. Use this method to find related
#' samples, check replicate similarity, or identify receptors shared across
#' many repertoires.
#'
#' @return
#'
#' ## 1) Shared-receptor counts (`repsim_intersection`)
#' A **symmetric numeric matrix** where rows/columns are repertoire labels derived
#' from `schema_repertoire`. Each cell is the count of shared unique receptors. The
#' diagonal contains per-repertoire richness (total unique receptors).
#'
#' @examples
#' #
#' # Count receptors shared by each pair of repertoires.
#' m_pub <- repsim_intersection(immdata)
#'
#' # Visualise the shared-receptor counts.
#' vis(m_pub)
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


#' @description
#' **2) Jaccard similarity (`repsim_jaccard`).** Compare the sets of unique
#' receptors in each pair of repertoires. The score accounts for the total
#' number of unique receptors in both repertoires. Use it to compare receptor
#' overlap when repertoire sizes differ.
#'
#' \deqn{J(A, B) = \frac{|A \cap B|}{|A \cup B|}}
#'
#' @inheritParams im_common_args
#'
#' @return
#'
#' ## 2) Jaccard similarity (`repsim_jaccard`)
#' A **symmetric numeric matrix** where rows/columns are repertoire labels derived
#' from `schema_repertoire`. Each cell is the Jaccard similarity in `[0, 1]`. The
#' diagonal is `1`.
#'
#' @examples
#' #
#' # Calculate Jaccard similarity.
#' m_jac <- repsim_jaccard(immdata)
#' vis(m_jac)
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_jaccard <- register_immunarch_method(repsim_jaccard_impl, "repsim", "jaccard")


#' @keywords internal
repsim_chao_jaccard_impl <- function(idata) {
  receptor_id_col <- immundata::imd_schema("receptor")
  repertoire_id_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  required_cols <- c(receptor_id_col, repertoire_id_col, count_col)

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

  # Singleton and doubleton status applies after counts for duplicate receptor
  # rows are combined. Keep this reusable projection inside DuckDB.
  receptor_abundances <- idata$annotations |>
    dplyr::select(dplyr::all_of(required_cols)) |>
    dplyr::summarise(
      abundance = sum(!!rlang::sym(count_col)),
      .by = dplyr::all_of(c(receptor_id_col, repertoire_id_col))
    ) |>
    dplyr::filter(.data$abundance > 0) |>
    dplyr::compute()

  rep_x <- paste0(repertoire_id_col, ".x")
  rep_y <- paste0(repertoire_id_col, ".y")
  abundance_x <- "abundance.x"
  abundance_y <- "abundance.y"

  pair_statistics <- receptor_abundances |>
    dplyr::inner_join(
      receptor_abundances,
      by = receptor_id_col,
      suffix = c(".x", ".y")
    ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::summarise(
      shared_x = sum(!!rlang::sym(abundance_x)),
      shared_y = sum(!!rlang::sym(abundance_y)),
      f1_x = sum(dplyr::if_else(!!rlang::sym(abundance_x) == 1, 1L, 0L)),
      f2_x = sum(dplyr::if_else(!!rlang::sym(abundance_x) == 2, 1L, 0L)),
      f1_y = sum(dplyr::if_else(!!rlang::sym(abundance_y) == 1, 1L, 0L)),
      f2_y = sum(dplyr::if_else(!!rlang::sym(abundance_y) == 2, 1L, 0L)),
      shared_x_y1 = sum(dplyr::if_else(
        !!rlang::sym(abundance_y) == 1,
        !!rlang::sym(abundance_x),
        0
      )),
      shared_y_x1 = sum(dplyr::if_else(
        !!rlang::sym(abundance_x) == 1,
        !!rlang::sym(abundance_y),
        0
      )),
      .by = dplyr::all_of(c(rep_x, rep_y))
    )

  repertoire_totals <- repertoire_ids_tbl |>
    dplyr::left_join(
      receptor_abundances |>
        dplyr::summarise(
          total = sum(.data$abundance),
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
    dplyr::left_join(pair_statistics, by = c(rep_x, rep_y)) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c(
          "shared_x", "shared_y", "f1_x", "f2_x", "f1_y", "f2_y",
          "shared_x_y1", "shared_y_x1"
        )),
        ~ dplyr::coalesce(.x, 0)
      ),
      u = .data$shared_x / !!rlang::sym(total_x) +
        ((!!rlang::sym(total_y) - 1) / !!rlang::sym(total_y)) *
          (.data$f1_y / (2 * dplyr::if_else(.data$f2_y > 0, .data$f2_y, 1))) *
          (.data$shared_x_y1 / !!rlang::sym(total_x)),
      v = .data$shared_y / !!rlang::sym(total_y) +
        ((!!rlang::sym(total_x) - 1) / !!rlang::sym(total_x)) *
          (.data$f1_x / (2 * dplyr::if_else(.data$f2_x > 0, .data$f2_x, 1))) *
          (.data$shared_y_x1 / !!rlang::sym(total_y)),
      u = dplyr::if_else(.data$u > 1, 1.0, .data$u),
      v = dplyr::if_else(.data$v > 1, 1.0, .data$v),
      denominator = .data$u + .data$v - .data$u * .data$v,
      value = dplyr::case_when(
        !!rlang::sym(total_x) <= 0 | !!rlang::sym(total_y) <= 0 ~ NA_real_,
        .data$shared_x <= 0 | .data$shared_y <= 0 ~ 0.0,
        TRUE ~ .data$u * .data$v / .data$denominator
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


#' @description
#' **3) Chao-Jaccard similarity (`repsim_chao_jaccard`).** Compare receptor
#' overlap using receptor counts. The method estimates similarity after allowing
#' for shared receptors that may be missing because one or both repertoires were
#' sampled at limited depth. Use it when sampling depth may affect the observed
#' overlap.
#'
#' @return
#'
#' ## 3) Chao-Jaccard similarity (`repsim_chao_jaccard`)
#' A **symmetric numeric matrix** where rows/columns are repertoire labels derived
#' from `schema_repertoire`. Each cell is the Chao-Jaccard similarity in `[0, 1]`
#' (`1` = maximum estimated overlap, `0` = no observed shared receptors). The
#' diagonal is `1`.
#'
#' @references
#' Chao A, Chazdon RL, Colwell RK, Shen TJ (2005). A new statistical approach
#' for assessing similarity of species composition with incidence and abundance
#' data. *Ecology Letters*, 8(2), 148-159. \doi{10.1111/j.1461-0248.2004.00707.x}
#'
#' @examples
#' #
#' # Calculate sampling-corrected Jaccard similarity.
#' m_chao_jac <- repsim_chao_jaccard(immdata)
#' vis(m_chao_jac)
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_chao_jaccard <- register_immunarch_method(
  repsim_chao_jaccard_impl,
  "repsim",
  "chao_jaccard"
)


#' @keywords internal
repsim_morisita_horn_impl <- function(idata) {
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

  # imd_proportion is repeated on barcode-level annotation rows. Keep one
  # abundance per receptor and repertoire, following repsim_bray_impl().
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

  cross_products <- receptor_abundances |>
    dplyr::inner_join(
      receptor_abundances,
      by = receptor_id_col,
      suffix = c(".x", ".y")
    ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::summarise(
      cross_product = sum(
        !!rlang::sym(prop_x) * !!rlang::sym(prop_y)
      ),
      .by = dplyr::all_of(c(rep_x, rep_y))
    )

  repertoire_moments <- repertoire_ids_tbl |>
    dplyr::left_join(
      receptor_abundances |>
        dplyr::summarise(
          total = sum(!!rlang::sym(prop_col)),
          sum_squares = sum(
            !!rlang::sym(prop_col) * !!rlang::sym(prop_col)
          ),
          .by = dplyr::all_of(repertoire_id_col)
        ),
      by = repertoire_id_col
    ) |>
    dplyr::mutate(
      total = dplyr::coalesce(.data$total, 0),
      sum_squares = dplyr::coalesce(.data$sum_squares, 0),
      concentration = dplyr::if_else(
        .data$total > 0,
        .data$sum_squares / (.data$total * .data$total),
        NA_real_
      )
    )

  total_x <- "total.x"
  total_y <- "total.y"
  concentration_x <- "concentration.x"
  concentration_y <- "concentration.y"

  off_diagonal <- dplyr::cross_join(
    repertoire_moments,
    repertoire_moments,
    suffix = c(".x", ".y")
  ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::left_join(cross_products, by = c(rep_x, rep_y)) |>
    dplyr::mutate(
      cross_product = dplyr::coalesce(.data$cross_product, 0),
      denominator =
        (!!rlang::sym(concentration_x) + !!rlang::sym(concentration_y)) *
          !!rlang::sym(total_x) * !!rlang::sym(total_y),
      value = dplyr::if_else(
        .data$denominator > 0,
        2 * .data$cross_product / .data$denominator,
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


#' @description
#' **4) Morisita-Horn similarity (`repsim_morisita_horn`).** Compare receptor
#' abundance profiles. Abundant receptors have more influence on the score than
#' rare receptors. Use it to determine whether the same expanded receptors
#' dominate two repertoires. The method uses `imd_proportion` as the abundance
#' measure.
#'
#' \deqn{C_{MH} = \frac{2 \sum_k x_{ik} x_{jk}}
#' {(\lambda_i + \lambda_j) N_i N_j}}
#'
#' @return
#'
#' ## 4) Morisita-Horn similarity (`repsim_morisita_horn`)
#' A **symmetric numeric matrix** whose rows/columns are repertoire labels derived
#' from `schema_repertoire`, with values in `[0, 1]` (`0` = no shared abundance,
#' `1` = identical abundance composition). The diagonal is `1`. Comparisons
#' involving an empty repertoire are `NA`.
#'
#' @references
#' Horn HS (1966). Measurement of overlap in comparative ecological studies.
#' *The American Naturalist*, 100(914), 419-424. \doi{10.1086/282436}
#'
#' @examples
#' #
#' # Calculate Morisita-Horn similarity.
#' m_mh <- repsim_morisita_horn(immdata)
#' vis(m_mh)
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_morisita_horn <- register_immunarch_method(
  repsim_morisita_horn_impl,
  "repsim",
  "morisita_horn"
)


#' @keywords internal
repsim_bray_impl <- function(
  idata,
  transform = c("none", "log1p")
) {
  transform <- match.arg(transform)

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
    dplyr::mutate(abundance = !!rlang::sym(prop_col)) |>
    dplyr::select(dplyr::all_of(c(
      receptor_id_col,
      repertoire_id_col,
      "abundance"
    )))

  if (transform == "log1p") {
    receptor_abundances <- receptor_abundances |>
      dplyr::mutate(abundance = dd$ln(1 + .data$abundance))
  }

  receptor_abundances <- receptor_abundances |>
    dplyr::compute()

  rep_x <- paste0(repertoire_id_col, ".x")
  rep_y <- paste0(repertoire_id_col, ".y")
  abundance_x <- "abundance.x"
  abundance_y <- "abundance.y"

  overlaps <- receptor_abundances |>
    dplyr::inner_join(
      receptor_abundances,
      by = receptor_id_col,
      suffix = c(".x", ".y")
    ) |>
    dplyr::filter(!!rlang::sym(rep_x) < !!rlang::sym(rep_y)) |>
    dplyr::mutate(
      min_prop = dplyr::if_else(
        !!rlang::sym(abundance_x) <= !!rlang::sym(abundance_y),
        !!rlang::sym(abundance_x),
        !!rlang::sym(abundance_y)
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
          total = sum(.data$abundance),
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

#' @description
#' **5) Bray-Curtis dissimilarity (`repsim_bray`).** Measure how different two
#' receptor abundance profiles are. A value of `0` means the profiles are
#' identical, and a value of `1` means they have no overlap. Use it when you want
#' to express differences directly instead of using a similarity score. The
#' method uses `imd_proportion` as the abundance measure.
#'
#' Set `transform = "log1p"` to reduce the influence of the most abundant
#' receptors and give lower-abundance receptors more influence on the result.
#' This can reveal broader repertoire similarity when a few expanded receptors
#' would otherwise dominate the comparison. The transformation is applied
#' inside the database query and does not modify `idata`.
#'
#' Immune-repertoire studies commonly calculate Bray-Curtis dissimilarity from
#' clonotype proportions without an additional transformation. Therefore,
#' `transform = "none"` remains the default. The `"log1p"` option is an
#' exploratory sensitivity analysis based on a common ecological pre-treatment
#' for reducing the effect of dominant taxa. After transformation, the values no
#' longer sum to one, and the result no longer has the simple interpretation of
#' the fraction of repertoire abundance that does not overlap.
#'
#' @param transform Transformation applied to `imd_proportion` before calculating
#'   dissimilarity. Use `"none"` (the default) for the original proportions or
#'   `"log1p"` for `log(1 + imd_proportion)`.
#'
#' \deqn{BC_{ij} = \frac{\sum_k |a_{ik} - a_{jk}|}
#' {\sum_k (a_{ik} + a_{jk})}}
#' where \eqn{a_{ik}} is the original or transformed abundance of receptor
#' \eqn{k} in repertoire \eqn{i}.
#'
#' @references
#' Clarke KR, Chapman MG, Somerfield PJ, Needham HR (2006). Dispersion-based
#' weighting of species counts in assemblage analyses. *Marine Ecology Progress
#' Series*, 320, 11-27. \doi{10.3354/meps320011}
#'
#' *The TCR assigns naive T cells to a preferred lymph node* (2024).
#' *Science Advances*. \doi{10.1126/sciadv.adl0796}
#'
#' *On the feasibility of using TCR sequencing to follow a vaccination response
#' - lessons learned* (2023). *Frontiers in Immunology*.
#' \doi{10.3389/fimmu.2023.1210168}
#'
#' @return
#'
#' ## 5) Bray-Curtis dissimilarity (`repsim_bray`)
#' A **symmetric numeric matrix** where rows/columns are repertoire labels derived
#' from `schema_repertoire`. Each cell is Bray-Curtis dissimilarity in `[0, 1]`
#' (`0` = identical, `1` = no overlap). The diagonal is `0`.
#'
#' @examples
#' #
#' # Calculate Bray-Curtis dissimilarity.
#' m_bray <- repsim_bray(immdata)
#' vis(m_bray)
#'
#' # Reduce the influence of highly expanded receptors.
#' m_bray_log <- repsim_bray(
#'   immdata,
#'   transform = "log1p"
#' )
#' vis(m_bray_log)
#'
#' @rdname repsim
#' @concept Repertoire similarity
#' @export
repsim_bray <- register_immunarch_method(repsim_bray_impl, "repsim", "bray")
