#' @title Clonality — receptor overabundance statistics for immune repertoires
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **receptor overabundance** per repertoire.
#'
#' ## Available functions:
#'
#' @param idata An `ImmunData` object.
#' @inheritParams airr_clonality_line
#' @inheritParams airr_clonality_rank
#' @inheritParams airr_clonality_prop
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @name airr_clonality
#' @concept Clonality
NULL


#' @keywords internal
airr_clonality_line_impl <- function(idata, limit = 100000) {
  checkmate::check_numeric(limit, lower = 10, len = 1)

  n_repertoires <- idata$repertoires |>
    distinct(!!immundata::imd_schema_sym("repertoire")) |>
    pull() |>
    length()

  idata$annotations |>
    select(all_of(c(
      immundata::imd_schema("repertoire"),
      immundata::imd_schema("receptor"),
      immundata::imd_schema("count")
    ))) |>
    distinct(!!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("receptor"),
      .keep_all = TRUE
    ) |>
    arrange(desc(!!immundata::imd_schema_sym("count"))) |>
    slice_head(n = limit * n_repertoires) |> # Optimization before compute - does it make sense, though?
    compute() |> # TODO: If we remove compute, the output breaks. Open an issue in duckplr - something wrong with row_number + mutate-by
    mutate(
      index = row_number(),
      .by = immundata::imd_schema("repertoire")
    ) |>
    filter(index <= limit) |>
    select(-!!immundata::imd_schema_sym("receptor")) |>
    left_join(
      idata$repertoires |> select(-any_of(c(
        immundata::imd_schema("n_barcodes"),
        immundata::imd_schema("n_receptors")
      ))),
      by = immundata::imd_schema("repertoire")
    ) |>
    arrange(index) |>
    collect()
}


#' @description `airr_clonality_line` — build ranked abundance lines: for each
#' repertoire, take the top `limit` receptors by `count` and attach repertoire
#' metadata. Useful for per-repertoire rank–abundance plots.
#'
#' @param limit Positive integer ≥ 10: maximum number of top receptors to keep
#'   **per repertoire** (default `100000`).
#'
#' @return
#'
#' ## `airr_clonality_line`
#' A tibble with columns:
#' * `repertoire_id` — repertoire identifier
#' * `index` — rank within repertoire (1 = most abundant)
#' * `count` — receptor count used for ranking
#' * plus any repertoire metadata columns carried from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_clonality_line
#' #
#' top_line <- airr_clonality_line(immdata, limit = 1000)
#'
#' @rdname airr_clonality
#' @concept Clonality
#' @export
airr_clonality_line <- register_immunarch_method(airr_clonality_line_impl, "airr_clonality", "line")


#' @keywords internal
airr_clonality_rank_impl <- function(idata,
                                     bins = c(10, 30, 100, 300, 1000, 10000, 100000),
                                     output = c("stat", "annot")) {
  checkmate::check_numeric(bins, lower = 1)

  output <- match.arg(output)

  bins <- sort(bins, decreasing = FALSE)

  sql_expr <- paste0(
    "CASE ",
    paste0(map_chr(
      bins,
      ~ cli::format_inline("WHEN ROW_NUMBER() OVER (PARTITION BY {immundata::imd_schema('repertoire')} ORDER BY {immundata::imd_schema('proportion')} DESC) <= {.x} THEN {.x}")
    ), collapse = " "), " ELSE NULL END"
  )

  clonality_df <- idata$annotations |>
    select(all_of(c(
      immundata::imd_schema("repertoire"),
      immundata::imd_schema("receptor"),
      immundata::imd_schema("proportion")
    ))) |>
    distinct(!!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("receptor"),
      .keep_all = TRUE
    ) |>
    arrange() |>
    as_tbl() |>
    mutate(clonal_rank_bin = dbplyr::sql(sql_expr)) |>
    as_duckdb_tibble() |>
    compute()

  if (output == "stat") {
    clonality_df |>
      summarise(
        .by = c(immundata::imd_schema("repertoire"), "clonal_rank_bin"),
        occupied_prop = sum(!!immundata::imd_schema_sym("proportion"), na.rm = TRUE)
      ) |>
      right_join(idata$repertoires, by = immundata::imd_schema("repertoire"))
  } else {
    ImmunData$new(
      schema = idata$schema_receptor,
      annotations = idata$annotations |>
        left_join(clonality_df,
          by = immundata::imd_schema("receptor")
        )
    )
  }
}


#' @description `airr_clonality_rank` — aggregate clonal space by **rank bins**.
#' Receptors are ordered by `proportion` within each repertoire; each receptor
#' is assigned to the smallest threshold in `bins` that contains its rank.
#'
#' @param bins Integer vector of rank thresholds (e.g., `c(10, 100, 1000)`).
#'   For each repertoire, receptors with ranks `<= bins[i]` contribute to bin
#'   `bins[i]`. Bins are sorted ascending internally.
#' @param output One of `"stat"` (default) to return per-repertoire bin
#'   aggregates, or `"annot"` to return an `ImmunData` with an added
#'   `clonal_rank_bin` column in `annotations`.
#'
#' @return
#'
#' ## `airr_clonality_rank`
#' If `output = "stat"`: a tibble with
#' * `repertoire_id`
#' * `clonal_rank_bin` — the rank threshold (e.g., `10`, `100`, …)
#' * `occupied_prop` — sum of `proportion` within the bin
#' * plus repertoire metadata columns from `idata$repertoires`
#'
#' If `output = "annot"`: an `ImmunData` object where `annotations` includes
#' `clonal_rank_bin`.
#'
#' @examples
#' #
#' # airr_clonality_rank
#' #
#' rank_stat <- airr_clonality_rank(immdata, bins = c(10, 100), output = "stat")
#' rank_annot <- airr_clonality_rank(immdata, bins = c(10, 100), output = "annot")
#'
#' @rdname airr_clonality
#' @concept Clonality
#' @export
airr_clonality_rank <- register_immunarch_method(airr_clonality_rank_impl, "airr_clonality", "rank")


#' @keywords internal
airr_clonality_prop_impl <- function(
    idata, bins = c(
      Hyperexpanded = 1e-2,
      Large = 1e-3,
      Medium = 1e-4,
      Small = 1e-5,
      Rare = 1e-6
    ),
    output = c("stat", "annot")) {
  checkmate::check_numeric(bins, lower = 0, min.len = 1)

  output <- match.arg(output)

  bins <- sort(bins, decreasing = TRUE)

  sql_expr <- paste0(
    "CASE ",
    paste0(map2_chr(
      bins, names(bins),
      ~ sprintf("WHEN %s >= %s THEN '%s'", immundata::imd_schema("proportion"), .x, .y)
    ), collapse = " "), " ELSE 'Ultra-rare' END"
  )

  clonality_df <- idata$annotations |>
    select(all_of(c(
      immundata::imd_schema("repertoire"),
      immundata::imd_schema("receptor"),
      immundata::imd_schema("proportion")
    ))) |>
    distinct(!!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("receptor"),
      .keep_all = TRUE
    ) |>
    duckplyr::as_tbl() |>
    mutate(clonal_prop_bin = dbplyr::sql(sql_expr)) |>
    duckplyr::as_duckdb_tibble() |>
    compute()

  if (output == "stat") {
    clonality_df |>
      summarise(
        .by = c(immundata::imd_schema("repertoire"), "clonal_prop_bin"),
        occupied_prop = sum(!!immundata::imd_schema_sym("proportion"), na.rm = TRUE)
      ) |>
      right_join(idata$repertoires, by = immundata::imd_schema("repertoire"))
  } else {
    ImmunData$new(
      schema = idata$schema_receptor,
      annotations = idata$annotations |>
        left_join(clonality_df,
          by = immundata::imd_schema("receptor")
        )
    )
  }
}


#' @description `airr_clonality_prop` — aggregate clonal space by **proportion bins**.
#' Each receptor is assigned to a named bin according to its `proportion`
#' (e.g., `Hyperexpanded ≥ 1e-2`, `Large ≥ 1e-3`, …). Thresholds are matched in
#' descending order; unmatched receptors fall into `"Ultra-rare"`.
#'
#' @param bins A **named** numeric vector of thresholds (e.g.,
#'   `c(Hyperexpanded = 1e-2, Large = 1e-3, ...)`). Names become bin labels and
#'   must be non-empty. Internally sorted in descending order.
#' @param output One of `"stat"` (default) to return per-repertoire bin
#'   aggregates, or `"annot"` to return an `ImmunData` with an added
#'   `clonal_prop_bin` column in `annotations`.
#'
#' @return
#'
#' ## `airr_clonality_prop`
#' If `output = "stat"`: a tibble with
#' * `repertoire_id`
#' * `clonal_prop_bin` — factor-like label from `names(bins)` or `"Ultra-rare"`
#' * `occupied_prop` — sum of `proportion` within the bin
#' * plus repertoire metadata columns from `idata$repertoires`
#'
#' If `output = "annot"`: an `ImmunData` object where `annotations` includes
#' `clonal_prop_bin`.
#'
#' @examples
#' #
#' # airr_clonality_prop
#' #
#' prop_stat <- airr_clonality_prop(immdata, output = "stat")
#' prop_annot <- airr_clonality_prop(immdata, output = "annot")
#'
#' @rdname airr_clonality
#' @concept Clonality
#' @export
airr_clonality_prop <- register_immunarch_method(airr_clonality_prop_impl, "airr_clonality", "prop")
