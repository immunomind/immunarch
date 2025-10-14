#' @title Clonality - receptor overabundance statistics for immune repertoires
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **receptor overabundance** per repertoire. Helps in deciphering the structure and partition the repertoire.
#'
#' ## Available functions
#'
#' Supported methods are the following.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams airr_clonality_line
#' @inheritParams airr_clonality_rank
#' @inheritParams airr_clonality_prop
#' @inheritParams im_common_args
#'
#' @seealso
#' * Per-repertoire summaries: [annotate_clonality]
#' * Data container: [immundata::ImmunData]
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this session.
#' # Change this only if you know what you're doing (e.g., multi-user machines, shared CI/servers).
#' db_exec("SET threads TO 1")
#'
#' # Load data
#' \dontrun{
#' immdata <- get_test_idata() |> agg_repertoires("Therapy")
#' }
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
    collect() |> # TODO: .by doesn't work in slice_head in duckplyr. What to do instead then?
    slice_head(n = limit * n_repertoires, by = !!immundata::imd_schema_sym("repertoire")) |>
    mutate(
      index = row_number(),
      .by = immundata::imd_schema("repertoire")
    ) |>
    select(-!!immundata::imd_schema_sym("receptor")) |>
    arrange(index)
}


#' @description `airr_clonality_line` - build ranked abundance lines: for each
#' repertoire, take the top `limit` receptors by `count` and attach repertoire
#' metadata. Useful for per-repertoire rank-abundance plots.
#'
#' @param limit Positive integer >= 10: maximum number of top receptors to keep
#'   **per repertoire** (default `100000`).
#'
#' @return
#'
#' ## `airr_clonality_line`
#' A tibble with columns:
#' * `repertoire_id` - repertoire identifier
#' * `index` - rank within repertoire (1 = most abundant)
#' * `count` - receptor count used for ranking
#' * plus any repertoire metadata columns carried from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_clonality_line
#' #
#' \dontrun{
#' top_line <- airr_clonality_line(immdata, limit = 1000)
#' }
#'
#' @rdname airr_clonality
#' @concept Clonality
#' @export
airr_clonality_line <- register_immunarch_method(
  core = airr_clonality_line_impl,
  family = "airr_clonality",
  name = "line",
  need_repertoires = TRUE
)


#' @keywords internal
airr_clonality_rank_impl <- function(idata,
                                     bins = c(10, 30, 100, 300, 1000, 10000, 100000)) {
  checkmate::check_numeric(bins, lower = 1)

  bins <- sort(bins, decreasing = FALSE)

  clonality_df <- base_clonality_rank(idata = idata, bins = bins)

  clonality_df |>
    summarise(
      .by = c(immundata::imd_schema("repertoire"), "clonal_rank_bin"),
      occupied_prop = sum(!!immundata::imd_schema_sym("proportion"), na.rm = TRUE)
    )
}


#' @description `airr_clonality_rank` - aggregate clonal space by **rank bins**.
#' Receptors are ordered by `proportion` within each repertoire; each receptor
#' is assigned to the smallest threshold in `bins` that contains its rank.
#'
#' @param bins Integer vector of rank thresholds (e.g., `c(10, 100, 1000)`).
#'   For each repertoire, receptors with ranks `<= bins[i]` contribute to bin
#'   `bins[i]`. Bins are sorted ascending internally.
#'
#' @return
#'
#' ## `airr_clonality_rank`
#' A tibble with
#' * `repertoire_id`
#' * `clonal_rank_bin` - the rank threshold (e.g., `10`, `100`, ...)
#' * `occupied_prop` - sum of `proportion` within the bin
#' * plus repertoire metadata columns from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_clonality_rank
#' #
#' \dontrun{
#' rank_stat <- airr_clonality_rank(immdata, bins = c(10, 100))
#' }
#'
#' @rdname airr_clonality
#' @concept Clonality
#' @export
airr_clonality_rank <- register_immunarch_method(
  core = airr_clonality_rank_impl,
  family = "airr_clonality",
  name = "rank",
  need_repertoires = TRUE
)


#' @keywords internal
airr_clonality_prop_impl <- function(
    idata, bins = c(
      Hyperexpanded = 1e-2,
      Large = 1e-3,
      Medium = 1e-4,
      Small = 1e-5,
      Rare = 1e-6
    )) {
  checkmate::check_numeric(bins, lower = 0, min.len = 1)

  bins <- sort(bins, decreasing = TRUE)

  clonality_df <- base_clonality_prop(idata = idata, bins = bins)

  clonality_df |>
    summarise(
      .by = c(immundata::imd_schema("repertoire"), "clonal_prop_bin"),
      occupied_prop = sum(!!immundata::imd_schema_sym("proportion"), na.rm = TRUE)
    )
}


#' @description `airr_clonality_prop` - aggregate clonal space by **proportion bins**.
#' Each receptor is assigned to a named bin according to its `proportion`
#' (e.g., `Hyperexpanded >= 1e-2`, `Large >= 1e-3`, ...). Thresholds are matched in
#' descending order; unmatched receptors fall into `"Ultra-rare"`.
#'
#' @param bins A **named** numeric vector of thresholds (e.g.,
#'   `c(Hyperexpanded = 1e-2, Large = 1e-3, ...)`). Names become bin labels and
#'   must be non-empty. Internally sorted in descending order.
#'
#' @return
#'
#' ## `airr_clonality_prop`
#' A tibble with
#' * `repertoire_id`
#' * `clonal_prop_bin` - factor-like label from `names(bins)` or `"Ultra-rare"`
#' * `occupied_prop` - sum of `proportion` within the bin
#' * plus repertoire metadata columns from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_clonality_prop
#' #
#' \dontrun{
#' prop_stat <- airr_clonality_prop(immdata)
#' }
#'
#' @rdname airr_clonality
#' @concept Clonality
#' @export
airr_clonality_prop <- register_immunarch_method(
  core = airr_clonality_prop_impl,
  family = "airr_clonality",
  name = "prop",
  need_repertoires = TRUE
)
