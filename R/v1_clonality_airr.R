#' @title Clonality - receptor overabundance statistics for immune repertoires
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **receptor overabundance** per repertoire.
#' It helps in deciphering the structure and partition the repertoire.
#' Higher clonality means that one or several receptors are over-expanded. Such effects
#' is usually associated with the ongoing immune system response or anomalous behaviour, e.g.,
#' in blood cancers.
#'
#' ## Available functions
#'
#' The following methods are available.
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
#' @section Visualisation:
#' All three `airr_clonality_*()` results can be passed directly to [vis()].
#' With `autojoin = TRUE` (the default), repertoire metadata is included in the
#' result and can be selected by the plotting arguments.
#'
#' ## 1) Rank-abundance lines (`airr_clonality_line`)
#'
#' `vis()` plots receptor rank (`index`) against
#' `imd_proportion` on a logarithmic y-axis by default. The line plot accepts:
#'
#' * `yval` selects either `"imd_proportion"` (the depth-normalised default) or
#'   `"imd_count"` (raw observed counts).
#' * `color` colours the lines by a repertoire or metadata column. Lines remain
#'   grouped by repertoire when a metadata column is selected.
#' * `log = FALSE` switches to a linear y-axis.
#' * `facet` splits the plot by one column, or creates a facet grid when given
#'   two columns; `dir = "h"` or `dir = "v"` controls the wrapping direction.
#' * `title` replaces the default plot title.
#'
#' ## 2) Rank-bin statistics (`airr_clonality_rank`)
#'
#' `vis()` produces a stacked column plot of
#' `occupied_prop` per repertoire, filled by `clonal_rank_bin`. Receptors beyond
#' the largest requested rank have a missing bin and are not drawn, so the total
#' bar height shows the repertoire space occupied by the displayed ranks.
#'
#' ## 3) Proportion-bin statistics (`airr_clonality_prop`)
#'
#' `vis()` produces a stacked column plot of
#' `occupied_prop` per repertoire, filled by `clonal_prop_bin`. The
#' `"Ultra-rare"` bin is included, so each complete repertoire normally sums to
#' 100%.
#'
#' Both stacked plots accept `xval`, `yval`, and `fill` to remap columns, plus
#' `facet`, `dir`, and `title` for layout and labeling.
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this example.
#' # Generally, you should NOT do this in your session.
#' db_exec("SET threads TO 1")
#'
#' # Load example data.
#' immdata <- get_test_idata()
#'
#' @name airr_clonality
#' @concept Clonality
NULL


#' @keywords internal
airr_clonality_line_impl <- function(idata, limit = 100000) {
  checkmate::assert_integerish(
    limit,
    lower = 10,
    len = 1,
    any.missing = FALSE
  )

  idata$annotations |>
    select(all_of(c(
      immundata::imd_schema("repertoire"),
      immundata::imd_schema("receptor"),
      immundata::imd_schema("count"),
      immundata::imd_schema("proportion")
    ))) |>
    distinct(!!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("receptor"),
      .keep_all = TRUE
    ) |>
    arrange(desc(!!immundata::imd_schema_sym("count"))) |>
    collect() |> # TODO: .by doesn't work in slice_head in duckplyr. What to do instead then?
    slice_head(n = limit, by = !!immundata::imd_schema_sym("repertoire")) |>
    mutate(
      index = row_number(),
      .by = immundata::imd_schema("repertoire")
    ) |>
    select(-!!immundata::imd_schema_sym("receptor")) |>
    arrange(index)
}


#' @description
#' **1) Rank-abundance lines (`airr_clonality_line`).** Prepare a plot for each
#' repertoire. It orders receptors from most to least abundant and keeps the top
#' `limit` receptors. Use it to see whether a small number of receptors dominate
#' the repertoire: a steep line suggests strong clonal expansion.
#'
#' @param limit Positive integer >= 10: maximum number of top receptors to keep
#'   **per repertoire** (default `100000`).
#'
#' @return
#'
#' ## 1) Rank-abundance lines (`airr_clonality_line`)
#' A tibble with columns:
#' * `imd_repertoire_id` - repertoire identifier
#' * `index` - rank within repertoire (1 = most abundant)
#' * `imd_count` - receptor count used for ranking
#' * `imd_proportion` - receptor proportion within the repertoire
#' * plus any repertoire metadata columns carried from `idata$repertoires`
#'
#' @examples
#' #
#' # Analyse receptor abundance by rank.
#' top_line <- airr_clonality_line(immdata, limit = 1000)
#'
#' # Visualise the depth-normalised rank-abundance curve.
#' vis(top_line)
#'
#' # Or visualise raw receptor counts.
#' vis(top_line, yval = "imd_count")
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
  checkmate::assert_integerish(bins, lower = 1, min.len = 1, any.missing = FALSE)

  bins <- sort(bins, decreasing = FALSE)

  clonality_df <- base_clonality_rank(idata = idata, bins = bins)

  clonality_df |>
    summarise(
      .by = c(immundata::imd_schema("repertoire"), "clonal_rank_bin"),
      occupied_prop = sum(!!immundata::imd_schema_sym("proportion"), na.rm = TRUE)
    )
}


#' @description
#' **2) Rank-bin statistics (`airr_clonality_rank`).** Summarise how much of a repertoire is
#' occupied by receptors in different rank ranges. Receptors are ordered by
#' proportion and grouped using rank thresholds such as the top 10, top 100, or
#' top 1,000 receptors. Use it to compare whether expansion is concentrated in
#' the highest-ranked receptors across repertoires.
#'
#' @param bins Integer vector of rank thresholds (e.g., `c(10, 100, 1000)`).
#'   For each repertoire, receptors with ranks `<= bins[i]` contribute to bin
#'   `bins[i]`. Bins are sorted ascending internally.
#'
#' @return
#'
#' ## 2) Rank-bin statistics (`airr_clonality_rank`)
#' A tibble with
#' * `repertoire_id`
#' * `clonal_rank_bin` - the rank threshold (e.g., `10`, `100`, ...)
#' * `occupied_prop` - sum of `proportion` within the bin
#' * plus repertoire metadata columns from `idata$repertoires`
#'
#' @examples
#' #
#' # Analyse the repertoire space occupied by different rank ranges.
#' rank_stat <- airr_clonality_rank(immdata, bins = c(10, 100))
#'
#' # Visualise the rank-bin summary.
#' vis(rank_stat)
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
  )
) {
  checkmate::assert_numeric(
    bins,
    lower = 0,
    min.len = 1,
    any.missing = FALSE,
    finite = TRUE,
    names = "named"
  )

  bins <- sort(bins, decreasing = TRUE)

  clonality_df <- base_clonality_prop(idata = idata, bins = bins)

  clonality_df |>
    summarise(
      .by = c(immundata::imd_schema("repertoire"), "clonal_prop_bin"),
      occupied_prop = sum(!!immundata::imd_schema_sym("proportion"), na.rm = TRUE)
    )
}


#' @description
#' **3) Proportion-bin statistics (`airr_clonality_prop`).** Summarise how much of a repertoire is
#' occupied by receptors at different abundance levels. It assigns each receptor
#' to a proportion-based group, such as `"Hyperexpanded"`, `"Large"`, or
#' `"Ultra-rare"`. Each bin contains receptors with the proportion in the specified limits --
#' just like in a histogram. Use it to compare the balance of expanded and rare receptors
#' across repertoires.
#'
#' @param bins A **named** numeric vector of thresholds (e.g.,
#'   `c(Hyperexpanded = 1e-2, Large = 1e-3, ...)`). Names become bin labels and
#'   must be non-empty. Internally sorted in descending order.
#'
#' @return
#'
#' ## 3) Proportion-bin statistics (`airr_clonality_prop`)
#' A tibble with
#' * `repertoire_id`
#' * `clonal_prop_bin` - factor-like label from `names(bins)` or `"Ultra-rare"`
#' * `occupied_prop` - sum of `proportion` within the bin
#' * plus repertoire metadata columns from `idata$repertoires`
#'
#' @examples
#' #
#' # Analyse the repertoire space occupied by abundance groups.
#' prop_stat <- airr_clonality_prop(immdata)
#'
#' # Visualise the abundance-group summary.
#' vis(prop_stat)
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
