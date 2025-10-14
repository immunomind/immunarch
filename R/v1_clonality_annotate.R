#' @title Annotate clonality - per-receptor labels for overabundance
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A small family of helpers that **add clonality labels to each receptor** in
#' an [immundata::ImmunData] object.
#'
#' ## Available functions
#' * `annotate_clonality_rank()` - label by **rank bins** within each repertoire.
#' * `annotate_clonality_prop()` - label by **proportion bins** (named thresholds).
#'
#' @param idata An [immundata::ImmunData] object.
#' @inheritParams im_common_args
#'
#' @seealso
#' * Per-repertoire summaries: [airr_clonality]
#' * Data container: [immundata::ImmunData]
#'
#' @examples
#' \dontrun{
#' idata <- get_test_idata() |> agg_repertoires("Therapy")
#' idata_rank <- annotate_clonality_rank(idata)
#' idata_prop <- annotate_clonality_prop(idata)
#' }
#'
#' @name annotate_clonality
#' @concept Clonality
NULL


#' @keywords internal
annotate_clonality_rank_impl <- function(idata,
                                         bins = c(10, 30, 100, 300, 1000, 10000, 100000)) {
  checkmate::check_numeric(bins, lower = 1)

  bins <- sort(bins, decreasing = FALSE)

  clonality_df <- base_clonality_rank(idata = idata, bins = bins)

  by_cols <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  names(by_cols) <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  annotate_immundata(idata,
    clonality_df |>
      select(
        immundata::imd_schema("receptor"),
        immundata::imd_schema("repertoire"),
        clonal_rank_bin
      ),
    by = by_cols
  )
}


#' @description
#' `annotate_clonality_rank()` - for each repertoire, receptors are ordered by
#' within-repertoire abundance (proportion) and assigned a **rank bin** label.
#'
#' @inheritParams airr_clonality_rank
#'
#' @return
#' An [immundata::ImmunData] whose `$annotations` gains:
#' * `clonal_rank_bin` - integer-like label with the applied rank threshold
#'   (outside all thresholds -> `NA`).
#'
#' @rdname annotate_clonality
#' @concept Clonality
#' @export
annotate_clonality_rank <- register_immunarch_method(
  core = annotate_clonality_rank_impl,
  family = "annotate_clonality",
  name = "rank",
  need_repertoires = TRUE
)


#' @keywords internal
annotate_clonality_prop_impl <- function(
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

  by_cols <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  names(by_cols) <- c(immundata::imd_schema("receptor"), immundata::imd_schema("repertoire"))
  annotate_immundata(idata,
    clonality_df |>
      select(
        immundata::imd_schema("receptor"),
        immundata::imd_schema("repertoire"),
        clonal_prop_bin
      ),
    by = by_cols
  )
}


#' @description
#' `annotate_clonality_prop()` - label each receptor by **proportion bin**
#' using named thresholds (matched in descending order; else `"Ultra-rare"`).
#'
#' @inheritParams airr_clonality_prop
#'
#' @return
#' An [immundata::ImmunData] whose `$annotations` gains:
#' * `clonal_prop_bin` - label from `names(bins)` or `"Ultra-rare"`.
#'
#' @rdname annotate_clonality
#' @concept Clonality
#' @export
annotate_clonality_prop <- register_immunarch_method(
  core = annotate_clonality_prop_impl,
  family = "annotate_clonality",
  name = "prop",
  need_repertoires = TRUE
)
