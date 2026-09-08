#' @title Add Immunoglobulin germline information
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' Annotate clonal families with germline sequences.
#'
#' @name annotate_ig_germlines
NULL

#' @keywords internal
annotate_ig_germlines_impl <- function(
  idata,
  germline_col = "germline_alignment",
  v_end_col = "v_alignment_end",
  j_start_col = "j_alignment_start",
  d_mask_col = "germline_alignment_d_mask"
) {
  germline <- rlang::sym(germline_col)
  v_end <- rlang::sym(v_end_col)
  j_start <- rlang::sym(j_start_col)
  idata |>
    dplyr::mutate(
      !!d_mask_col := dplyr::if_else(
        is.na(!!germline) | is.na(!!v_end) | is.na(!!j_start),
        NA_character_,
        dd$concat(
          dd$substr(!!germline, 1L, !!v_end),
          # AIRR coordinates are commonly read as BIGINT, while rpad requires INTEGER.
          dd$rpad("", as.integer(!!j_start - !!v_end - 1L), "N"),
          dd$substr(!!germline, !!j_start)
        )
      )
    )
}


#' @rdname annotate_ig_germlines
#' @export
annotate_ig_germlines <- register_method_annotate(
  core = annotate_ig_germlines_impl,
  family = "annotate_ig",
  name = "germlines",
  required_cols = c("germline_col", "v_end_col", "j_start_col"),
  need_repertoires = FALSE
)
