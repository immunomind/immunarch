#' @keywords internal
resolve_argument_by <- function(idata, by) {
  if (is.null(by)) {
    return(character())
  }

  annotation_cols <- colnames(idata$annotations)

  if (length(by) == 1L && is.na(by)) {
    locus_col <- immundata::imd_schema("locus")
    return(if (locus_col %in% annotation_cols) locus_col else character())
  }

  checkmate::assert_character(by, min.len = 1L, any.missing = FALSE, unique = TRUE)

  missing_by <- setdiff(by, annotation_cols)
  if (length(missing_by)) {
    cli::cli_abort(
      "Grouping column(s) [{.field {missing_by}}] are missing from {.field idata$annotations}."
    )
  }

  by
}
