#' Annotate ImmunData with an external sequence database (e.g. VDJdb)
#'
#' Add annotations to ImmunData from an external database such as VDJdb.
#'
#' @param idata An [immundata::ImmunData] object.
#' @param db A data frame, tibble, `duckplyr_df`, or a single character string
#'   containing a path to a delimited or Parquet file. Files are read via
#'   duckplyr.
#' @param by A named character vector of join columns, of the form
#'   `c(idata_col = "db_col")`.
#' @param label_col A single character string naming the column from `db` to
#'   add to `idata` (for VDJdb, this could be `"species"`).
#' @param ... Optional dplyr or duckplyr filter expressions evaluated on
#'   `db` before joining (e.g. `species == "HomoSapiens"`, `gene == "TRB"`).
#'
#' @return A new [immundata::ImmunData] object with the selected database label
#'   added as an annotation column.
#' @export
annotate_with_database <- function(idata, db, by, label_col = "species", ...) {
  checkmate::assert_r6(idata, "ImmunData")
  checkmate::assert_character(by, min.len = 1, names = "named")
  checkmate::assert_string(label_col)

  checkmate::assert(
    checkmate::check_data_frame(db),
    checkmate::check_character(db, len = 1)
  )

  if (is.character(db)) {
    checkmate::assert_file_exists(db)

    ext <- tolower(tools::file_ext(db))
    if (ext %in% c("parquet", "parq")) {
      db_df <- duckplyr::read_parquet_duckdb(db)
    } else {
      db_df <- duckplyr::read_csv_duckdb(db)
    }
  } else {
    # convert to duckplyr? is this necessary?
    db_df <- db
  }

  idata_cols <- names(by)
  db_cols <- unname(by)

  missing_in_idata <- setdiff(idata_cols, colnames(idata$annotations))
  if (length(missing_in_idata)) {
    cli::cli_abort(
      "Column(s) {missing_in_idata} specified in {.arg by} are not found in ImmunData. Please double-check the column names: {.code colnames(idata$annotations)}."
    )
  }

  missing_in_db <- setdiff(db_cols, colnames(db_df))
  if (length(missing_in_db)) {
    cli::cli_abort(
      "Column(s) {missing_in_db} specified in {.arg by} are not found in the database."
    )
  }

  if (!label_col %in% colnames(db_df)) {
    cli::cli_abort(
      "Column '{label_col}' specified as {.arg label_col} is not found in the database."
    )
  }

  filters <- rlang::enquos(...)
  if (length(filters)) {
    db_df <- db_df |>
      filter(!!!filters)
  }

  db_df <- db_df |>
    dplyr::select(dplyr::all_of(unique(c(db_cols, label_col))))

  immundata::annotate(
    idata = idata,
    annotations = db_df,
    by = by
  )
}
