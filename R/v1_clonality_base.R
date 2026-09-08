#' @keywords internal
base_clonality_rank <- function(idata, bins) {
  checkmate::assert_integerish(bins, lower = 1, min.len = 1, any.missing = FALSE)

  bins <- sort(bins, decreasing = FALSE)
  repertoire_col <- immundata::imd_schema("repertoire")
  proportion_col <- immundata::imd_schema("proportion")

  source <- idata$annotations |>
    select(all_of(c(
      repertoire_col,
      immundata::imd_schema("receptor"),
      proportion_col
    ))) |>
    distinct(!!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("receptor"),
      .keep_all = TRUE
    ) |>
    arrange() |>
    as_tbl()
  con <- dbplyr::remote_con(source)

  conditions <- glue::glue_sql(
    "WHEN ROW_NUMBER() OVER (
      PARTITION BY {`repertoire_col`} ORDER BY {`proportion_col`} DESC
    ) <= {bins} THEN {bins}",
    .con = con
  )
  conditions <- glue::glue_sql_collapse(conditions, sep = " ")
  sql_expr <- glue::glue_sql(
    "CASE {conditions} ELSE NULL END", .con = con
  )

  clonality_df <- source |>
    mutate(clonal_rank_bin = dbplyr::sql(sql_expr)) |>
    as_duckdb_tibble() |>
    compute()

  clonality_df
}


#' @keywords internal
base_clonality_prop <- function(idata, bins) {
  repertoire_col <- immundata::imd_schema("repertoire")
  proportion_col <- immundata::imd_schema("proportion")

  source <- idata$annotations |>
    select(all_of(c(
      repertoire_col,
      immundata::imd_schema("receptor"),
      proportion_col
    ))) |>
    distinct(!!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("receptor"),
      .keep_all = TRUE
    ) |>
    duckplyr::as_tbl()
  con <- dbplyr::remote_con(source)

  thresholds <- unname(bins)
  labels <- names(bins)
  conditions <- glue::glue_sql(
    "WHEN {`proportion_col`} >= {thresholds} THEN {labels}",
    .con = con
  )
  conditions <- glue::glue_sql_collapse(conditions, sep = " ")
  fallback <- "Ultra-rare"
  sql_expr <- glue::glue_sql(
    "CASE {conditions} ELSE {fallback} END", .con = con
  )

  clonality_df <- source |>
    mutate(clonal_prop_bin = dbplyr::sql(sql_expr)) |>
    duckplyr::as_duckdb_tibble() |>
    compute()

  clonality_df
}
