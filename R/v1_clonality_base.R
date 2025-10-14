#' @keywords internal
base_clonality_rank <- function(idata, bins) {
  checkmate::check_numeric(bins, lower = 1)

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

  clonality_df
}


#' @keywords internal
base_clonality_prop <- function(idata, bins) {
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

  clonality_df
}
