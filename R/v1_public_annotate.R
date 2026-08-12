#' Annotate receptors with publicness statistics
#'
#' `r lifecycle::badge("experimental")`
#'
#' Adds per-receptor publicness metrics to an [immundata::ImmunData] object.
#' Global metrics are always computed. If `$repertoires` contains `imd_strata_id`
#' (created by `immundata::agg_strata()`), strata-specific metrics are added as well.
#' If strata is missing, the function reports that only global metrics are produced.
#'
#' For large datasets, this can be a computationally heavy operation.
#' After annotation, consider persisting the result as a snapshot with
#' [immundata::write_immundata()] to avoid recomputation.
#'
#' @param idata An [immundata::ImmunData] object with aggregated repertoires.
#'
#' @return
#' An [immundata::ImmunData] with added annotation columns:
#' * global metrics:
#'   `imd_public_incidence`, `imd_public_incidence_prop`,
#'   `imd_public_count_min`, `imd_public_count_max`,
#'   `imd_public_count_mean`, `imd_public_count_median`,
#'   `imd_public_prop_min`, `imd_public_prop_max`,
#'   `imd_public_prop_mean`, `imd_public_prop_median`;
#' * if `imd_strata_id` exists in `$repertoires`, per-strata columns with suffix
#'   `_strata_{id}` for the same metric family.
#'
#' @examples
#' \dontrun{
#' idata_global <- get_test_immundata() |>
#'   agg_repertoires(c("Response", "Therapy")) |>
#'   annotate_public()
#'
#' idata_strata <- get_test_immundata() |>
#'   agg_repertoires(c("Response", "Therapy")) |>
#'   agg_strata(schema = "Response") |>
#'   annotate_public()
#'
#' write_immundata(idata_global, "snapshots/public-annotated-global")
#' write_immundata(idata_strata, "snapshots/public-annotated-strata")
#' }
#'
#' @concept Publicness
#' @export
annotate_public <- function(idata) {
  checkmate::assert_r6(idata, "ImmunData")

  if (is.null(idata$repertoires) || is.null(idata$schema_repertoire)) {
    cli::cli_abort(
      "Repertoire aggregation is required for {.fn annotate_public}. Run {.code agg_repertoires()} first."
    )
  }

  repertoire_col <- immundata::imd_schema("repertoire")
  receptor_col <- immundata::imd_schema("receptor")
  count_col <- immundata::imd_schema("count")
  prop_col <- immundata::imd_schema("proportion")
  rep_tbl <- idata$repertoires
  strata_col <- immundata::imd_schema("strata")

  required_cols <- c(repertoire_col, receptor_col, count_col, prop_col)
  missing_required_cols <- setdiff(required_cols, colnames(idata$annotations))
  if (length(missing_required_cols)) {
    cli::cli_abort(
      "Required annotation column(s) {missing_required_cols} are missing in {.field idata$annotations}."
    )
  }

  n_repertoires_total <- nrow(rep_tbl)
  if (n_repertoires_total == 0) {
    cli::cli_abort("No repertoires found in {.field idata$repertoires}.")
  }

  receptor_per_repertoire <- idata$annotations |>
    select(all_of(required_cols)) |>
    summarise(
      .public_count = max(!!rlang::sym(count_col), na.rm = TRUE),
      .public_prop = max(!!rlang::sym(prop_col), na.rm = TRUE),
      .by = all_of(c(receptor_col, repertoire_col))
    )

  public_ann <- receptor_per_repertoire |>
    summarise(
      imd_public_incidence = n(),
      imd_public_count_min = min(.data$.public_count, na.rm = TRUE),
      imd_public_count_max = max(.data$.public_count, na.rm = TRUE),
      imd_public_count_mean = mean(.data$.public_count, na.rm = TRUE),
      imd_public_count_median = stats::median(.data$.public_count, na.rm = TRUE),
      imd_public_prop_min = min(.data$.public_prop, na.rm = TRUE),
      imd_public_prop_max = max(.data$.public_prop, na.rm = TRUE),
      imd_public_prop_mean = mean(.data$.public_prop, na.rm = TRUE),
      imd_public_prop_median = stats::median(.data$.public_prop, na.rm = TRUE),
      .by = all_of(receptor_col)
    ) |>
    mutate(
      imd_public_incidence_prop =
        !!rlang::sym("imd_public_incidence") / n_repertoires_total
    )

  if (!(strata_col %in% colnames(rep_tbl))) {
    cli::cli_inform(
      "No strata found in {.field idata$repertoires}. Computing global publicness metrics only. Run {.code immundata::agg_strata()} to add per-strata metrics."
    )
  } else {
    rep_strata_map <- rep_tbl |>
      select(all_of(c(repertoire_col, strata_col))) |>
      distinct() |>
      filter(!is.na(!!rlang::sym(strata_col)))

    if (nrow(rep_strata_map) > 0) {
      rep_strata_map_db <- duckplyr::as_duckdb_tibble(rep_strata_map)
      strata_sizes <- rep_strata_map |>
        count(!!rlang::sym(strata_col), name = "imd_strata_n_repertoires")

      strata_ids <- strata_sizes |>
        pull(!!rlang::sym(strata_col))

      for (sid in strata_ids) {
        n_rep_in_strata <- strata_sizes |>
          filter(!!rlang::sym(strata_col) == sid) |>
          pull("imd_strata_n_repertoires")

        sid_chr <- as.character(sid)

        strata_stats <- receptor_per_repertoire |>
          inner_join(
            rep_strata_map_db |> filter(!!rlang::sym(strata_col) == sid),
            by = repertoire_col
          ) |>
          summarise(
            imd_public_incidence = n(),
            imd_public_count_min = min(.data$.public_count, na.rm = TRUE),
            imd_public_count_max = max(.data$.public_count, na.rm = TRUE),
            imd_public_count_mean = mean(.data$.public_count, na.rm = TRUE),
            imd_public_count_median = stats::median(.data$.public_count, na.rm = TRUE),
            imd_public_prop_min = min(.data$.public_prop, na.rm = TRUE),
            imd_public_prop_max = max(.data$.public_prop, na.rm = TRUE),
            imd_public_prop_mean = mean(.data$.public_prop, na.rm = TRUE),
            imd_public_prop_median = stats::median(.data$.public_prop, na.rm = TRUE),
            .by = all_of(receptor_col)
          ) |>
          mutate(
            imd_public_incidence_prop =
              !!rlang::sym("imd_public_incidence") / n_rep_in_strata
          )

        strata_metric_cols <- setdiff(colnames(strata_stats), receptor_col)
        strata_rename_map <- stats::setNames(
          strata_metric_cols,
          paste0(strata_metric_cols, "_strata_", sid_chr)
        )

        strata_stats <- strata_stats |>
          rename(!!!strata_rename_map)

        public_ann <- public_ann |>
          left_join(strata_stats, by = receptor_col)
      }
    }
  }

  immundata::annotate_receptors(
    idata = idata,
    annotations = public_ann,
    conflicts = "replace",
    remove_limit = TRUE
  )
}
