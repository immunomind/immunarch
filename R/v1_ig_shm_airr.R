#' Summarize immunoglobulin somatic hypermutation per repertoire
#'
#' Summarize the mutation counts and frequencies produced by
#' [annotate_ig_shm()]. Aggregation runs in DuckDB on the duckplyr-backed input.
#'
#' @param idata An [immundata::ImmunData] object with SHM annotations and clones.
#' @param by Character vector, `NA`, or `NULL`. The default uses the canonical
#'   locus column when available. Use `NULL` to pool loci. Repertoires are always
#'   summarized separately; clone IDs cannot be included in `by`.
#' @param clone_col Column containing clone assignments. Missing clone IDs are
#'   included in sequence summaries but excluded from clone summaries.
#' @inheritParams method_args_airr
#' @param format Inherited from the common AIRR interface. This method always
#'   returns long-format summaries.
#'
#' @details
#' This function uses the mutation counts and frequencies already calculated
#' by [annotate_ig_shm()]. Here, a sequence means one immunoglobulin chain
#' sequence in the input.
#'
#' The function summarizes these counts and frequencies in two ways:
#'
#' * Across sequences: each sequence has equal weight. A clone with more
#'   sequences therefore contributes more to the result.
#' * Across clones: first calculate a mean or median within each clone, then
#'   summarize these clone means or medians. Each clone has equal weight,
#'   regardless of how many sequences it contains.
#'
#' These calculations do not use read counts or UMI counts as weights.
#'
#' @section Output statistics:
#' Each repertoire has several output rows. The columns `region`, `mutation`,
#' and `measure` tell you which mutation counts or frequencies the row describes:
#'
#' | Column | Meaning |
#' |:-------|:--------|
#' | `region` | The region being analyzed: `overall` for all available annotated regions; `v` for FWR1 through FWR3; `cdr` for all available CDRs; `fwr` for all available framework regions; or one region such as `cdr3`. |
#' | `mutation` | `replacement` for mutations that change the amino acid, `silent` for mutations that do not change the amino acid, or `total` for replacement plus silent mutations. See [annotate_ig_shm()] for counting rules. |
#' | `measure` | `count` for the number of mutations, or `freq` for mutation frequency. Frequencies are fractions: `0.05` means 5%. |
#'
#' The other columns contain the statistics for that row. For example,
#' consider four sequences in one repertoire. Clone A has total CDR3 mutation
#' counts of 0, 0, and 9. Clone B has one sequence with 6 total CDR3 mutations.
#' Part of the output would be:
#'
#' | imd_repertoire_id | region | mutation | measure | sequence_mean | clone_mean_mean | clone_median_mean |
#' |:------------------|:-------|:---------|:--------|--------------:|----------------:|------------------:|
#' | 1 | cdr3 | total | count | 3.75 | 4.5 | 3 |
#'
#' All three statistics in this row describe total CDR3 mutation counts:
#'
#' * `sequence_mean`: average the four counts: `(0 + 0 + 9 + 6) / 4 = 3.75`.
#' * `clone_mean_mean`: the clone means are 3 and 6. Their mean is 4.5.
#' * `clone_median_mean`: the clone medians are 0 and 6. Their mean is 3.
#'
#' A row with `measure = "freq"` contains the same kinds of statistics for
#' mutation frequencies instead of mutation counts. Results are calculated
#' separately for each repertoire and each group selected with `by`, such as
#' locus. Repertoire metadata is added when `autojoin = TRUE`.
#'
#' The full list of 24 statistics is below. In every description, the mutation
#' count or frequency is the one specified by `region`, `mutation`, and
#' `measure` in that row. `q25` and `q75` are the lower and upper quartiles
#' (25th and 75th percentiles).
#'
#' | Column | Meaning |
#' |:-------|:--------|
#' | `n_sequences` | Number of input sequences, including those with a missing mutation count or frequency, or a missing clone ID. |
#' | `n_valid` | Number of sequences with an available mutation count or frequency. |
#' | `n_missing` | Number of sequences whose mutation count or frequency is `NA`: `n_sequences - n_valid`. |
#' | `n_valid_without_clone` | Number of sequences with an available mutation count or frequency but no clone ID. |
#' | `n_zero` | Number of sequences with a mutation count or frequency of zero. |
#' | `zero_prop` | Proportion of sequences with zero mutations of this type in this region, among sequences with an available count or frequency: `n_zero / n_valid`. |
#' | `sequence_mean` | Mean mutation count or frequency across sequences. Each sequence has equal weight. |
#' | `sequence_median` | Median mutation count or frequency across sequences. |
#' | `sequence_q25` | Lower quartile of mutation counts or frequencies across sequences. |
#' | `sequence_q75` | Upper quartile of mutation counts or frequencies across sequences. |
#' | `n_clones` | Number of clones with a known clone ID, including clones for which all mutation counts or frequencies are missing. |
#' | `n_valid_clones` | Number of clones with at least one sequence whose mutation count or frequency is available. |
#' | `n_valid_with_clone` | Number of sequences with both an available mutation count or frequency and a clone ID. |
#' | `clone_mean_weighted_mean` | Mean mutation count or frequency across sequences with a clone ID. This is also the mean of clone means weighted by the number of sequences with available counts or frequencies in each clone. |
#' | `clone_mean_mean` | Calculate the mean mutation count or frequency in each clone, then take the mean of these clone means. Each clone has equal weight. |
#' | `clone_mean_median` | Calculate the mean mutation count or frequency in each clone, then take the median of these clone means. |
#' | `clone_mean_q25` | Lower quartile of the clone means. |
#' | `clone_mean_q75` | Upper quartile of the clone means. |
#' | `clone_median_mean` | Calculate the median mutation count or frequency in each clone, then take the mean of these clone medians. Each clone has equal weight. |
#' | `clone_median_median` | Calculate the median mutation count or frequency in each clone, then take the median of these clone medians. |
#' | `clone_median_q25` | Lower quartile of the clone medians. |
#' | `clone_median_q75` | Upper quartile of the clone medians. |
#' | `n_zero_clone_means` | Number of clones whose mean mutation count or frequency is zero. |
#' | `n_zero_clone_medians` | Number of clones whose median mutation count or frequency is zero. |
#'
#' Missing mutation counts or frequencies (`NA`) are excluded from means,
#' medians, quartiles, and the denominator of `zero_prop`. Zero is included.
#' Sequences without clone IDs are included in sequence statistics but excluded
#' from clone statistics. Clones with no available mutation counts or frequencies
#' are included in `n_clones` but excluded from the clone means, medians, and
#' quartiles. If none of the required counts or frequencies are available, the
#' mean, median, quartiles, and `zero_prop` are `NA`. Quartiles use linear
#' interpolation.
#'
#' Combined region frequencies come from [annotate_ig_shm()], which divides
#' the total mutation count by the total number of comparable bases in those
#' regions. This function does not average the frequencies of different regions
#' to obtain a combined frequency. Regions without the required input columns
#' are skipped.
#'
#' @return A long-format table of repertoire statistics. Each row contains the
#'   repertoire ID, any `by` columns, `region`, `mutation`, `measure`, and the
#'   24 statistics listed above. Repertoire metadata is added by default.
#' @name airr_ig_shm
NULL


#' @keywords internal
airr_ig_shm_impl <- function(idata, by = NA_character_, clone_col = "clone_id") {
  checkmate::assert_string(clone_col)
  by <- resolve_argument_by(idata, by)
  repertoire <- immundata::imd_schema("repertoire")
  groups <- unique(c(repertoire, by))
  columns <- colnames(idata)
  required <- c(groups, clone_col, "mu_count", "mu_freq",
                "mu_count_r", "mu_count_s", "mu_freq_r", "mu_freq_s")
  if (length(setdiff(required, columns))) {
    cli::cli_abort("Missing SHM column(s): {.field {setdiff(required, columns)}}. Run {.code annotate_ig_shm()} and supply clone assignments first.")
  }
  if (clone_col %in% groups) cli::cli_abort("{.arg by} must not include {.arg clone_col}; this function returns repertoire summaries.")
  source <- idata$annotations |> duckplyr::as_tbl()
  con <- dbplyr::remote_con(source)
  qi <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  qs <- function(x) as.character(DBI::dbQuoteString(con, x))
  # Each entry selects an existing measurement; regional totals are R + S.
  regions <- c(overall = "", v = "_v", cdr = "_cdr", fwr = "_fwr",
               fwr1 = "_fwr1", cdr1 = "_cdr1", fwr2 = "_fwr2", cdr2 = "_cdr2",
               fwr3 = "_fwr3", cdr3 = "_cdr3", fwr4 = "_fwr4")
  entries <- character()
  for (region in names(regions)) for (measure in c("count", "freq")) {
    base <- paste0("mu_", measure, regions[[region]])
    if (!all(paste0(base, c("_r", "_s")) %in% columns)) next
    for (mutation in c("total", "replacement", "silent")) {
      column <- paste0(base, switch(mutation, total = "", replacement = "_r", silent = "_s"))
      value <- if (column %in% columns) qi(column) else {
        paste0("(", qi(paste0(base, "_r")), " + ", qi(paste0(base, "_s")), ")")
      }
      entries <- c(entries, paste0("(", qs(region), ", ", qs(mutation), ", ",
                                  qs(measure), ", CAST(", value, " AS DOUBLE))"))
    }
  }
  keys <- c(groups, "region", "mutation", "measure")
  keys_sql <- paste(qi(keys), collapse = ", ")
  # All aggregation below is executed by DuckDB on the duckplyr-backed input.
  statistics <- function(value, prefix) {
    paste0(c("avg(", "median(", "quantile_cont(", "quantile_cont("), value,
           c(")", ")", ", 0.25)", ", 0.75)"), " AS ", prefix,
           c("mean", "median", "q25", "q75"), collapse = ", ")
  }
  join <- paste0("s.", qi(keys), " IS NOT DISTINCT FROM c.", qi(keys), collapse = " AND ")
  sql <- paste0(
    "WITH measurements AS (SELECT ", paste(qi(groups), collapse = ", "), ", ",
    qi(clone_col), " AS clone_key, m.* FROM (", dbplyr::sql_render(source), ") input, ",
    "LATERAL (VALUES ", paste(entries, collapse = ", "), ") m(region, mutation, measure, value)), ",
    "sequences AS (SELECT ", keys_sql, ", count(*) AS n_sequences, count(value) AS n_valid, ",
    "count(*) - count(value) AS n_missing, ",
    "count(value) FILTER (WHERE clone_key IS NULL) AS n_valid_without_clone, ",
    "count(*) FILTER (WHERE value = 0) AS n_zero, ",
    "count(*) FILTER (WHERE value = 0)::DOUBLE / nullif(count(value), 0) AS zero_prop, ",
    statistics("value", "sequence_"), " FROM measurements GROUP BY ", keys_sql, "), ",
    "clones AS (SELECT ", keys_sql, ", clone_key, count(value) AS clone_n, ",
    "avg(value) AS clone_mean, median(value) AS clone_median FROM measurements ",
    "WHERE clone_key IS NOT NULL GROUP BY ", keys_sql, ", clone_key), ",
    "clone_summary AS (SELECT ", keys_sql, ", count(*) AS n_clones, ",
    "count(clone_mean) AS n_valid_clones, sum(clone_n) AS n_valid_with_clone, ",
    "sum(clone_mean * clone_n) / nullif(sum(clone_n), 0) AS clone_mean_weighted_mean, ",
    statistics("clone_mean", "clone_mean_"), ", ", statistics("clone_median", "clone_median_"), ", ",
    "count(*) FILTER (WHERE clone_mean = 0) AS n_zero_clone_means, ",
    "count(*) FILTER (WHERE clone_median = 0) AS n_zero_clone_medians ",
    "FROM clones GROUP BY ", keys_sql, ") ",
    "SELECT s.*, c.* EXCLUDE (", keys_sql, ") FROM sequences s ",
    "LEFT JOIN clone_summary c ON ", join
  )
  result <- dplyr::tbl(con, dbplyr::sql(sql)) |>
    duckplyr::as_duckdb_tibble() |>
    dplyr::collect()
  count_columns <- c("n_clones", "n_valid_clones", "n_valid_with_clone",
                     "n_zero_clone_means", "n_zero_clone_medians")
  result[count_columns] <- lapply(result[count_columns], function(x) replace(x, is.na(x), 0))
  result
}


#' @rdname airr_ig_shm
#' @export
airr_ig_shm <- register_method_airr(
  core = airr_ig_shm_impl,
  family = "airr_ig",
  name = "shm",
  required_cols = "clone_col"
)
