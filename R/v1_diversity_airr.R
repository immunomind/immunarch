#' @title Diversity - estimating the heterogeneity of immune repertoires
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **receptor diversity** per repertoire. A characteristic of a whole repertoire.
#'
#' ## Available functions
#'
#' Supported methods are the following.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams airr_diversity_dxx
#' @inheritParams airr_diversity_chao1
#' @inheritParams airr_diversity_shannon
#' @inheritParams airr_diversity_pielou
#' @inheritParams airr_diversity_hill
#' @inheritParams airr_diversity_index
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this session.
#' # Change this only if you know what you're doing (e.g., multi-user machines, shared CI/servers).
#' db_exec("SET threads TO 1")
#' # Load data
#' \dontrun{
#' immdata <- get_test_idata() |> agg_repertoires("Therapy")
#' }
#'
#' @name airr_diversity
#' @concept Diversity
NULL


#' @keywords internal
airr_diversity_dxx_impl <- function(idata, perc = 50) {
  checkmate::assert_numeric(perc, any.missing = FALSE)
  if (!all(perc > 0 & perc <= 100)) {
    cli::cli_abort("{.code perc} must be in (0, 100].")
  }

  rep_str <- immundata::imd_schema("repertoire")
  rep_sym <- immundata::imd_schema_sym("repertoire")
  rec_sym <- immundata::imd_schema_sym("receptor")
  prop_str <- immundata::imd_schema("proportion")
  prop_sym <- immundata::imd_schema_sym("proportion")

  base_tbl <- idata$annotations |>
    dplyr::select(!!rec_sym, !!rep_sym, !!prop_sym) |>
    dplyr::distinct(!!rec_sym, !!rep_sym, .keep_all = TRUE) |>
    dplyr::arrange()

  k_sql <- sprintf(
    "ROW_NUMBER() OVER (PARTITION BY %s ORDER BY %s DESC)",
    rep_str, prop_str
  )
  cum_sql <- sprintf(
    "SUM(%s) OVER (PARTITION BY %s ORDER BY %s DESC ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)",
    prop_str, rep_str, prop_str
  )

  ranked <- base_tbl |>
    duckplyr::as_tbl() |>
    dplyr::mutate(
      k   = dbplyr::sql(k_sql),
      cum = dbplyr::sql(cum_sql)
    ) |>
    duckplyr::as_duckdb_tibble()

  res <- purrr::map_dfr(perc, function(p) {
    ranked |>
      dplyr::filter(.data$cum >= p / 100) |>
      dplyr::group_by(!!rep_sym) |>
      dplyr::summarise(dxx = min(.data$k), .groups = "drop") |>
      dplyr::mutate(perc = p)
  }) |>
    dplyr::select(!!rep_sym, .data$perc, .data$dxx) |>
    dplyr::arrange(!!rep_sym, .data$perc) |>
    collect()

  res
}


#' @description `airr_diversity_dxx` - **coverage diversity**: minimal number of
#' top receptors needed to reach `perc%` of clonal space (by `proportion`).
#' Great for spotting dominance/overexpansion and for quick, interpretable dashboards
#' (e.g., D50 = receptors to cover half of the repertoire).
#'
#' @param perc A number or numeric vector in `(0, 100]` (default `50`), e.g.
#'   `50` for D50, `20` for D20.
#'
#' @return
#'
#' ## `airr_diversity_dxx`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `perc`
#' * `dxx` - minimal count of top receptors to reach `perc%`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_dxx
#' #
#' \dontrun{
#' d50 <- airr_diversity_dxx(immdata, perc = 50)
#' d_multi <- airr_diversity_dxx(immdata, perc = c(20, 50, 80))
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_dxx <- register_immunarch_method(airr_diversity_dxx_impl, "airr_diversity", "dxx")


#' @keywords internal
airr_diversity_chao1_impl <- function(idata) {
  rep_col <- immundata::imd_schema("repertoire")
  rep_sym <- immundata::imd_schema_sym("repertoire")
  cnt_sym <- immundata::imd_schema_sym("count")

  # TODO: optimize this please, loading all the data in R is not good.
  # TODO: check if no integer overflow
  idata$annotations |>
    select(!!rep_sym, !!cnt_sym) |>
    collect() |>
    summarise(counts = list(!!cnt_sym), .by = !!rep_sym) |>
    mutate(ch = lapply(counts, chao1)) |>
    transmute(
      !!rep_col := !!rep_sym,
      Estimator = vapply(ch, function(x) unname(x["Estimator.1"]), numeric(1)),
      SD = vapply(ch, function(x) unname(x["SD.2"]), numeric(1)),
      `Conf.95.lo` = vapply(ch, function(x) unname(x["Conf.95.lo.1"]), numeric(1)),
      `Conf.95.hi` = vapply(ch, function(x) unname(x["Conf.95.hi.1"]), numeric(1))
    ) |>
    collect()
}

#' @description `airr_diversity_chao1` - Chao1 estimator is a nonparameteric
#'  asymptotic estimator of species richness (number of species in a population).
#'  One of the most used methods for estimating immune repertoire diversity.
#'
#' @return
#'
#' ## `airr_diversity_chao1`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `Estimator` - number of species
#' * `SD` - standard deviation for the estimator value
#' * `Conf.95.lo` - CI 0.025
#' * `Conf.95.hi` - CI 0.975
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_chao1
#' #
#' \dontrun{
#' chao <- airr_diversity_chao1(immdata)
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_chao1 <- register_immunarch_method(airr_diversity_chao1_impl, "airr_diversity", "chao1")


#' @keywords internal
airr_diversity_shannon_impl <- function(idata) {
  idata$annotations |>
    select(
      !!immundata::imd_schema_sym("receptor"),
      !!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("proportion")
    ) |>
    distinct(!!immundata::imd_schema_sym("receptor"),
      !!immundata::imd_schema_sym("repertoire"),
      .keep_all = TRUE
    ) |>
    summarise(
      .by = !!immundata::imd_schema_sym("repertoire"),
      shannon = -sum(!!immundata::imd_schema_sym("proportion") * dd$log2(!!immundata::imd_schema_sym("proportion")))
    ) |>
    collect()
}


#' @description `airr_diversity_shannon` - Shannon entropy (base 2) per repertoire
#' computed from `proportion`. Ideal when you want a single evenness-aware
#' diversity score; pair with Pielou/Hill for samples with very different richness.
#'
#' @return
#'
#' ## `airr_diversity_shannon`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `shannon` - entropy in bits
#'
#' @examples
#' #
#' # airr_diversity_shannon
#' #
#' \dontrun{
#' sh <- airr_diversity_shannon(immdata)
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_shannon <- register_immunarch_method(airr_diversity_shannon_impl, "airr_diversity", "shannon")


#' @keywords internal
airr_diversity_pielou_impl <- function(idata) {
  shannon_values <- airr_diversity_shannon(idata, autojoin = FALSE)

  idata$repertoires |>
    select(c(immundata::imd_schema("repertoire"), immundata::imd_schema("n_receptors"))) |>
    left_join(shannon_values,
      by = imd_schema("repertoire")
    ) |>
    mutate(richness = dd$log2(!!immundata::imd_schema_sym("n_receptors")), pielou = shannon / richness) |>
    collect()
}


#' @description `airr_diversity_pielou` - Pielou's evenness `H / log2(S)` with
#' richness `S`. Best when you need a **size-normalized** evenness score that's
#' comparable across repertoires with different receptor counts.
#'
#' @return
#'
#' ## `airr_diversity_pielou`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `shannon`
#' * `n_receptors`
#' * `pielou` - evenness in `[0, 1]` (NA if `S <= 1`)
#'
#' @examples
#' #
#' # airr_diversity_pielou
#' #
#' \dontrun{
#' pj <- airr_diversity_pielou(immdata)
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_pielou <- register_immunarch_method(airr_diversity_pielou_impl, "airr_diversity", "pielou")


#' @keywords internal
airr_diversity_index_impl <- function(idata) {
  airr_diversity_hill(idata, q = 1)
}


#' @description `airr_diversity_index` - convenience alias for Hill number with
#' `q = 1` (`exp(Shannon)` using natural log). A solid **default single metric**
#' that's relatively robust to rare-count noise and easy to compare across samples.
#'
#' @return
#'
#' ## `airr_diversity_index`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `q = 1`
#' * `hill_number`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_index
#' #
#' \dontrun{
#' idx <- airr_diversity_index(immdata)
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_index <- register_immunarch_method(airr_diversity_index_impl, "airr_diversity", "index")


#' @keywords internal
airr_diversity_hill_impl <- function(idata, q = 0:5) {
  checkmate::check_numeric(q, lower = 0, sorted = TRUE)

  receptors <- idata$annotations |>
    select(
      !!immundata::imd_schema_sym("receptor"),
      !!immundata::imd_schema_sym("repertoire"),
      !!immundata::imd_schema_sym("proportion")
    ) |>
    distinct(!!immundata::imd_schema_sym("receptor"),
      !!immundata::imd_schema_sym("repertoire"),
      .keep_all = TRUE
    )

  result <- NULL

  # TODO: Join by value of q-s and run? if_else in case of different q-s
  for (q_val in q) {
    if (q_val == 0) {
      q_val_tbl <- idata$repertoires |>
        summarise(
          .by = !!immundata::imd_schema_sym("repertoire"),
          q = 0,
          hill_number = as.numeric(!!immundata::imd_schema_sym("n_receptors"))
        )
    } else if (q_val == 1) {
      q_val_tbl <- receptors |>
        summarise(
          .by = !!immundata::imd_schema_sym("repertoire"),
          q = q_val,
          hill_number = dd$exp(-sum(!!immundata::imd_schema_sym("proportion") * dd$ln(!!immundata::imd_schema_sym("proportion"))))
        )
    } else {
      q_val_tbl <- receptors |>
        summarise(
          .by = !!immundata::imd_schema_sym("repertoire"),
          q = q_val,
          hill_number = dd$pow(sum(dd$pow(!!immundata::imd_schema_sym("proportion"), q_val)), 1 / (1 - q_val))
        )
    }

    if (is.null(result)) {
      result <- q_val_tbl
    } else {
      result <- result |> union_all(q_val_tbl)
    }
  }

  idata$metadata |>
    left_join(result, by = imd_schema("repertoire")) |>
    collect()
}


#' @description `airr_diversity_hill` - Hill numbers ("true diversity") for
#' orders `q \eqn{\in}{in} {0, 1, 2, ...}`: `q=0` richness, `q=1` exp(Shannon), `q>1`
#' emphasizes abundant receptors. Perfect when you want a **diversity profile**
#' that tunes sensitivity to rare vs. abundant clonotypes.
#'
#' @inheritParams im_common_args
#' @param q A scalar or vector of non-negative orders. Defaults to `0:5`.
#'
#' @return
#'
#' ## `airr_diversity_hill`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `q` - Hill order
#' * `hill_number` - true diversity of order `q`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_hill
#' #
#' \dontrun{
#' hill <- airr_diversity_hill(immdata, q = c(0, 1, 2))
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_hill <- register_immunarch_method(airr_diversity_hill_impl, "airr_diversity", "hill")
