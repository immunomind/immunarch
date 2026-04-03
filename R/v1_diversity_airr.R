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
#' @inheritParams airr_diversity_rarefaction
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
  rec_sym <- immundata::imd_schema_sym("receptor")
  cnt_sym <- immundata::imd_schema_sym("count")

  # TODO: optimize this please, loading all the data in R is not good.
  # TODO: check if no integer overflow
  idata$annotations |>
    select(!!rep_sym, !!rec_sym, !!cnt_sym) |>
    distinct(!!rec_sym, !!rep_sym, .keep_all = TRUE) |>
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
airr_diversity_rarefaction_collect_rep_counts <- function(idata) {
  rep_sym <- immundata::imd_schema_sym("repertoire")
  rec_sym <- immundata::imd_schema_sym("receptor")
  cnt_sym <- immundata::imd_schema_sym("count")

  idata$annotations |>
    dplyr::select(!!rep_sym, !!rec_sym, !!cnt_sym) |>
    dplyr::distinct(!!rec_sym, !!rep_sym, .keep_all = TRUE) |>
    dplyr::select(!!rep_sym, !!cnt_sym) |>
    dplyr::collect()
}


#' @keywords internal
airr_diversity_rarefaction_chao1_stats <- function(counts_vec) {
  counts_vec <- as.numeric(counts_vec)
  counts_vec <- counts_vec[is.finite(counts_vec) & counts_vec >= 0]

  if (length(counts_vec) == 0 || sum(counts_vec) <= 0) {
    return(c(Estimator = 0, SD = 0, `Conf.95.lo` = 0, `Conf.95.hi` = 0))
  }

  counts <- table(counts_vec)
  n <- sum(counts_vec)
  D <- length(counts_vec)
  f1 <- counts["1"]
  f2 <- counts["2"]

  if (is.na(f1) && is.na(f2)) {
    e <- D
    i <- unique(counts_vec)
    v <- sum(sapply(i, function(j) sum(counts_vec == j) * (exp(-j) - exp(-2 * j)))) -
      (sum(sapply(i, function(j) j * exp(-j) * sum(counts_vec == j))))^2 / n
    P <- sum(sapply(i, function(j) sum(counts_vec == j) * exp(-j) / D))
    lo <- max(D, D / (1 - P) - stats::qnorm(1 - .05 / 2) * sqrt(v) / (1 - P))
    hi <- D / (1 - P) + stats::qnorm(1 - .05 / 2) * sqrt(v) / (1 - P)
  } else if (is.na(f2)) {
    e <- D + f1 * (f1 - 1) / 2 * (n - 1) / n
    v <- (n - 1) / n * f1 * (f1 - 1) / 2 +
      ((n - 1) / n)^2 * f1 * (2 * f1 - 1)^2 / 4 -
      ((n - 1) / n)^2 * f1^4 / 4 / e
    t_val <- e - D
    K <- exp(stats::qnorm(1 - .05 / 2) * sqrt(log(1 + v / t_val^2)))
    lo <- D + t_val / K
    hi <- D + t_val * K
  } else {
    const <- (n - 1) / n
    e <- D + f1^2 / (2 * f2) * const
    f12 <- f1 / f2
    v <- f2 * (const * f12^2 / 2 + const^2 * f12^3 + const^2 * f12^4 / 4)
    t_val <- e - D
    K <- exp(stats::qnorm(1 - .05 / 2) * sqrt(log(1 + v / t_val^2)))
    lo <- D + t_val / K
    hi <- D + t_val * K
  }

  c(
    Estimator = as.numeric(e),
    SD = as.numeric(sqrt(v)),
    `Conf.95.lo` = as.numeric(lo),
    `Conf.95.hi` = as.numeric(hi)
  )
}


#' @keywords internal
airr_diversity_rarefaction_impl <- function(idata, step = NA, quantile = c(.025, .975),
                                            extrapolation = NA, norm = TRUE, verbose = TRUE) {
  checkmate::assert_logical(norm, len = 1)
  checkmate::assert_logical(verbose, len = 1)
  checkmate::assert_numeric(quantile, len = 2, any.missing = FALSE, lower = 0, upper = 1)

  if (!is.na(step)) {
    checkmate::assert_number(step, lower = 1, finite = TRUE)
  }

  if (!is.na(extrapolation)) {
    checkmate::assert_number(extrapolation, lower = 0, finite = TRUE)
  }

  quantile <- sort(quantile)
  if (quantile[1] == quantile[2]) {
    cli::cli_abort("{.code quantile} values must be different.")
  }

  rep_col <- immundata::imd_schema("repertoire")
  cnt_col <- immundata::imd_schema("count")
  counts_tbl <- airr_diversity_rarefaction_collect_rep_counts(idata)

  if (!nrow(counts_tbl)) {
    out <- tibble::tibble(
      size = numeric(0),
      q_low = numeric(0),
      mean = numeric(0),
      q_high = numeric(0),
      type = character(0)
    )
    out[[rep_col]] <- numeric(0)
    out <- out |>
      dplyr::select(all_of(c(rep_col, "size", "q_low", "mean", "q_high", "type")))
    return(out)
  }

  counts_split <- counts_tbl |>
    dplyr::summarise(counts = list(.data[[cnt_col]]), .by = all_of(rep_col))

  total_counts <- vapply(counts_split$counts, function(x) sum(as.numeric(x), na.rm = TRUE), numeric(1))

  if (is.na(step)) {
    min_total <- min(total_counts[total_counts > 0], na.rm = TRUE)
    step <- if (is.finite(min_total)) floor(min_total / 50) else 1
  }
  step <- max(1L, as.integer(step))

  if (is.na(extrapolation)) {
    extrapolation <- max(total_counts, na.rm = TRUE) * 20
  }
  extrapolation <- as.numeric(extrapolation)

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Computing rarefaction in RAM for {nrow(counts_split)} repertoire(s).")
  }

  res_list <- lapply(seq_len(nrow(counts_split)), function(i) {
    rep_id <- counts_split[[rep_col]][i]
    bc_vec <- as.numeric(counts_split$counts[[i]])
    bc_vec <- bc_vec[is.finite(bc_vec) & bc_vec > 0]

    if (!length(bc_vec)) {
      return(NULL)
    }

    Sobs <- length(bc_vec)
    n <- sum(bc_vec)

    if (n <= 0) {
      return(NULL)
    }

    ch_stats <- airr_diversity_rarefaction_chao1_stats(bc_vec)
    Sest <- unname(ch_stats["Estimator"])
    if (!is.finite(Sest)) {
      Sest <- Sobs
    }

    sizes <- seq(step, n, step)
    if (!length(sizes) || tail(sizes, 1) != n) {
      sizes <- c(sizes, n)
    }
    sizes <- sort(unique(sizes))

    count_freq <- table(bc_vec)
    freqs <- as.numeric(names(count_freq))
    freq_mult <- as.numeric(count_freq)
    z_value <- stats::qnorm(quantile[2])

    interpolation_df <- do.call(rbind, lapply(sizes, function(sz) {
      alpha <- (1 - sz / n)^freqs
      Sind <- sum((1 - alpha) * freq_mult)

      if (Sest == Sobs) {
        SD <- 0
      } else {
        var_est <- sum((1 - alpha)^2 * freq_mult) - Sind^2 / Sest
        SD <- sqrt(max(var_est, 0))
      }

      t_val <- Sind - Sobs
      if (t_val != 0) {
        K <- exp(z_value * sqrt(log(1 + (SD / t_val)^2)))
        ci_1 <- Sobs + t_val / K
        ci_2 <- Sobs + t_val * K
        low <- min(ci_1, ci_2)
        high <- max(ci_1, ci_2)
      } else {
        low <- Sind
        high <- Sind
      }

      tibble::tibble(
        size = as.numeric(sz),
        q_low = as.numeric(low),
        mean = as.numeric(Sind),
        q_high = as.numeric(high),
        type = "interpolation"
      )
    }))

    extrapolation_df <- NULL
    if (extrapolation > 0) {
      extrap_sizes <- seq(max(sizes) + step, extrapolation, step)
      if (length(extrap_sizes)) {
        f0 <- Sest - Sobs
        f1 <- unname(count_freq["1"])
        extrapolation_df <- do.call(rbind, lapply(extrap_sizes, function(sz) {
          if (is.na(f1) || f0 == 0) {
            Sind <- Sobs
          } else {
            Sind <- Sobs + f0 * (1 - exp(-(sz - n) / n * f1 / f0))
          }
          tibble::tibble(
            size = as.numeric(sz),
            q_low = as.numeric(Sind),
            mean = as.numeric(Sind),
            q_high = as.numeric(Sind),
            type = "extrapolation"
          )
        }))
      }
    }

    out <- dplyr::bind_rows(interpolation_df, extrapolation_df)
    if (isTRUE(norm)) {
      out <- out |>
        dplyr::mutate(
          size = .data$size / n,
          q_low = .data$q_low / Sobs,
          mean = .data$mean / Sobs,
          q_high = .data$q_high / Sobs
        )
    }

    out[[rep_col]] <- rep_id
    out |>
      dplyr::select(all_of(c(rep_col, "size", "q_low", "mean", "q_high", "type")))
  })

  dplyr::bind_rows(res_list) |>
    dplyr::arrange(.data[[rep_col]], .data$size, .data$type)
}


#' @description `airr_diversity_rarefaction` - interpolation/extrapolation curves
#' for receptor richness as a function of sampled clones. Computation is done in
#' RAM after one DB query that fetches only repertoire id and receptor counts.
#'
#' @param step Rarefaction step size. Defaults to `floor(min_total_clones / 50)`,
#'   lower-bounded by `1`.
#' @param quantile Numeric vector of length 2 with confidence interval bounds.
#' @param extrapolation Upper size limit for extrapolation. Use `0` to disable
#'   extrapolation. Defaults to `max_total_clones * 20`.
#' @param norm Logical; if `TRUE`, size is divided by total repertoire size and
#'   richness estimates are divided by observed richness.
#' @param verbose Logical; show a concise progress message.
#'
#' @return
#'
#' ## `airr_diversity_rarefaction`
#' A tibble with:
#' * `imd_repertoire_id`
#' * `size` - sample size (absolute or normalized)
#' * `q_low` - lower confidence bound
#' * `mean` - expected richness
#' * `q_high` - upper confidence bound
#' * `type` - `interpolation` or `extrapolation`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_rarefaction
#' #
#' \dontrun{
#' raref <- airr_diversity_rarefaction(immdata, step = 2000, extrapolation = 0)
#' vis(raref)
#' }
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_rarefaction <- register_immunarch_method(airr_diversity_rarefaction_impl, "airr_diversity", "rarefaction")


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
