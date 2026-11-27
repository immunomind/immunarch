#' @title Diversity - estimating the heterogeneity of immune repertoires
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **receptor diversity** in each repertoire.
#' Diversity describes how many different receptors are present and how evenly
#' their abundances are distributed. These methods help you compare repertoire
#' structure, detect clonal expansion, and assess sampling depth.
#'
#' ## Available functions
#'
#' The following methods are available.
#'
#' @param idata An [immundata::ImmunData] object.
#' @inheritParams airr_diversity_dxx
#' @inheritParams airr_diversity_chao1
#' @inheritParams airr_diversity_shannon
#' @inheritParams airr_diversity_pielou
#' @inheritParams airr_diversity_hill
#' @inheritParams airr_diversity_index
#' @inheritParams airr_diversity_rarefaction
#' @inheritParams method_args_airr
#'
#' @seealso [immundata::ImmunData]
#'
#' @section Visualisation:
#' Most `airr_diversity_*()` results can be passed directly to [vis()]. With
#' `autojoin = TRUE` (the default), repertoire metadata is included in the
#' result and can be selected by the plotting arguments.
#'
#' ## 1) Coverage diversity (`airr_diversity_dxx`)
#'
#' `vis()` plots `dxx` for each repertoire. When several `perc` values are
#' present, use `fill = "perc"` to show a separate column for each coverage
#' threshold.
#'
#' ## 2) Chao1 richness (`airr_diversity_chao1`)
#'
#' `vis()` produces a column plot of the estimated richness (`Estimator`) for
#' each repertoire.
#'
#' ## 3) Rarefaction and extrapolation (`airr_diversity_rarefaction`)
#'
#' `vis()` plots estimated richness (`mean`) against sample size (`size`). Use
#' `color` to group curves, `show_ci = FALSE` to hide confidence intervals, and
#' `log = TRUE` to use a logarithmic x-axis.
#'
#' ## 4) Shannon entropy (`airr_diversity_shannon`)
#'
#' `vis()` produces a column plot of `shannon` for each repertoire.
#'
#' ## 5) Pielou evenness (`airr_diversity_pielou`)
#'
#' `vis()` produces a column plot of `pielou` for each repertoire.
#'
#' ## 6) Diversity index (`airr_diversity_index`)
#'
#' `vis()` produces a column plot of `hill_number` at `q = 1` for each
#' repertoire.
#'
#' ## 7) Hill numbers (`airr_diversity_hill`)
#'
#' This result does not currently have a dedicated `vis()` method. Use the
#' returned `q` and `hill_number` columns to create a diversity profile.
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this example.
#' # Generally, you should NOT do this in your session.
#' db_exec("SET threads TO 1")
#'
#' # Load example data.
#' immdata <- get_test_idata()
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
    dplyr::arrange() |>
    duckplyr::as_tbl()
  con <- dbplyr::remote_con(base_tbl)

  k_sql <- glue::glue_sql(
    "ROW_NUMBER() OVER (
      PARTITION BY {`rep_str`} ORDER BY {`prop_str`} DESC
    )",
    .con = con
  )
  cum_sql <- glue::glue_sql(
    "SUM({`prop_str`}) OVER (
      PARTITION BY {`rep_str`} ORDER BY {`prop_str`} DESC
      ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
    )",
    .con = con
  )
  n_sql <- glue::glue_sql("COUNT(*) OVER (PARTITION BY {`rep_str`})", .con = con)

  ranked <- base_tbl |>
    dplyr::mutate(
      k   = dbplyr::sql(k_sql),
      cum = dbplyr::sql(cum_sql),
      n = dbplyr::sql(n_sql)
    ) |>
    duckplyr::as_duckdb_tibble()

  thresholds <- tibble::tibble(
    perc = perc,
    perc_id = seq_along(perc)
  ) |>
    duckplyr::as_duckdb_tibble()

  res <- ranked |>
    dplyr::mutate(.join_key = 1L) |>
    dplyr::inner_join(
      thresholds |> dplyr::mutate(.join_key = 1L),
      by = ".join_key"
    ) |>
    dplyr::filter(
      (perc == 100 & k == n) |
        (perc < 100 & cum >= perc / 100 - 8 * .Machine$double.eps)
    ) |>
    dplyr::group_by(!!rep_sym, .data$perc, .data$perc_id) |>
    dplyr::summarise(dxx = min(.data$k), .groups = "drop") |>
    dplyr::select(dplyr::all_of(c(rep_str, "perc", "dxx"))) |>
    dplyr::arrange(!!rep_sym, .data$perc) |>
    collect()

  res
}


#' @description
#' **1) Coverage diversity (`airr_diversity_dxx`).** Calculate the minimum
#' number of the most abundant receptors needed to cover `perc%` of a
#' repertoire. For example, D50 is the number of receptors that cover half of
#' the repertoire. Use it to identify repertoires dominated by expanded
#' receptors.
#'
#' @param perc A number or numeric vector in `(0, 100]` (default `50`), e.g.
#'   `50` for D50, `20` for D20.
#'
#' @return
#'
#' ## 1) Coverage diversity (`airr_diversity_dxx`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `perc`
#' * `dxx` - minimal count of top receptors to reach `perc%`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # Calculate D50 and several coverage thresholds.
#' d50 <- airr_diversity_dxx(immdata, perc = 50)
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' d_multi <- airr_diversity_dxx(immdata, perc = c(20, 50, 80))
#'
#' # Visualise the coverage diversity.
#' vis(d_multi)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_dxx <- register_method_airr(airr_diversity_dxx_impl, "airr_diversity", "dxx")


chao1 <- function(.data) {
  counts <- table(.data)
  e <- NA
  v <- NA
  lo <- NA
  hi <- NA
  n <- sum(.data)
  D <- length(.data)
  f1 <- counts["1"]
  f2 <- counts["2"]
  # No singletons
  if (is.na(f1)) {
    e <- D
    i <- unique(.data)
    v <- sum(sapply(i, function(i) sum(.data == i) * (exp(-i) - exp(-2 * i)))) - (sum(sapply(i, function(i) i * exp(-i) * sum(.data == i))))^2 / n
    P <- sum(sapply(i, function(i) sum(.data == i) * exp(-i) / D))
    lo <- max(D, D / (1 - P) - stats::qnorm(1 - .05 / 2) * sqrt(v) / (1 - P))
    hi <- D / (1 - P) + stats::qnorm(1 - .05 / 2) * sqrt(v) / (1 - P)
  }
  # f1 != 0 && f2 == 0
  else if (is.na(f2)) {
    e <- D + f1 * (f1 - 1) / 2 * (n - 1) / n
    v <- (n - 1) / n * f1 * (f1 - 1) / 2 + ((n - 1) / n)^2 * f1 * (2 * f1 - 1)^2 / 4 - ((n - 1) / n)^2 * f1^4 / 4 / e
    t <- e - D
    K <- exp(stats::qnorm(1 - .05 / 2) * sqrt(log(1 + v / t^2)))
    lo <- D + t / K
    hi <- D + t * K
  }
  # f1 != && f2 != 0
  else {
    const <- (n - 1) / n
    e <- D + f1^2 / (2 * f2) * const
    f12 <- f1 / f2
    v <- f2 * (const * f12^2 / 2 + const^2 * f12^3 + const^2 * f12^4 / 4)
    t <- e - D
    K <- exp(stats::qnorm(.975) * sqrt(log(1 + v / t^2)))
    lo <- D + t / K
    hi <- D + t * K
  }
  c(Estimator = e, SD = sqrt(v), Conf.95.lo = lo, Conf.95.hi = hi)
}


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
      Estimator = vapply(ch, function(x) unname(x[[1]]), numeric(1)),
      SD = vapply(ch, function(x) unname(x[[2]]), numeric(1)),
      `Conf.95.lo` = vapply(ch, function(x) unname(x[[3]]), numeric(1)),
      `Conf.95.hi` = vapply(ch, function(x) unname(x[[4]]), numeric(1))
    ) |>
    collect()
}

#' @description
#' **2) Chao1 richness (`airr_diversity_chao1`).** Estimate the total number of
#' receptors, including receptors that may be missing because of limited
#' sampling. Use this non-parametric estimator when rare receptors suggest that
#' the observed repertoire is incomplete.
#'
#' @return
#'
#' ## 2) Chao1 richness (`airr_diversity_chao1`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `Estimator` - estimated receptor richness
#' * `SD` - standard deviation for the estimator value
#' * `Conf.95.lo` - CI 0.025
#' * `Conf.95.hi` - CI 0.975
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' #
#' # Estimate receptor richness.
#' chao <- airr_diversity_chao1(immdata)
#'
#' # Visualise the richness estimates.
#' vis(chao)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_chao1 <- register_method_airr(airr_diversity_chao1_impl, "airr_diversity", "chao1")


#' @keywords internal
airr_diversity_rarefaction_impl <- function(idata, step = NA, quantile = c(.025, .975),
                                            extrapolation = NA, nboot = 50,
                                            norm = TRUE, verbose = TRUE) {
  estimate_unseen_richness <- function(counts) {
    n <- sum(counts)
    singleton_count <- sum(counts == 1)
    doubleton_count <- sum(counts == 2)

    if (doubleton_count > 0) {
      (n - 1) / n * singleton_count^2 / (2 * doubleton_count)
    } else {
      (n - 1) / n * singleton_count * (singleton_count - 1) / 2
    }
  }

  estimate_richness <- function(counts, sizes) {
    n <- sum(counts)
    observed_richness <- length(counts)
    unseen_richness <- estimate_unseen_richness(counts)
    singleton_count <- sum(counts == 1)

    count_frequency <- table(counts)
    frequencies <- as.numeric(names(count_frequency))
    frequency_multiplicity <- as.numeric(count_frequency)

    vapply(sizes, function(size) {
      if (size < n) {
        log_p_absent <- lchoose(n - frequencies, size) - lchoose(n, size)
        p_present <- pmin(pmax(-expm1(log_p_absent), 0), 1)
        estimate <- sum(p_present * frequency_multiplicity)
        return(min(max(estimate, 1), size, observed_richness))
      }

      if (size == n || unseen_richness <= 0 || singleton_count == 0) {
        return(as.numeric(observed_richness))
      }

      additional_size <- size - n
      chao_a <- n * unseen_richness / (n * unseen_richness + singleton_count)
      estimate <- observed_richness +
        unseen_richness * (-expm1(additional_size * log(chao_a)))
      min(max(estimate, observed_richness), size, observed_richness + unseen_richness)
    }, numeric(1))
  }

  make_bootstrap_sampler <- function(counts) {
    n <- sum(counts)
    singleton_count <- sum(counts == 1)
    unseen_richness <- estimate_unseen_richness(counts)

    if (unseen_richness <= 0 || singleton_count == 0) {
      observed_probabilities <- counts / n
      unseen_mass <- 0
      unseen_species <- 0
    } else {
      chao_a <- n * unseen_richness / (n * unseen_richness + singleton_count)
      unseen_mass <- singleton_count / n * chao_a
      empirical_probabilities <- counts / n
      non_detection <- exp(n * log1p(-empirical_probabilities))
      adjustment <- unseen_mass / sum(empirical_probabilities * non_detection)
      observed_probabilities <- empirical_probabilities *
        (1 - adjustment * non_detection)
      observed_probabilities <- pmax(observed_probabilities, 0)

      probability_total <- sum(observed_probabilities) + unseen_mass
      observed_probabilities <- observed_probabilities / probability_total
      unseen_mass <- unseen_mass / probability_total
      unseen_species <- ceiling(unseen_richness)
    }

    function() {
      unseen_draws <- if (unseen_mass > 0) {
        stats::rbinom(1, size = n, prob = unseen_mass)
      } else {
        0
      }
      observed_draws <- n - unseen_draws

      observed_counts <- numeric(0)
      if (observed_draws > 0) {
        probability_sum <- sum(observed_probabilities)
        if (probability_sum > 0) {
          observed_counts <- stats::rmultinom(
            1,
            size = observed_draws,
            prob = observed_probabilities / probability_sum
          )[, 1]
          observed_counts <- observed_counts[observed_counts > 0]
        }
      }

      unseen_counts <- numeric(0)
      if (unseen_draws > 0) {
        unseen_ids <- sample.int(unseen_species, unseen_draws, replace = TRUE)
        unseen_counts <- as.numeric(table(unseen_ids))
      }

      c(observed_counts, unseen_counts)
    }
  }

  checkmate::assert_logical(norm, len = 1)
  checkmate::assert_logical(verbose, len = 1)
  checkmate::assert_numeric(quantile, len = 2, any.missing = FALSE, finite = TRUE)
  quantile <- sort(as.numeric(quantile))
  if (quantile[1] <= 0 || quantile[2] >= 1) {
    cli::cli_abort("{.arg quantile} values must lie strictly between 0 and 1.")
  }
  if (quantile[1] == quantile[2]) {
    cli::cli_abort("{.arg quantile} values must be different.")
  }

  checkmate::assert_number(
    nboot,
    lower = 0,
    upper = .Machine$integer.max,
    finite = TRUE
  )
  if (nboot != round(nboot) || nboot == 1) {
    cli::cli_abort("{.arg nboot} must be 0 or an integer greater than or equal to 2.")
  }
  nboot <- as.integer(nboot)

  if (!is.na(step)) {
    checkmate::assert_number(step, lower = 1, finite = TRUE)
    if (step != round(step)) {
      cli::cli_abort("{.arg step} must be an integer.")
    }
    step <- as.numeric(round(step))
  }

  if (!is.na(extrapolation)) {
    checkmate::assert_number(extrapolation, lower = 0, finite = TRUE)
    if (extrapolation != round(extrapolation)) {
      cli::cli_abort("{.arg extrapolation} must be an integer.")
    }
    extrapolation <- as.numeric(round(extrapolation))
  }

  rep_col <- immundata::imd_schema("repertoire")
  cnt_col <- immundata::imd_schema("count")
  rep_sym <- immundata::imd_schema_sym("repertoire")
  rec_sym <- immundata::imd_schema_sym("receptor")
  cnt_sym <- immundata::imd_schema_sym("count")

  counts_tbl <- idata$annotations |>
    dplyr::select(!!rep_sym, !!rec_sym, !!cnt_sym) |>
    dplyr::distinct(!!rec_sym, !!rep_sym, .keep_all = TRUE) |>
    dplyr::select(!!rep_sym, !!cnt_sym)

  # Keep projection and deduplication lazy. The numerical estimator requires
  # abundance vectors, so only the two required columns cross the RAM boundary.
  counts_tbl <- counts_tbl |>
    dplyr::collect()

  if (!nrow(counts_tbl)) {
    out <- tibble::tibble(
      size = numeric(0),
      q_low = numeric(0),
      mean = numeric(0),
      q_high = numeric(0),
      type = character(0)
    )
    out[[rep_col]] <- numeric(0)
    return(out |>
      dplyr::select(all_of(c(rep_col, "size", "q_low", "mean", "q_high", "type"))) |>
      structure(rarefaction_norm = isTRUE(norm)))
  }

  counts_split <- counts_tbl |>
    dplyr::summarise(counts = list(!!rlang::sym(cnt_col)), .by = all_of(rep_col))

  counts_split$counts <- lapply(seq_len(nrow(counts_split)), function(i) {
    repertoire_id <- counts_split[[rep_col]][i]
    counts <- as.numeric(counts_split$counts[[i]])

    if (!length(counts) || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts <= 0)) {
      cli::cli_abort(
        "Counts for repertoire {.value {repertoire_id}} must be finite, positive clonotype counts."
      )
    }
    if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
      cli::cli_abort(
        "Counts for repertoire {.value {repertoire_id}} must be integer clonotype counts."
      )
    }

    counts <- as.numeric(round(counts))
    if (!is.finite(sum(counts))) {
      cli::cli_abort(
        "The total clone count for repertoire {.value {repertoire_id}} must be finite."
      )
    }
    counts
  })

  total_counts <- vapply(counts_split$counts, sum, numeric(1))
  default_sizes <- function(from, to) {
    unique(floor(seq(from, to, length.out = min(50, to - from + 1))))
  }
  if (is.na(step)) {
    shared_depth <- min(total_counts)
    shared_sizes <- default_sizes(1, shared_depth)
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info(
      "Computing rarefaction with {nboot} bootstrap replicate{?s} for {nrow(counts_split)} repertoire(s)."
    )
  }

  res_list <- lapply(seq_len(nrow(counts_split)), function(i) {
    repertoire_id <- counts_split[[rep_col]][i]
    counts <- counts_split$counts[[i]]
    n <- total_counts[i]
    observed_richness <- length(counts)

    max_bootstrap_size <- .Machine$integer.max
    if (nboot > 0 && n > max_bootstrap_size) {
      cli::cli_abort(
        "Bootstrap rarefaction supports at most {.value {max_bootstrap_size}} total clones; use {.code nboot = 0} for this repertoire."
      )
    }

    interpolation_sizes <- 1
    if (is.na(step)) {
      interpolation_sizes <- c(shared_sizes, default_sizes(shared_depth, n)[-1L])
    } else if (step <= n) {
      interpolation_sizes <- c(interpolation_sizes, seq(step, n, by = step))
    }
    interpolation_sizes <- sort(unique(c(interpolation_sizes, n)))

    extrapolation_limit <- if (is.na(extrapolation)) 2 * n else extrapolation
    extrapolation_sizes <- numeric(0)
    if (extrapolation_limit > n) {
      if (is.na(step)) {
        extrapolation_sizes <- default_sizes(n, extrapolation_limit)[-1L]
      } else if (n + step <= extrapolation_limit) {
        extrapolation_sizes <- seq(n + step, extrapolation_limit, by = step)
      }
      extrapolation_sizes <- sort(unique(c(
        extrapolation_sizes,
        extrapolation_limit
      )))
    }

    sizes <- c(interpolation_sizes, extrapolation_sizes)
    type <- c(
      rep("interpolation", length(interpolation_sizes)),
      rep("extrapolation", length(extrapolation_sizes))
    )
    estimate <- estimate_richness(counts, sizes)

    q_low <- rep(NA_real_, length(sizes))
    q_high <- rep(NA_real_, length(sizes))
    if (nboot > 0) {
      draw_bootstrap_counts <- make_bootstrap_sampler(counts)
      bootstrap_estimates <- matrix(
        vapply(seq_len(nboot), function(bootstrap_index) {
          estimate_richness(draw_bootstrap_counts(), sizes)
        }, numeric(length(sizes))),
        nrow = length(sizes),
        ncol = nboot
      )
      bootstrap_se <- apply(bootstrap_estimates, 1, stats::sd)

      q_low <- estimate + stats::qnorm(quantile[1]) * bootstrap_se
      q_high <- estimate + stats::qnorm(quantile[2]) * bootstrap_se
      q_low <- pmin(pmax(q_low, 1), sizes)
      q_high <- pmin(pmax(q_high, 1), sizes)
    }

    out <- tibble::tibble(
      size = as.numeric(sizes),
      q_low = as.numeric(q_low),
      mean = as.numeric(estimate),
      q_high = as.numeric(q_high),
      type = type
    )

    if (isTRUE(norm)) {
      out <- out |>
        dplyr::mutate(
          size = .data$size / n,
          q_low = .data$q_low / observed_richness,
          mean = .data$mean / observed_richness,
          q_high = .data$q_high / observed_richness
        )
    }

    out[[rep_col]] <- repertoire_id
    out |>
      dplyr::select(all_of(c(rep_col, "size", "q_low", "mean", "q_high", "type")))
  })

  dplyr::bind_rows(res_list) |>
    dplyr::arrange(!!rlang::sym(rep_col), .data$size, .data$type) |>
    structure(rarefaction_norm = isTRUE(norm))
}

#' @description
#' **3) Rarefaction and extrapolation (`airr_diversity_rarefaction`).** Estimate
#' how receptor richness changes with sample size. Use interpolation to compare
#' repertoires at a common sampling depth and extrapolation to estimate how many
#' additional receptors may be found with deeper sampling.
#'
#' Computation is done in RAM after one database query that fetches only the
#' repertoire identifier and receptor counts.
#' Interpolation uses the exact hypergeometric expectation for sampling without
#' replacement. Extrapolation uses the finite-sample Chao estimator, and
#' confidence bounds use bootstrap standard errors from an estimated abundance
#' distribution. Clonotype counts must be positive integers.
#'
#' @param step A whole number greater than or equal to `1`, or `NA`. The default,
#'   `NA`, uses up to 50 evenly spaced points through the smallest observed
#'   depth. Each repertoire gets up to 50 points for each remaining interpolation
#'   and extrapolation interval.
#' @param quantile A numeric vector of length 2 containing confidence
#'   probabilities strictly between `0` and `1`.
#' @param extrapolation A nonnegative whole number or `NA` giving the upper size
#'   limit for extrapolation. Use `0` to disable extrapolation. The default,
#'   `NA`, uses twice each repertoire's observed clone count.
#' @param nboot A nonnegative whole number that is either `0` or at least `2`.
#'   It gives the number of bootstrap replicates used to estimate confidence
#'   bounds. Defaults to `50`; use `0` to return `NA` confidence bounds.
#' @param norm A single logical value; if `TRUE`, size is divided by total
#'   repertoire size and richness estimates are divided by observed richness.
#' @param verbose A single logical value controlling whether a concise progress
#'   message is shown.
#'
#' @references
#' Chao, A. et al. (2014). Rarefaction and extrapolation with Hill numbers:
#' a framework for sampling and estimation in species diversity studies.
#' *Ecological Monographs*, 84, 45-67. \doi{10.1890/13-0133.1}
#'
#' @return
#'
#' ## 3) Rarefaction and extrapolation (`airr_diversity_rarefaction`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `size` - sample size (absolute or normalised)
#' * `q_low` - lower confidence bound
#' * `mean` - expected richness
#' * `q_high` - upper confidence bound
#' * `type` - `interpolation` or `extrapolation`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' #
#' # Calculate a rarefaction curve without extrapolation.
#' raref <- airr_diversity_rarefaction(immdata, step = 2000)
#'
#' # Visualise the rarefaction curve.
#' vis(raref)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_rarefaction <- register_method_airr(airr_diversity_rarefaction_impl, "airr_diversity", "rarefaction")


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


#' @description
#' **4) Shannon entropy (`airr_diversity_shannon`).** Calculate Shannon entropy
#' from receptor proportions in each repertoire. The value increases when a
#' repertoire contains more receptors or has a more even abundance
#' distribution. Use it as one summary that reflects both richness and
#' evenness.
#'
#' @return
#'
#' ## 4) Shannon entropy (`airr_diversity_shannon`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `shannon` - entropy in bits
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' #
#' # Calculate Shannon entropy.
#' sh <- airr_diversity_shannon(immdata)
#' vis(sh)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_shannon <- register_method_airr(airr_diversity_shannon_impl, "airr_diversity", "shannon")


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


#' @description
#' **5) Pielou evenness (`airr_diversity_pielou`).** Calculate how evenly
#' receptor abundance is distributed within each repertoire. The method divides
#' Shannon entropy by the maximum entropy for the observed richness. Use it to
#' compare evenness between repertoires with different numbers of receptors.
#'
#' @return
#'
#' ## 5) Pielou evenness (`airr_diversity_pielou`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `shannon`
#' * `n_receptors`
#' * `pielou` - evenness in `[0, 1]` (`NA` if `S <= 1`)
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' #
#' # Calculate Pielou evenness.
#' pj <- airr_diversity_pielou(immdata)
#' vis(pj)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_pielou <- register_method_airr(airr_diversity_pielou_impl, "airr_diversity", "pielou")


#' @keywords internal
airr_diversity_index_impl <- function(idata) {
  airr_diversity_hill_impl(idata, q = 1)
}


#' @description
#' **6) Diversity index (`airr_diversity_index`).** Calculate the Hill number at
#' `q = 1`, which is the exponential of Shannon entropy calculated with natural
#' logarithms. Use it as an effective number of receptors that is easy to
#' compare between repertoires.
#'
#' @return
#'
#' ## 6) Diversity index (`airr_diversity_index`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `q = 1`
#' * `hill_number`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' #
#' # Calculate the default diversity index.
#' idx <- airr_diversity_index(immdata)
#' vis(idx)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_index <- register_method_airr(airr_diversity_index_impl, "airr_diversity", "index")


#' @keywords internal
airr_diversity_hill_impl <- function(idata, q = 0:5) {
  checkmate::assert_numeric(
    q,
    lower = 0,
    min.len = 1,
    any.missing = FALSE,
    finite = TRUE,
    sorted = TRUE
  )

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

  result |>
    collect()
}


#' @description
#' **7) Hill numbers (`airr_diversity_hill`).** Calculate a diversity profile
#' for one or more orders `q`. At `q = 0`, the result is receptor richness; at
#' `q = 1`, it is the exponential of Shannon entropy; higher values of `q` give
#' more weight to abundant receptors. Use several orders to compare the
#' influence of rare and abundant receptors.
#'
#' @inheritParams method_args_airr
#' @param q A numeric scalar or vector of nonnegative orders. Defaults to `0:5`.
#'
#' @return
#'
#' ## 7) Hill numbers (`airr_diversity_hill`)
#' A tibble with:
#' * `imd_repertoire_id`
#' * `q` - Hill order
#' * `hill_number` - true diversity of order `q`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' #
#' # Calculate a diversity profile for three Hill orders.
#' hill <- airr_diversity_hill(immdata, q = c(0, 1, 2))
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_hill <- register_method_airr(airr_diversity_hill_impl, "airr_diversity", "hill")
