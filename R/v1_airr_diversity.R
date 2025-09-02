#' @title Diversity — estimating the heterogeneity of immune repertoires
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions to quantify **receptor diversity** per repertoire.
#'
#' ## Available functions:
#'
#' @param idata An `ImmunData` object.
#' @inheritParams airr_diversity_dxx
#' @inheritParams airr_diversity_shannon
#' @inheritParams airr_diversity_pielou
#' @inheritParams airr_diversity_hill
#' @inheritParams airr_diversity_index
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @name airr_diversity
#' @concept Diversity
NULL


#' @keywords internal
airr_diversity_dxx_impl <- function(idata, perc = 50) {
  checkmate::assert_numeric(perc, null.ok = TRUE)
}


#' @description `airr_diversity_dxx` — **coverage diversity**: minimal number of
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
#' * `repertoire_id`
#' * `perc`
#' * `dxx` — minimal count of top receptors to reach `perc%`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_dxx
#' #
#' d50 <- airr_diversity_dxx(immdata, perc = 50)
#' d_multi <- airr_diversity_dxx(immdata, perc = c(20, 50, 80))
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_dxx <- register_immunarch_method(airr_diversity_dxx_impl, "airr_diversity", "dxx")


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
    )
}


#' @description `airr_diversity_shannon` — Shannon entropy (base 2) per repertoire
#' computed from `proportion`. Ideal when you want a single evenness-aware
#' diversity score; pair with Pielou/Hill for samples with very different richness.
#'
#' @return
#'
#' ## `airr_diversity_shannon`
#' A tibble with:
#' * `repertoire_id`
#' * `shannon` — entropy in bits
#'
#' @examples
#' #
#' # airr_diversity_shannon
#' #
#' sh <- airr_diversity_shannon(immdata)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_shannon <- register_immunarch_method(airr_diversity_shannon_impl, "airr_diversity", "shannon")


#' @keywords internal
airr_diversity_pielou_impl <- function(idata) {
  shannon_values <- airr_diversity_shannon(idata)

  idata$repertoires |>
    left_join(shannon_values,
      by = imd_schema("repertoire")
    ) |>
    mutate(richness = dd$log2(!!immundata::imd_schema_sym("n_receptors")), pielou = shannon / richness)
}


#' @description `airr_diversity_pielou` — Pielou’s evenness `H / log2(S)` with
#' richness `S`. Best when you need a **size-normalized** evenness score that’s
#' comparable across repertoires with different receptor counts.
#'
#' @return
#'
#' ## `airr_diversity_pielou`
#' A tibble with:
#' * `repertoire_id`
#' * `shannon`
#' * `n_receptors`
#' * `pielou` — evenness in `[0, 1]` (NA if `S ≤ 1`)
#'
#' @examples
#' #
#' # airr_diversity_pielou
#' #
#' pj <- airr_diversity_pielou(immdata)
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_pielou <- register_immunarch_method(airr_diversity_pielou_impl, "airr_diversity", "pielou")


#' @keywords internal
airr_diversity_index_impl <- function(idata) {
  airr_diversity_hill(idata, q = 1)
}


#' @description `airr_diversity_index` — convenience alias for Hill number with
#' `q = 1` (`exp(Shannon)` using natural log). A solid **default single metric**
#' that’s relatively robust to rare-count noise and easy to compare across samples.
#'
#' @return
#'
#' ## `airr_diversity_index`
#' A tibble with:
#' * `repertoire_id`
#' * `q = 1`
#' * `hill_number`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_index
#' #
#' idx <- airr_diversity_index(immdata)
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

  idata$repertoires |> left_join(result, by = imd_schema("repertoire"))
}


#' @description `airr_diversity_hill` — Hill numbers (“true diversity”) for
#' orders `q ∈ {0, 1, 2, …}`: `q=0` richness, `q=1` exp(Shannon), `q>1`
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
#' * `repertoire_id`
#' * `q` — Hill order
#' * `hill_number` — true diversity of order `q`
#' * plus repertoire metadata from `idata$repertoires`
#'
#' @examples
#' #
#' # airr_diversity_hill
#' #
#' hill <- airr_diversity_hill(immdata, q = c(0, 1, 2))
#'
#' @rdname airr_diversity
#' @concept Diversity
#' @export
airr_diversity_hill <- register_immunarch_method(airr_diversity_hill_impl, "airr_diversity", "hill")
