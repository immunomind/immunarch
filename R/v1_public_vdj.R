#' @keywords internal
vdj_public_impl <- function(idata, ...) {
  rules <- list(...)
  if (!length(rules)) rules <- list(make_public_filter(min_incidence = 2L))
  repertoires <- dplyr::collect(idata$repertoires)
  occurrences <- idata$annotations |>
    dplyr::select(imd_receptor_id, imd_repertoire_id, imd_count, imd_proportion) |>
    dplyr::distinct()
  selected <- occurrences |> dplyr::distinct(imd_receptor_id)

  for (rule in rules) {
    checkmate::assert_list(rule)
    rule <- do.call(make_public_filter, rule)
    population <- repertoires
    if (!is.null(rule$repertoires)) {
      unknown <- setdiff(rule$repertoires, repertoires$imd_repertoire_id)
      if (length(unknown)) cli::cli_abort("Unknown repertoire IDs: {.val {unknown}}.")
      population <- population[population$imd_repertoire_id %in% rule$repertoires, ]
    }
    if (!is.null(rule$strata) || rule$incidence_by == "strata") {
      if (!"imd_strata_id" %in% names(repertoires)) {
        cli::cli_abort("Strata must be defined with {.fn immundata::agg_strata} first.")
      }
    }
    if (!is.null(rule$strata)) {
      column <- if (is.numeric(rule$strata)) "imd_strata_id" else "strata_name"
      unknown <- setdiff(rule$strata, repertoires[[column]])
      if (length(unknown)) cli::cli_abort("Unknown strata: {.val {unknown}}.")
      population <- population[population[[column]] %in% rule$strata, ]
    }
    if (!nrow(population)) cli::cli_abort("Empty comparison population.")

    unit <- if (rule$incidence_by == "strata") "imd_strata_id" else "imd_repertoire_id"
    population <- population |>
      dplyr::select(dplyr::all_of(unique(c("imd_repertoire_id", unit)))) |>
      duckplyr::as_duckdb_tibble()
    qualifying <- occurrences |> dplyr::inner_join(population, by = "imd_repertoire_id")
    if (!is.null(rule$min_count)) qualifying <- qualifying |> dplyr::filter(imd_count >= !!rule$min_count)
    if (!is.null(rule$max_count)) qualifying <- qualifying |> dplyr::filter(imd_count <= !!rule$max_count)
    if (!is.null(rule$min_proportion)) qualifying <- qualifying |> dplyr::filter(imd_proportion >= !!rule$min_proportion)
    if (!is.null(rule$max_proportion)) qualifying <- qualifying |> dplyr::filter(imd_proportion <= !!rule$max_proportion)
    incidence <- qualifying |>
      dplyr::distinct(imd_receptor_id, !!rlang::sym(unit)) |>
      dplyr::count(imd_receptor_id, name = ".incidence")
    selected <- selected |>
      dplyr::left_join(incidence, by = "imd_receptor_id") |>
      dplyr::mutate(.incidence = dplyr::coalesce(.incidence, 0L))
    minimum <- rule$min_incidence
    if (is.null(minimum) && is.null(rule$max_incidence)) minimum <- 1L
    if (!is.null(minimum)) selected <- selected |> dplyr::filter(.incidence >= !!minimum)
    if (!is.null(rule$max_incidence)) selected <- selected |> dplyr::filter(.incidence <= !!rule$max_incidence)
    selected <- selected |> dplyr::select(imd_receptor_id)
  }

  keys <- c("imd_receptor_id", "imd_repertoire_id")
  abundances <- occurrences |>
    dplyr::semi_join(selected, by = "imd_receptor_id") |>
    dplyr::rename(imd_orig_count = imd_count, imd_orig_proportion = imd_proportion)
  immundata::annotate_immundata(
    idata,
    annotations = abundances, by = stats::setNames(keys, keys),
    conflicts = "replace", remove_limit = TRUE
  ) |>
    immundata::filter_immundata(!is.na(imd_orig_count))
}


#' @title Select receptors by their incidence across repertoires or strata
#'
#' @description
#' Select shared receptors, receptors private to a strata, or receptors absent
#' above an abundance threshold. Define one or more rules to select receptors
#' by their sharing and abundance in comparison populations using [make_public_filter()].
#'
#' @param idata An [immundata::ImmunData] object with defined repertoires.
#' @param ... One or more rules created by [make_public_filter()]. With no rules, require presence
#'   in at least two input repertoires.
#' @param output A character scalar: `"immundata"` returns a new `ImmunData`,
#'   `"df"` returns a data frame convenient to immediately pass to  [vis()].
#'
#' @details
#' Input counts and proportions are saved
#' as `imd_orig_count` and `imd_orig_proportion`.
#' The input is unchanged. No per-rule incidence columns or custom attributes
#' are added.
#'
#' @section Repeated calls:
#' **Warning:** Applying `vdj_public()` to an already filtered result overwrites
#' `imd_orig_count` and `imd_orig_proportion` with abundances from its immediate
#' input. Since proportions are recalculated after filtering, the saved
#' proportions no longer refer to the repertoire before the first call.
#' Keep the original `ImmunData` as a reference when calculating percentages
#' relative to the original repertoires. To apply several rules to the same
#' original input, pass them together in one `vdj_public()` call.
#'
#' @section Visualisation:
#' Plot the selected receptors with
#' `vis(selected, plot = "public", reference = idata)`, where `idata` is the
#' original input to `vdj_public()`. Alternatively, pipe
#' `vdj_public(idata, output = "df")` to `vis(reference = idata)`; its result
#' class selects this plot automatically. The ordinary observation rows in
#' the collected output are unchanged.
#'
#' Bars show **abundance**, the sum of the selected identities' saved original
#' repertoire proportions, read against the left percentage axis. Dots show
#' **coverage**, the number of selected identities detected in each repertoire,
#' read against the right receptor-count axis. Dot labels give the
#' detected/selected counts; bar labels give abundance percentages.
#' Observed occurrences below selection thresholds still contribute to both
#' measures. Coverage describes the selected set, not the fraction of all
#' receptor identities in a repertoire. Individual receptor trajectories can
#' be visualised with [vdj_track()].
#'
#' The required `reference` defines the complete input repertoire population,
#' retaining repertoires with no selected receptors at zero. Use the selected
#' result and its own original reference: repertoire metadata is matched using
#' the reference repertoire schema because filtering can rebuild repertoire
#' IDs. The reference does not restore abundances overwritten by repeated
#' `vdj_public()` calls; the saved original proportions must refer to this
#' baseline. An empty selection cannot be plotted.
#'
#' `xval` selects a repertoire ID or metadata column for the x-axis; each value
#' must identify one repertoire. `fill` selects a repertoire metadata column
#' for colour and defaults to the defined strata when available. Repertoires
#' are ordered by colour group and x value; factor levels control their order.
#' `title` changes the title or removes it with `NULL`.
#'
#' @return A new `ImmunData` or collected observation table.
#'
#' @seealso [make_public_filter()], [annotate_public()]
#'
#' @examples
#' idata <- get_test_idata() |> agg_repertoires("Therapy")
#'
#' # Keep receptors shared across at least two repertoires.
#' shared <- vdj_public(idata)
#' shared$repertoires
#' vis(shared, plot = "public", reference = idata)
#'
#' \dontrun{
#' # Work with repertoire and strata rules.
#' idata <- get_test_idata() |>
#'   agg_repertoires(c("Response", "Therapy")) |>
#'   agg_strata("Therapy")
#'
#' vdj_public(idata, output = "df") |>
#'   vis(reference = idata, fill = "Therapy")
#'
#' # Inspect strata and find the IDs for the two therapies.
#' groups <- idata$strata
#' groups
#' ici <- groups$imd_strata_id[groups$Therapy == "ICI"]
#' car_t <- groups$imd_strata_id[groups$Therapy == "CAR-T"]
#'
#' # The test data has one repertoire per therapy, so use incidence >= 1.
#' in_ici <- vdj_public(idata, make_public_filter(
#'   strata = ici, min_incidence = 1L
#' ))
#'
#' # Count distinct strata instead of repertoires.
#' across_therapies <- vdj_public(idata, make_public_filter(
#'   incidence_by = "strata", min_incidence = 2L
#' ))
#'
#' # Present in ICI and absent in CAR-T. Collect the retained observations.
#' ici_only <- vdj_public(
#'   idata,
#'   make_public_filter(strata = ici, min_incidence = 1L),
#'   make_public_filter(strata = car_t, max_incidence = 0L),
#'   output = "df"
#' )
#' ici_only |>
#'   dplyr::select(
#'     cdr3_aa, v_call, Therapy, imd_orig_count,
#'     imd_orig_proportion
#'   ) |>
#'   head()
#'
#' # Select by repertoire IDs and require count >= 5 in an ICI repertoire.
#' ids <- idata$repertoires |>
#'   dplyr::filter(Therapy == "ICI") |>
#'   dplyr::pull(imd_repertoire_id)
#' abundant <- vdj_public(idata, make_public_filter(
#'   repertoires = ids, min_count = 5
#' ))
#' # All observations of passing receptors are retained throughout idata.
#' }
#' @export
vdj_public <- register_method_vdj(
  core = vdj_public_impl,
  family = "vdj",
  name = "public",
  result_class = TRUE
)


#' @title Describe one publicness selection rule
#'
#' @description
#' Define a comparison population and inclusive incidence and abundance bounds
#' for [vdj_public()]. The constructor stores a rule without evaluating data.
#'
#' @param repertoires A numeric vector of repertoire IDs, or `NULL` for all.
#' @param strata Numeric strata IDs or existing `strata_name` labels, or `NULL`
#'   for all. When both selectors are supplied, use their intersection.
#' @param incidence_by A character scalar: `"repertoire"` counts distinct
#'   repertoires; `"strata"` counts existing strata.
#' @param min_incidence,max_incidence Non-negative whole-number bounds, or
#'   `NULL`. Without either bound, the minimum is one. An explicit maximum
#'   alone imposes no minimum, so zero incidence can pass.
#' @param min_count,max_count Numeric count bounds, or `NULL`.
#' @param min_proportion,max_proportion Numeric within-repertoire proportion
#'   bounds, or `NULL`.
#'
#' @details
#' Abundance bounds qualify individual receptor occurrences before incidence
#' is counted. Repeated chain/cell abundance values are deduplicated, never
#' summed. For example, `max_count = 10, min_incidence = 2` requires qualifying
#' occurrences in two repertoires, allowing larger occurrences elsewhere.
#' `min_count = 5, max_incidence = 0` allows observed counts below five.
#' Unknown selectors and empty comparison populations are errors during
#' evaluation. Selectors do not match arbitrary metadata values.
#'
#' @return A list describing one rule for [vdj_public()].
#'
#' @seealso [vdj_public()] for strata preparation and applying several rules.
#'
#' @examples
#' idata <- get_test_idata() |>
#'   agg_repertoires(c("Response", "Therapy"))
#'
#' # Count >= 5 in at least one repertoire.
#' abundant <- make_public_filter(min_count = 5)
#' vdj_public(idata, abundant, output = "df") |>
#'   dplyr::select(cdr3_aa, v_call, Therapy, imd_orig_count) |>
#'   head()
#'
#' # At least one occurrence with count <= 10. Larger occurrences are allowed.
#' low_count <- make_public_filter(max_count = 10)
#'
#' # No occurrence with count >= 5. Smaller counts are allowed.
#' below_threshold <- make_public_filter(min_count = 5, max_incidence = 0L)
#'
#' @export
make_public_filter <- function(
  repertoires = NULL,
  strata = NULL,
  incidence_by = c("repertoire", "strata"),
  min_incidence = NULL,
  max_incidence = NULL,
  min_count = NULL,
  max_count = NULL,
  min_proportion = NULL,
  max_proportion = NULL
) {
  checkmate::assert_integerish(repertoires, any.missing = FALSE, null.ok = TRUE)
  if (is.numeric(strata)) {
    checkmate::assert_integerish(strata, any.missing = FALSE)
  } else {
    checkmate::assert_character(strata, any.missing = FALSE, null.ok = TRUE)
  }
  for (bound in list(min_incidence, max_incidence)) {
    checkmate::assert_count(bound, null.ok = TRUE)
  }
  for (bound in list(min_count, max_count)) {
    checkmate::assert_number(bound, lower = 0, finite = TRUE, null.ok = TRUE)
  }
  for (bound in list(min_proportion, max_proportion)) {
    checkmate::assert_number(bound, lower = 0, upper = 1, null.ok = TRUE)
  }
  for (bounds in list(
    c(min_incidence, max_incidence), c(min_count, max_count),
    c(min_proportion, max_proportion)
  )) {
    if (length(bounds) == 2L && bounds[1] > bounds[2]) {
      cli::cli_abort("Minimum bounds must not exceed maximum bounds.")
    }
  }
  list(
    repertoires = repertoires, strata = strata,
    incidence_by = match.arg(incidence_by),
    min_incidence = min_incidence, max_incidence = max_incidence,
    min_count = min_count, max_count = max_count,
    min_proportion = min_proportion, max_proportion = max_proportion
  )
}
