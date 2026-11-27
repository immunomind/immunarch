#' One function to visualise them all
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom tidyr drop_na
#'
#' @description
#'
#' `r lifecycle::badge('experimental')`
#'
#' **`vis()` is a lightweight, quick-look plotting helper.**
#' It's designed to help you **visualise results fast** with sensible defaults.
#' It automatically detects the input type and chooses an appropriate visualisation
#' (e.g., output from `airr_desc_genes()` is recognised as gene usage values and plotted
#' without extra arguments).
#'
#' **vis() is not intended for publication-quality figures.**
#' For serious, highly customised, or publication-ready plots, I recommend
#' building your graphics directly with **ggplot2**.
#'
#' Information on how visualise a specific result from a domain, e.g., clonality or diversity,
#' is available at each domain's help page in the visualisation section.
#' Take a look at the "See also" section down below to *see* links to the domains, which support visualisations.
#'
#' @param .data The output from any immunarch analysis function. The function automatically resolves to a correct visualisation.
#' @param ... Any other arguments, see the "Details" section for specific visualisation functions.
#'
#' @details
#' `vis()` chooses a plot based on the analysis result. Start with `vis(result)`,
#' then consult the documentation for the analysis that created it to see its
#' plotting defaults and supported arguments.
#'
#' E.g., for clonality results, see [airr_clonality] by running `?airr_clonality`.
#' Or see [repsim] For repertoire similarity by running `?repsim`.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso
#' * [airr_desc] for descriptive statistics visualisations.
#' * [airr_clonality] for clonality visualisations.
#' * [airr_diversity] for diversity visualisations.
#' * [repsim] for repertoire similarity visualisations.
#' * [dist] for receptor distance visualisations.
#' * [vdj_track] for receptor tracking visualisations with `output = "df"`.
#' * [vdj_public] for selected receptor coverage and abundance. Pass its
#'   `output = "df"` result and the original `reference` to `vis()`, or use
#'   `plot = "public"` for its `ImmunData` result.
#' * [annotate_public] for incidence distributions with `plot = "incidence"`,
#'   stacked by stratum when strata are defined and unstacked otherwise.
#'
#' @examples
#' \dontrun{
#' immdata <- get_test_idata() |> agg_repertoires("Therapy")
#' airr_desc_genes(immdata, gene_col = "v_call") |> vis()
#' repsim_jaccard(immdata) |> vis()
#' airr_diversity_pielou(immdata) |> vis()
#' airr_diversity_chao1(immdata) |> vis()
#' # Tracking summaries include explicit zeros for unobserved combinations.
#' vdj_track(immdata, n = 5, output = "df") |> vis()
#' vdj_public(immdata, output = "df") |> vis(reference = immdata)
#' annotate_public(immdata) |> vis(plot = "incidence")
#' }
#' @export
vis <- function(.data, ...) {
  UseMethod("vis")
}
