#' @keywords internal
vdj_track_impl <- function(idata, from = NULL, n = NULL, output = c("immundata", "df")) {
  output <- match.arg(output)
  repertoire_source <- is.null(from) || is.numeric(from)
  if (!is.null(n)) {
    checkmate::assert_count(n)
    if (!repertoire_source) {
      cli::cli_abort("{.arg n} applies only to repertoire sources; select source receptors upstream instead.")
    }
  }

  if (repertoire_source) {
    source <- idata
    if (!is.null(from)) {
      checkmate::assert_integerish(from, min.len = 1L, any.missing = FALSE)
      unknown <- setdiff(from, idata$repertoires$imd_repertoire_id)
      if (length(unknown)) {
        cli::cli_abort("Unknown repertoire IDs in {.arg from}: {.val {unknown}}.")
      }
      source <- immundata::filter_immundata(source, imd_repertoire_id %in% !!from)
    }
    if (!is.null(n)) source <- immundata::top_n(source, n = n)
    selected <- source$receptors |> dplyr::distinct(imd_receptor_id)
  } else if (inherits(from, "ImmunData")) {
    if (!identical(idata$schema_receptor, from$schema_receptor)) {
      cli::cli_abort("Source and target must use the same receptor features and chain schema.")
    }
    target_receptors <- idata$receptors
    features <- setdiff(names(target_receptors), "imd_receptor_id")
    selected <- target_receptors |>
      dplyr::semi_join(from$receptors |> dplyr::select(-imd_receptor_id), by = features) |>
      dplyr::distinct(imd_receptor_id)
  } else if (is.data.frame(from)) {
    if (!ncol(from)) cli::cli_abort("The source query must contain at least one column.")
    checkmate::assert_names(names(from), type = "unique")
    unknown <- setdiff(names(from), names(idata$annotations))
    if (length(unknown)) {
      cli::cli_abort("Unknown query columns: {.field {unknown}}.")
    }
    query <- duckplyr::as_duckdb_tibble(from)
    missing <- query |>
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ any(is.na(.x)))) |>
      dplyr::collect()
    if (any(unlist(missing, use.names = FALSE), na.rm = TRUE)) {
      cli::cli_abort("Missing query values are not allowed; omit a column to leave it unconstrained.")
    }
    selected <- idata$annotations |>
      dplyr::semi_join(query, by = names(query)) |>
      dplyr::distinct(imd_receptor_id)
  } else {
    cli::cli_abort("{.arg from} must be NULL, repertoire IDs, ImmunData, or a query table.")
  }

  keys <- c("imd_receptor_id", "imd_repertoire_id")
  abundances <- idata$annotations |>
    dplyr::semi_join(selected, by = "imd_receptor_id") |>
    dplyr::select(dplyr::all_of(keys),
      imd_orig_count = imd_count, imd_orig_proportion = imd_proportion
    ) |>
    dplyr::distinct()
  if (output == "df") {
    receptors <- idata$receptors |>
      dplyr::semi_join(selected, by = "imd_receptor_id")
    features <- setdiff(names(receptors), "imd_receptor_id")
    if (length(immundata::imd_receptor_chains(idata$schema_receptor)) == 2L) {
      first <- c("locus.x", setdiff(grep("\\.x$", features, value = TRUE), "locus.x"))
      second <- c("locus.y", setdiff(grep("\\.y$", features, value = TRUE), "locus.y"))
      label <- rlang::expr(dd$concat_ws(
        " + ",
        dd$concat_ws(" ", !!!rlang::syms(first)),
        dd$concat_ws(" ", !!!rlang::syms(second))
      ))
    } else {
      label <- rlang::expr(dd$concat_ws(" ", !!!rlang::syms(features)))
    }
    receptors <- receptors |> dplyr::mutate(imd_receptor_label = !!label)
    # Receptor identities and the full repertoire axis come from the input.
    # Complete abundances before collecting; no cell/chain observations are added.
    return(receptors |>
      # A constant-key join keeps the Cartesian product in DuckDB;
      # duckplyr does not yet implement cross_join() relationally.
      dplyr::mutate(imd_track_join = 1L) |>
      dplyr::inner_join(
        duckplyr::as_duckdb_tibble(idata$repertoires) |>
          dplyr::mutate(imd_track_join = 1L),
        by = "imd_track_join", suffix = c("_receptor", "")
      ) |>
      dplyr::select(-imd_track_join) |>
      dplyr::left_join(abundances, by = keys) |>
      dplyr::mutate(
        imd_orig_count = dplyr::coalesce(imd_orig_count, 0),
        imd_orig_proportion = dplyr::coalesce(imd_orig_proportion, 0)
      ) |>
      dplyr::arrange(imd_repertoire_id, imd_receptor_id) |>
      dplyr::collect())
  }
  immundata::annotate_immundata(
    idata,
    annotations = abundances, by = stats::setNames(keys, keys),
    conflicts = "replace", remove_limit = TRUE
  ) |>
    immundata::filter_immundata(!is.na(imd_orig_count))
}


#' @title Track selected immune receptors across time and space
#'
#' @description
#' Use `vdj_track()` to follow receptors across samples, time points, tissues,
#' or cell populations. For example, select abundant receptors before treatment
#' and find their abundances in later samples. You can also follow a set of
#' sequences identified in another experiment.
#'
#' Select receptors with the argument `from`, or use `n` to select the most abundant
#' receptors in each source repertoire. Follow them throughout the target
#' `idata`, including repertoires where they are rare or unobserved.
#' The input `ImmunData` object is not changed.
#'
#' @param idata An [immundata::ImmunData] object with defined repertoires.
#'   These are the repertoires through which receptors will be followed.
#' @param from `NULL`, a numeric vector of repertoire IDs, an
#'   [immundata::ImmunData] object, or a data frame.
#'
#'   `NULL`, the default, uses every target repertoire as a source.
#'   Numeric values select source `imd_repertoire_id` values in `idata`.
#'
#'   An `ImmunData` source supplies its observed receptor identities.
#'
#'   A data frame supplies exact column combinations to match, such as `junction`
#'   alone or `junction` together with `v_call`. Tibbles and duckplyr tables
#'   are also accepted. [tibble::tibble] is recommended.
#' @param n A non-negative whole number or `NULL`. Select up to `n` receptors
#'   with the highest counts in each source repertoire, then follow their
#'   combined set of identities. Equal counts are resolved by receptor ID, as
#'   in [immundata::top_n()]. `NULL`, the default, selects all source receptors.
#'   This argument applies only when `from` is `NULL` or contains repertoire
#'   IDs.
#' @param output A character scalar. `"immundata"`, the default, returns a new
#'   `ImmunData` object containing the selected receptors' observed cells and
#'   chains. `"df"` returns an abundance-summary tibble with one row per
#'   selected receptor and input repertoire, including explicit zeros where
#'   unobserved. Both outputs preserve the original target abundances.
#'   The summary retains input receptor and repertoire IDs, receptor identity
#'   features, and repertoire metadata as columns, and ready for plotting with [vis()].
#'
#' @section Common use cases:
#' * **Follow receptors from one or several samples.** Supply their repertoire
#'   IDs through `from`. Receptors are followed in every target repertoire,
#'   including those outside the source selection.
#' * **Follow the top receptors from one sample.** Supply a repertoire ID and
#'   `n`, for example `from = 5L, n = 3L` for the top three from repertoire ID 5.
#' * **Follow the top receptors from every sample.** Supply `n` and leave
#'   `from = NULL`. A receptor selected in any source is followed everywhere,
#'   even where it is below the top `n`. The combined set can contain more than
#'   `n` identities.
#' * **Follow an upstream selection.** Pass an `ImmunData` object containing
#'   receptors of interest, such as expanded clones or a selected cell
#'   population. The source can also come from an independent import.
#' * **Follow a sequence panel.** Pass a table containing sequences and optional
#'   gene or locus columns. All columns in a row must match the same chain;
#'   rows are alternatives. Matching one chain keeps the complete receptor,
#'   including its partner for paired receptors.
#'
#' @details
#' An `ImmunData` source must use the same receptor features and chain schema
#' as the target. Complete feature combinations, including both chains of a
#' paired receptor, are matched
#'
#' Column names same as in the `idata` must exist in the target observations in the `from` argument.
#' Matching is exact and case-sensitive, with no sequence or gene-call normalization.
#' Missing query values in `from` are not allowed. Unmatched
#' queries, empty query tables, and `n = 0` produce an empty result.
#'
#' @section Original abundances:
#' `imd_orig_count` and `imd_orig_proportion` columns in the output store the target
#' receptors counts and proportions before `vdj_track` filters the data. Use these columns to compare
#' abundances across the original samples. The ordinary `imd_count` and
#' `imd_proportion` describe the filtered object, following normal `immundata`
#' behavior, i.e., normalising the proportions.
#' Please be careful - another `vdj_track` call replaces the original-abundance columns
#' `imd_orig_count` and `imd_orig_proportion`.
#'
#' The `ImmunData` result contains observed chains only and follows normal
#' `immundata` filtering, including rebuilding IDs for the remaining repertoires.
#' In other words, the repertoire ID could change after calling `vdj_track` with
#' `output = "immundata"`. This is temporary inconvenience, sorry for that. The good news
#' is that no repertoire metadata is missing so you still analyse the data.
#' The `"df"` summary uses the input IDs and contains zero-count rows for
#' unobserved combinations. These are abundance rows, without barcode or
#' observation-specific fields. A receptor not observed in a sample is not
#' necessarily biologically absent.
#'
#' @section Visualisation:
#' Pass `vdj_track(..., output = "df")` directly to [vis()]. The summary
#' tibble carries a result class that selects the tracking plot. The default
#' `ImmunData` output is intended for further analysis, so pass `output = "df"`
#' whenever you want to plot the results.
#'
#' `vis()` draws an alluvial plot of `imd_orig_proportion` by
#' `imd_repertoire_id` by default. Receptors are distinguished by colour and
#' labelled by their receptor features in the legend, for example `AAA V01`.
#' For paired receptors, each chain's locus and features are combined, for
#' example `IGH AAA V01 + IGL BBB V02`. Receptor IDs remain the grouping keys,
#' so pairs sharing one chain remain separate trajectories. The plot accepts:
#'
#' * `plot = "alluvial"` (the default) draws stacked bars connected by smooth
#'   ribbons and requires the optional package **ggalluvial**.
#' * `plot = "area"` draws stacked areas.
#' * `plot = "line"` draws individual trajectories with
#'   points.
#' * `xval` selects the repertoire ID or metadata column for the x-axis. The
#'   default is to use the numeric repertoire IDs. This allows you to select
#'   metadata columns such as "Timepoint".
#' * `yval` selects `"imd_orig_proportion"` (the default) or `"imd_orig_count"`.
#'   These preserve abundances from the original target repertoires. Proportions
#'   are not renormalised to the selected receptors and need not sum to 100%.
#' * `facet` selects one metadata column for separate panels. Each x value
#'   within a panel must identify exactly one repertoire. Define repertoires
#'   upstream to pool samples; plotting does not aggregate them.
#' * `title` replaces the default title, `"Receptor tracking"`.
#'
#' The summary already contains zero abundances for unobserved receptor-repertoire
#' combinations. Every original target repertoire remains on the plot, including samples
#' with no tracked observations.
#'
#' Repertoire metadata is returned in regular columns for every target sample,
#' including samples with no tracked observations. Use those columns for
#' `xval` and `facet`, or prepare additional display columns in the summary.
#'
#' @return With `output = "immundata"`, a new [immundata::ImmunData] object
#'   containing the selected receptors' observed rows, with `imd_orig_count`
#'   and `imd_orig_proportion`.
#'
#'   With `output = "df"`, an immunarch result tibble containing:
#'   * `imd_receptor_id` and receptor identity features from the input receptor
#'     dictionary. Paired-chain features are shown side by side.
#'   * `imd_receptor_label`, a readable label combining the receptor features.
#'     Paired-chain labels include both loci and use ` + ` between chains.
#'   * `imd_repertoire_id` and metadata from the input repertoire table.
#'   * `imd_orig_count` and `imd_orig_proportion`, including zeros where the
#'     receptor was unobserved. Proportions use the original repertoire totals.
#'
#'   Each selected receptor has exactly one row in every input repertoire.
#'   Receptor features sharing a name with repertoire metadata receive the
#'   suffix `_receptor`. Barcode, chain-observation, and other row-level fields
#'   are excluded. The tibble is recognised by [vis()] and needs no attached
#'   repertoire metadata. Empty selections produce a zero-row summary with
#'   the same column structure.
#'
#' @seealso
#' * [vis()] for plotting the tracking abundance summary.
#' * [immundata::top_n()] for selecting abundant receptors upstream.
#' * [immundata::filter_immundata()] for filtering observations upstream.
#'
#' @examples
#' library(immundata)
#'
#' # A small bulk dataset sampled at four time points.
#' # At T3, only TTT is observed; neither AAA nor CCC is observed.
#' rows <- tibble::tibble(
#'   sample = c(rep(c("T0", "T1", "T2"), each = 3), "T3"),
#'   junction = c(
#'     "AAA", "CCC", "GGG", "AAA", "CCC", "TTT",
#'     "AAA", "GGG", "ACG", "TTT"
#'   ),
#'   v_call = c("V1", "V2", "V1", "V1", "V2", "V3", "V1", "V1", "V4", "V3"),
#'   abundance = c(8L, 2L, 1L, 1L, 9L, 4L, 2L, 3L, 5L, 5L)
#' )
#' input <- tempfile(fileext = ".tsv")
#' readr::write_tsv(rows, input)
#' idata <- read_repertoires(
#'   input,
#'   schema = c("junction", "v_call"), count_col = "abundance",
#'   repertoire_schema = "sample", output_folder = tempfile("vdj-track-"),
#'   preprocess = NULL, postprocess = NULL, rename_columns = NULL,
#'   verbose = FALSE
#' )
#' source_id <- idata$repertoires |>
#'   dplyr::filter(sample == "T0") |>
#'   dplyr::pull(imd_repertoire_id)
#'
#' # Follow every receptor observed at T0, or only its two most abundant ones.
#' from_t0 <- vdj_track(idata, from = source_id)
#'
#' @examplesIf identical(Sys.getenv("IMMUNARCH_RUN_EXAMPLES"), "true")
#' top_from_t0 <- vdj_track(idata, from = source_id, n = 2L)
#'
#' # Follow any receptor that is the most abundant in at least one repertoire.
#' top_from_each <- vdj_track(idata, n = 1L)
#'
#' # An upstream selection can also supply the receptors to follow.
#' selected <- immundata::top_n(idata, n = 1L)
#' from_selected <- vdj_track(idata, from = selected)
#'
#' # Match sequences alone, or exact sequence/gene combinations.
#' from_sequences <- vdj_track(idata, from = tibble::tibble(
#'   junction = c("AAA", "GGG")
#' ))
#' from_panel <- vdj_track(idata, from = tibble::tibble(
#'   junction = c("AAA", "CCC"), v_call = c("V1", "V2")
#' ), output = "df")
#'
#' # The summary already has one row per receptor and repertoire: 2 x 4 rows.
#' from_panel
#'
#' # Both selected receptors have explicit zero counts and proportions at T3.
#' from_panel |>
#'   dplyr::filter(sample == "T3") |>
#'   dplyr::select(sample, junction, imd_orig_count, imd_orig_proportion)
#'
#' # Select baseline receptors, then visualise their original abundances.
#' tracked <- vdj_track(idata, from = source_id, n = 2L, output = "df")
#' tracked$sample <- factor(tracked$sample, levels = c("T0", "T1", "T2", "T3"))
#' # T3 remains on each plot, with both receptor abundances at zero.
#' # Legend labels combine receptor features automatically, for example AAA V1.
#' vis(tracked, xval = "sample")
#' vis(tracked, xval = "sample", plot = "area")
#' vis(tracked, xval = "sample", plot = "line")
#'
#' # Plot original counts instead of proportions.
#' vis(tracked, xval = "sample", plot = "line", yval = "imd_orig_count")
#'
#' @export
vdj_track <- register_method_vdj(
  core = vdj_track_impl,
  family = "vdj",
  name = "track",
  result_class = TRUE
)
