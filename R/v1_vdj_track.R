#' @keywords internal
vdj_track_impl <- function(idata, from = NULL, n = NULL) {
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
                  imd_orig_count = imd_count, imd_orig_proportion = imd_proportion) |>
    dplyr::distinct()
  immundata::annotate_immundata(
    idata, annotations = abundances, by = stats::setNames(keys, keys),
    conflicts = "replace", remove_limit = TRUE
  ) |>
    immundata::filter_immundata(!is.na(imd_orig_count))

}


#' @title Follow selected receptors across repertoires
#'
#' @description
#' Use `vdj_track()` to follow receptors across samples, time points, tissues,
#' or cell populations. For example, select abundant receptors before treatment
#' and find their abundances in later samples. You can also follow a set of
#' sequences identified in another experiment.
#'
#' Select receptors with `from`, or use `n` to select the most abundant
#' receptors in each source repertoire. The result contains their observations
#' throughout the target `idata`, including repertoires where they are rare.
#' The original object is not changed.
#'
#' @param idata An [immundata::ImmunData] object with defined repertoires.
#'   These are the repertoires through which receptors will be followed.
#' @param from `NULL`, a numeric vector of repertoire IDs, an
#'   [immundata::ImmunData] object, or a data frame. `NULL`, the default, uses
#'   every target repertoire as a source. Numeric values select source
#'   `imd_repertoire_id` values in `idata`; they are IDs, not row positions.
#'   An `ImmunData` source supplies its observed receptor identities. A data
#'   frame supplies exact column combinations to match, such as `junction`
#'   alone or `junction` together with `v_call`. Tibbles and duckplyr tables
#'   are also accepted.
#' @param n A non-negative whole number or `NULL`. Select up to `n` receptors
#'   with the highest counts in each source repertoire, then follow their
#'   combined set of identities. Equal counts are resolved by receptor ID, as
#'   in [immundata::top_n()]. `NULL`, the default, selects all source receptors.
#'   This argument applies only when `from` is `NULL` or contains repertoire
#'   IDs. Use upstream selection for an `ImmunData` source or a query table.
#' @param output A character scalar. `"immundata"`, the default, returns a new
#'   `ImmunData` object. `"data.frame"` collects the same observations into a
#'   tibble. Both outputs contain the original target abundances.
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
#' paired receptor, are matched. Numeric receptor IDs from another object are
#' not treated as target IDs. Source abundances and metadata are not imported.
#' A source object does not need repertoire definitions.
#'
#' Query-table column names must exist in the target observations. Matching is
#' exact and case-sensitive, with no sequence or gene-call normalization.
#' A query with fewer features than the receptor schema can select several
#' different receptors; their identities remain separate. Missing query values
#' are not allowed. Duplicate query rows do not duplicate the result. Unmatched
#' queries, empty query tables, and `n = 0` produce an empty result.
#'
#' @section Original abundances:
#' `imd_orig_count` and `imd_orig_proportion` store the target counts and
#' proportions before this call filters the data. Use these columns to compare
#' abundances across the original samples. The ordinary `imd_count` and
#' `imd_proportion` describe the filtered object, following normal immundata
#' behavior. Another tracking call replaces the original-abundance columns
#' using that call's input population.
#'
#' Counts and proportions can repeat across cell and chain observations. Use
#' `distinct()` to obtain one value per receptor and repertoire before adding
#' abundances or preparing a trajectory table.
#'
#' The result contains observed chains only. No artificial zero-count rows are
#' added for receptors not observed in a sample. A receptor not observed in a
#' sample is not necessarily biologically absent.
#'
#' @return A new [immundata::ImmunData] object or a collected tibble containing
#'   all target observations of the selected receptors, with `imd_orig_count`
#'   and `imd_orig_proportion`.
#'
#' @seealso [immundata::top_n()], [immundata::filter_immundata()]
#'
#' @examples
#' library(immundata)
#'
#' # A small bulk dataset sampled at three time points.
#' rows <- tibble::tibble(
#'   sample = rep(c("T0", "T1", "T2"), each = 3),
#'   junction = c("AAA", "CCC", "GGG", "AAA", "CCC", "TTT",
#'                "AAA", "GGG", "ACG"),
#'   v_call = c("V1", "V2", "V1", "V1", "V2", "V3", "V1", "V1", "V4"),
#'   abundance = c(8L, 2L, 1L, 1L, 9L, 4L, 2L, 3L, 5L)
#' )
#' input <- tempfile(fileext = ".tsv")
#' readr::write_tsv(rows, input)
#' idata <- read_repertoires(
#'   input, schema = c("junction", "v_call"), count_col = "abundance",
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
#' ), output = "data.frame")
#' from_panel |>
#'   dplyr::distinct(sample, junction, imd_orig_count, imd_orig_proportion)
#'
#' @export
vdj_track <- register_method_vdj(
  core = vdj_track_impl,
  family = "vdj",
  name = "track"
)
