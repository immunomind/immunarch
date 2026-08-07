#' @title Descriptive immune repertoire statistics
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions that extract **core descriptive statistics** from an `ImmunData` object.
#'
#' ## Available functions
#'
#' Supported methods are the following.
#'
#' @param idata An `ImmunData` object.
#' @inheritParams airr_desc_chains
#' @inheritParams airr_desc_lengths
#' @inheritParams airr_desc_genes
#' @inheritParams im_common_args
#'
#' @seealso [immundata::ImmunData]
#'
#' @examples
#' # Limit the number of threads used by the underlying DB for this session.
#' # Change this only if you know what you're doing (e.g., multi-user machines, shared CI/servers).
#' db_exec("SET threads TO 2")
#'
#' # Load data
#' \dontrun{
#' immdata <- get_test_idata() |> agg_repertoires("Therapy")
#' }
#'
#' @name airr_desc
#' @concept Key AIRR statistics
NULL


#' @keywords internal
airr_desc_chains_impl <- function(idata, locus_col = NA) {
  checkmate::assert_character(locus_col, null.ok = TRUE)

  if (is.null(idata$repertoires)) {
    cli::cli_abort("No repertoires in the input ImmunData. Run {.code agg_repertoires} first.")
  }

  if (is.na(locus_col)) {
    locus_col <- immundata::imd_schema("locus")

    if (!locus_col %in% colnames(idata$annotations)) {
      locus_col <- NULL
    }
  } else if (!is.null(locus_col)) {
    if (!locus_col %in% colnames(idata$annotations)) {
      cli::cli_alert_warning("No locus column {.code locus_col} found.")
    }
  }

  repertoire_id_col <- immundata::imd_schema("repertoire")

  by_cols <- c(repertoire_id_col, locus_col)

  chain_stats <- idata$annotations |>
    summarise(
      .by = all_of(by_cols),
      n_chains = n()
    ) |>
    collect()

  chain_stats <- idata$repertoires |>
    select(-idata$schema_repertoire) |>
    left_join(chain_stats, by = repertoire_id_col) |>
    collect()

  if (is.null(locus_col)) {
    chain_stats |> mutate(locus = NA)
  } else {
    chain_stats |> rename(locus = all_of(locus_col))
  }
}


#' @description `airr_desc_chains` --- count V(D)J *chains* per repertoire
#'   (optionally split by locus). Quickly gauges capture depth per repertoire
#'   and, when split by locus, reveals TRA/TRB/IGH balance. Use it for QC,
#'   library-size checks, and to spot locus-specific dropouts or
#'   over-representation.
#'
#' @section Visualising descriptive statistics:
#' All three `airr_desc_*()` results can be passed directly to [vis()]. With
#' `autojoin = TRUE` (the default), repertoire summaries and user metadata are
#' included in the result and can be selected by the plotting arguments.
#'
#' **Chain statistics.** `vis()` plots `n_chains` per repertoire by default.
#' The plot can be remapped to other result columns:
#'
#' * `xval` selects the column shown on the x-axis.
#' * `yval` selects the numeric column shown on the y-axis.
#' * `fill` colours and groups observations by a column.
#' * `facet` splits the plot by one column, or creates a facet grid when given
#'   two columns.
#' * `title` replaces the default plot title.
#'
#' For example, use `yval = "n_receptors"` to plot the auto-joined receptor
#' count instead of the chain count.
#'
#' **Length statistics.** `vis()` plots `seq_len` on the x-axis and `prop` on
#' the y-axis, filled by repertoire. Set `fill` to a metadata or grouping
#' column to compare its groups with box-and-violin distributions. `facet`
#' accepts one column for a wrapped layout or two columns for a facet grid;
#' `dir = "h"` or `dir = "v"` controls the wrapping direction.
#'
#' **Gene statistics.** `vis()` produces a dot plot with gene segments as rows,
#' repertoires as columns, and `n` mapped to dot size and colour. Use `row` to
#' select the gene column, `col` to select one or more repertoire, grouping, or
#' metadata columns, and `value` to select the numeric measure. `size_max_mm`
#' controls the largest dot; `row_order` and `col_order` set display order.
#'
#' @param locus_col Column in `idata$annotations` that stores the locus (e.g.
#'   `"locus"`). If `NULL` or missing, the result is not split by locus.
#'
#' @return
#'
#' ## `airr_desc_chains` Returns a tibble with columns:
#' * `repertoire_id` -- repertoire identifier
#' * `locus` -- TRA, TRB, IGH, ... (present only if `locus_col` is supplied)
#' * `n_chains` -- number of chains
#'
#' @examples
#' #
#' # airr_desc_chains
#' #
#'
#' \dontrun{
#' chain_stats <- airr_desc_chains(immdata)
#'
#' # Default: chain counts per repertoire
#' vis(chain_stats)
#'
#' # Plot another numeric result column
#' vis(
#'   chain_stats,
#'   yval = "n_receptors",
#'   title = "No. receptors per sample"
#' )
#'
#' # Use auto-joined metadata and locus information for comparison
#' vis(chain_stats, xval = "Therapy", fill = "locus")
#' vis(chain_stats, xval = "Therapy", facet = "locus")
#' }
#'
#' @rdname airr_desc
#' @concept Key AIRR statistics
#' @export
airr_desc_chains <- register_immunarch_method(
  core = airr_desc_chains_impl,
  family = "airr_desc",
  name = "chains"
)


#' @keywords internal
airr_desc_lengths_impl <- function(
  idata,
  seq_col = "cdr3_aa",
  by = NA_character_
) {
  checkmate::assert_string(seq_col)
  by_cols <- resolve_argument_by(idata, by)

  repertoire_col <- immundata::imd_schema("repertoire")
  output_group_cols <- c(repertoire_col, by_cols, "seq_len")
  proportion_group_cols <- c(repertoire_col, by_cols)

  idata$annotations |>
    dplyr::select(dplyr::all_of(c(
      repertoire_col,
      by_cols,
      seq_col
    ))) |>
    dplyr::mutate(seq_len = dd$length(!!rlang::sym(seq_col))) |>
    dplyr::summarise(
      n = dplyr::n(),
      .by = dplyr::all_of(output_group_cols)
    ) |>
    dplyr::mutate(
      prop = n / sum(n, na.rm = TRUE),
      pct = 100 * prop,
      .by = dplyr::all_of(proportion_group_cols)
    ) |>
    dplyr::select(dplyr::all_of(c(
      repertoire_col, by_cols, "seq_len", "n", "prop", "pct"
    ))) |>
    dplyr::arrange(!!!rlang::syms(c(repertoire_col, by_cols, "seq_len"))) |>
    collect()
}


#' @description `airr_desc_lengths` --- count the number of sequence lengths per
#' repertoire, optionally split by annotation columns. Summarizes the CDR3
#' length distribution, a sensitive QC fingerprint of repertoire prep and
#' selection. Helpful for detecting primer/UMI biases, comparing cohorts, and
#' deriving length-based features for models.
#'
#' @param seq_col Name of the column containing sequences.
#' @param by Grouping columns from `idata$annotations`. The default `NA`
#'   automatically groups by the canonical `locus` column when it is available.
#'   Supply `NULL` to pool loci, or a character vector to group by and return
#'   those columns explicitly.
#'
#' @return
#'
#' ## `airr_desc_lengths` Returns a tibble with columns:
#' * `imd_repertoire_id` -- internal repertoire identifier
#' * grouping columns requested through `by` (for example, `locus`)
#' * `seq_len` -- lengths of sequences
#' * `n` -- number of sequence rows
#' * `prop`, `pct` -- proportion and percentage within each repertoire and
#'   combination of grouping columns
#'
#' @examples
#' #
#' # airr_desc_lengths
#' #
#'
#' \dontrun{
#' length_stats <- airr_desc_lengths(immdata)
#'
#' # Default: CDR3 length proportions by repertoire
#' vis(length_stats)
#'
#' # Compare distributions between metadata groups and split them by locus
#' vis(length_stats, fill = "Therapy", facet = "locus")
#'
#' airr_desc_lengths(immdata, by = "locus")
#' airr_desc_lengths(immdata, by = NULL)
#' }
#'
#' @rdname airr_desc
#' @concept Key AIRR statistics
#' @export
airr_desc_lengths <- register_immunarch_method(
  core = airr_desc_lengths_impl,
  family = "airr_desc",
  name = "lengths",
  required = "seq_col"
)


#' @keywords internal
airr_desc_genes_impl <- function(
  idata,
  gene_col = "v_call",
  level = c("receptor", "barcode"),
  by = NA_character_
) {
  checkmate::assert_string(gene_col)
  level <- match.arg(level)
  by_cols <- resolve_argument_by(idata, by)

  repertoire_col <- immundata::imd_schema("repertoire")
  receptor_col <- immundata::imd_schema("receptor")
  count_col <- immundata::imd_schema("count")
  observation_cols <- c(repertoire_col, receptor_col, by_cols, gene_col)
  output_group_cols <- c(repertoire_col, by_cols, gene_col)

  receptors <- idata$annotations |>
    dplyr::select(dplyr::all_of(c(
      repertoire_col,
      receptor_col,
      if (level == "barcode") count_col,
      by_cols,
      gene_col
    )))

  if (level == "receptor") {
    genes <- receptors |>
      dplyr::distinct() |>
      dplyr::summarise(
        n = dplyr::n(),
        .by = dplyr::all_of(output_group_cols)
      )
  } else {
    genes <- receptors |>
      dplyr::summarise(
        .imd_count = dplyr::first(!!rlang::sym(count_col)),
        .by = dplyr::all_of(observation_cols)
      ) |>
      dplyr::summarise(
        n = sum(.data$.imd_count, na.rm = TRUE),
        .by = dplyr::all_of(output_group_cols)
      )
  }

  genes |>
    dplyr::select(dplyr::all_of(c(gene_col, repertoire_col, by_cols, "n"))) |>
    dplyr::arrange(
      !!!rlang::syms(c(repertoire_col, by_cols)),
      dplyr::desc(.data$n)
    ) |>
    collect()
}

#' @description `airr_desc_genes` - count V(D)J gene segments per repertoire,
#'   optionally split by locus and using either receptor counts or barcode/UMI
#'   counts as the measure. Profiles V/D/J gene usage to characterize repertoire
#'   composition and germline biases, with optional locus split. Useful for
#'   cohort comparisons, flagging clonal expansions, and producing ML-ready
#'   features for repertoire-level ML tasks.
#'
#' @param gene_col A single column name in `idata$annotations` with gene segment
#'   calls (e.g., `"v_call"`, `"d_call"`, `"j_call"`, `"c_call"`). Default is
#'   `"v_call"`.
#' @param level One of `"receptor"` or `"barcode"`. If `"receptor"` (default),
#'   the function counts **unique receptors** (one per receptor ID) that carry a
#'   given gene segment. If `"barcode"`, the function **sums counts** (e.g.,
#'   cells/UMIs) per gene segment using the column defined by
#'   `immundata::imd_schema("count")`.
#' @param by Grouping columns from `idata$annotations`. The default `NA`
#'   automatically groups by the canonical `locus` column when it is available.
#'   Supply `NULL` to pool loci, or a character vector to group by and return
#'   those columns explicitly.
#'
#' @return
#'
#' ## `airr_desc_genes` A tibble with columns:
#' * `<gene_col>` - the gene segment value (e.g., V gene)
#' * `imd_repertoire_id` - internal repertoire identifier
#' * grouping columns requested through `by` (for example, `locus`)
#' * `n` - the measure:
#'   - if `level = "receptor"`: number of receptors carrying the gene segment
#'     within the grouping
#'   - if `level = "barcode"`: sum of counts across receptors for the segment
#'
#' @examples
#' #
#' # airr_desc_genes
#' #
#'
#' \dontrun{
#' # V gene usage by receptor count
#' gene_stats <- airr_desc_genes(
#'   immdata,
#'   gene_col = "v_call",
#'   level = "receptor"
#' )
#'
#' # Default: genes by repertoire, with counts mapped to dot size and colour
#' vis(gene_stats)
#'
#' # Compare gene usage across metadata groups and loci
#' vis(
#'   gene_stats,
#'   row = "v_call",
#'   col = c("Therapy", "locus"),
#'   value = "n",
#'   size_max_mm = 5
#' )
#'
#' # V gene usage by summed cell/UMI counts (if a count column is present)
#' airr_desc_genes(immdata, gene_col = "v_call", level = "barcode")
#'
#' # Split by locus (TRA/TRB/... if locus column exists)
#' airr_desc_genes(immdata, gene_col = "v_call", level = "receptor", by = "locus")
#'
#' # Pool loci explicitly
#' airr_desc_genes(immdata, gene_col = "v_call", by = NULL)
#' }
#'
#' @rdname airr_desc
#' @concept Key AIRR statistics
#' @export
airr_desc_genes <- register_immunarch_method(
  core = airr_desc_genes_impl,
  family = "airr_desc",
  name = "genes",
  required = "gene_col"
)
