#' @title Compute key immune repertoire statistics
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
#' @inheritParams airr_stats_chains
#' @inheritParams airr_stats_lengths
#' @inheritParams airr_stats_genes
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
#' @name airr_stats
#' @concept Key AIRR statistics
NULL


#' @keywords internal
airr_stats_chains_impl <- function(idata, locus_col = NA) {
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
    chain_stats |> rename(locus = locus_col)
  }
}


#' @description `airr_stats_chains` --- count V(D)J *chains* per repertoire
#'   (optionally split by locus). Quickly gauges capture depth per repertoire
#'   and, when split by locus, reveals TRA/TRB/IGH balance. Use it for QC,
#'   library-size checks, and to spot locus-specific dropouts or
#'   over-representation.
#'
#' @param locus_col Column in `idata$annotations` that stores the locus (e.g.
#'   `"locus"`). If `NULL` or missing, the result is not split by locus.
#'
#' @return
#'
#' ## `airr_stats_chains` Returns a tibble with columns:
#' * `repertoire_id` -- repertoire identifier
#' * `locus` -- TRA, TRB, IGH, ... (present only if `locus_col` is supplied)
#' * `n_chains` -- number of chains
#'
#' @examples
#' #
#' # airr_stats_chains
#' #
#'
#' \dontrun{
#' airr_stats_chains(immdata)
#' }
#'
#' @rdname airr_stats
#' @concept Key AIRR statistics
#' @export
airr_stats_chains <- register_immunarch_method(
  core = airr_stats_chains_impl,
  family = "airr_stats",
  name = "chains"
)


#' @keywords internal
airr_stats_lengths_impl <- function(idata, seq_col = "cdr3_aa") {
  idata$annotations |>
    dplyr::select(dplyr::all_of(c(immundata::imd_schema("repertoire"), seq_col))) |>
    dplyr::mutate(seq_len = dd$length(!!rlang::sym(seq_col))) |>
    dplyr::summarise(
      n = dplyr::n(),
      .by = dplyr::all_of(c(immundata::imd_schema("repertoire"), "seq_len"))
    ) |>
    collect() |>
    dplyr::mutate(
      prop = n / sum(n, na.rm = TRUE), # proportion within repertoire
      pct = 100 * prop,
      .by = immundata::imd_schema("repertoire")
    )
}


#' @description `airr_stats_lengths` --- count the number of sequence lengths
#' per repertoire. Summarizes the CDR3 length distribution, a sensitive QC
#' fingerprint of repertoire prep and selection. Helpful for detecting
#' primer/UMI biases, comparing cohorts, and deriving length-based features for
#' models.
#'
#' @param seq_col Character vector with names of the columns containing
#'   sequences.
#'
#' @return
#'
#' ## `airr_stats_lengths` Returns a tibble with columns:
#' * `repertoire_id` -- repertoire identifier
#' * `seq_len` -- lengths of sequences
#' * `n` -- number of receptors
#'
#' @examples
#' #
#' # airr_stats_lengths
#' #
#'
#' \dontrun{
#' airr_stats_lengths(immdata)
#' }
#'
#' @rdname airr_stats
#' @concept Key AIRR statistics
#' @export
airr_stats_lengths <- register_immunarch_method(
  core = airr_stats_lengths_impl,
  family = "airr_stats",
  name = "lengths",
  required = "seq_col"
)


#' @keywords internal
airr_stats_genes_impl <- function(idata, gene_col = "v_call", level = c("receptor", "barcode"), by = c(NA, "locus")) {
  checkmate::assert_logical(gene_col %in% colnames(idata$annotations))
  level <- match.arg(level)
  by <- match.arg(by)
  if (!is.na(by)) {
    checkmate::assert_logical(by %in% colnames(idata$annotations))
  }

  receptors <- idata$annotations |>
    distinct(!!rlang::sym(immundata::imd_schema("receptor")), !!rlang::sym(immundata::imd_schema("repertoire")), .keep_all = TRUE)

  if (level == "receptor") {
    genes <- receptors |>
      summarise(.by = all_of(c(gene_col, immundata::imd_schema("repertoire"))), n = n())
  } else {
    genes <- receptors |>
      summarise(.by = all_of(c(gene_col, immundata::imd_schema("repertoire"))), n = sum(!!rlang::sym(immundata::imd_schema("count"))))
  }

  genes |>
    arrange(!!rlang::sym(immundata::imd_schema("repertoire")), desc(n)) |>
    collect()
}

#' @description `airr_stats_genes` - count V(D)J gene segments per repertoire,
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
#' @param by Either `NULL` (no split) or `"locus"`. When `"locus"`, the result
#'   is further split by the locus column if present (as given by
#'   `immundata::imd_schema("locus")`); otherwise a warning is emitted and the
#'   split is ignored.
#'
#' @return
#'
#' ## `airr_stats_genes` A tibble with columns:
#' * `repertoire_id` - repertoire identifier
#' * *(optional)* `locus` - TRA, TRB, IGH, ... (present only when `by = "locus"`
#' and the locus column exists)
#' * `<gene_col>` - the gene segment value (e.g., V gene)
#' * `n` - the measure:
#'   - if `level = "receptor"`: number of receptors carrying the gene segment
#'   - if `level = "barcode"`: sum of counts across receptors for the segment
#'
#' @examples
#' #
#' # airr_stats_genes
#' #
#'
#' \dontrun{
#' # V gene usage by receptor count
#' airr_stats_genes(immdata, gene_col = "v_call", level = "receptor")
#'
#' # V gene usage by summed cell/UMI counts (if a count column is present)
#' airr_stats_genes(immdata, gene_col = "v_call", level = "barcode")
#'
#' # Split by locus (TRA/TRB/... if locus column exists)
#' airr_stats_genes(immdata, gene_col = "v_call", level = "receptor", by = "locus")
#' }
#'
#' @rdname airr_stats
#' @concept Key AIRR statistics
#' @export
airr_stats_genes <- register_immunarch_method(
  core = airr_stats_genes_impl,
  family = "airr_stats",
  name = "genes",
  required = "gene_col"
)
