#' @title Annotate immunoglobulin somatic hypermutation
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' Calculate overall and AIRR-region somatic-hypermutation counts and
#' frequencies from aligned observed and D-masked germline sequences.
#'
#' The input must contain the columns selected by `sequence_col` and
#' `germline_col`, `v_sequence_start`, `v_alignment_start`, and the canonical
#' AIRR `*_start` and `*_end` columns for FWR1, CDR1, FWR2, CDR2, FWR3, CDR3,
#' and FWR4. The default D-masked germline column is normally created with
#' [annotate_ig_germlines()].
#'
#' @param idata [immundata::ImmunData]. IgBLAST AIRR data to annotate.
#' @param sequence_col Character scalar. Aligned observed-sequence column.
#' @param germline_col Character scalar. Aligned D-masked germline column.
#'   `sequence_col` and `germline_col` may select alternative alignment columns,
#'   but those columns must use the same alignment and coordinate system as the
#'   canonical AIRR region-coordinate columns.
#' @param by Character vector, `"<auto>"`, or `NULL`. Columns used to calculate
#'   group-level `mu_freq` summaries. The default resolves to `family_id` when
#'   available, otherwise to `subject_id` plus `clone_id`; `locus` is added to
#'   either automatic grouping when available. `NULL` skips grouped summaries.
#'
#' @details
#' AIRR region coordinates are interpreted as 1-based closed intervals on the
#' unaligned query sequence. `v_alignment_start` identifies the aligned
#' position corresponding to `v_sequence_start`. While mapping query positions,
#' `N` consumes a query coordinate, whereas `.` and `-` do not.
#'
#' A position is comparable only when both aligned characters are unambiguous
#' A, C, G, or T nucleotides. Thus `N`, `.`, and `-` positions contribute to
#' neither mutation counts nor frequency denominators. `mu_count` is the sum of
#' classifiable replacement and silent mutations across all available regions;
#' `mu_freq` divides that count by the number of comparable positions in those
#' regions. Overall and combined CDR/FWR replacement and silent frequencies use
#' the sum of comparable positions across their selected regions, not aggregated
#' regional frequencies. Individual-region frequencies use the comparable
#' positions in that region as their denominator.
#'
#' Replacement and silent consequences are classified independently for each
#' mutation against the germline codon. Stop-producing substitutions are not
#' counted as either class. If one boundary of a region is missing, all four
#' outputs for that region are missing while other regions and overall outputs
#' remain available. A valid region with no comparable positions has zero
#' counts and missing frequencies. Combined CDR and framework metrics use the
#' available regions in their respective sets. If no region in a set has both
#' boundaries, all combined outputs for that set are missing.
#'
#' When grouping is enabled, `mu_freq_median` and `mu_freq_mean` summarize
#' non-missing `mu_freq` values without weighting, and `mu_freq_n` counts their
#' contributors. Rows missing any grouping key receive missing median and mean
#' values and a zero count. Explicit `by` columns are used exactly as supplied;
#' include `locus` explicitly to keep heavy and light chains separate. Use
#' `dplyr::distinct()` on the grouping columns and summaries when a group-level
#' rather than row-level result is needed.
#'
#' @return [immundata::ImmunData]. The input rows with combined `mu_count` and
#'   `mu_freq` columns plus replacement and silent counts and frequencies
#'   overall, across all CDRs, across all framework regions, and for each
#'   individual region. When grouping is enabled, the result also includes
#'   `mu_freq_median`, `mu_freq_mean`, and `mu_freq_n` on every row.
#'
#' @name annotate_ig_shm
NULL


#' @keywords internal
annotate_ig_shm_impl <- function(
  idata,
  sequence_col = "sequence_alignment",
  germline_col = "germline_alignment_d_mask",
  by = "<auto>"
) {
  regions <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
  annotation_cols <- names(idata$annotations)
  locus_col <- immundata::imd_schema("locus")

  checkmate::assert_character(
    by,
    min.len = 1L,
    any.missing = FALSE,
    unique = TRUE,
    null.ok = TRUE
  )

  if (identical(by, "<auto>")) {
    if ("family_id" %in% annotation_cols) {
      resolved_by <- c("family_id", intersect(locus_col, annotation_cols))
    } else if (all(c("subject_id", "clone_id") %in% annotation_cols)) {
      resolved_by <- c(
        "subject_id", "clone_id", intersect(locus_col, annotation_cols)
      )
    } else {
      resolved_by <- NULL
    }

    if (is.null(resolved_by)) {
      cli::cli_inform(
        "{.arg by = \"<auto>\"} found no supported grouping columns; grouped SHM statistics were skipped."
      )
    } else {
      cli::cli_inform(
        "{.arg by = \"<auto>\"} resolved to {.field {resolved_by}}."
      )
    }
  } else {
    if ("<auto>" %in% by) {
      cli::cli_abort(
        "{.code \"<auto>\"} must be supplied alone through {.arg by}."
      )
    }
    resolved_by <- by
  }

  missing_by <- setdiff(resolved_by, annotation_cols)
  if (length(missing_by)) {
    cli::cli_abort(
      "Column(s) {.field {missing_by}} selected by {.arg by} are missing from {.field idata$annotations}."
    )
  }

  coordinate_cols <- c(
    "v_sequence_start", "v_alignment_start",
    unlist(lapply(regions, function(region) {
      paste0(region, c("_start", "_end"))
    }), use.names = FALSE)
  )

  missing_cols <- setdiff(coordinate_cols, names(idata$annotations))
  if (length(missing_cols)) {
    cli::cli_abort(
      "Passed column name(s) [{.code {missing_cols}}] is not in the input ImmunData."
    )
  }

  con <- idata$annotations |>
    duckplyr::as_tbl() |>
    dbplyr::remote_con()

  genetic_code <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )

  codons <- names(genetic_code)
  amino_acids <- unname(genetic_code)
  genetic_code_values <- glue::glue_sql_collapse(
    glue::glue_sql("({codons}, {amino_acids})", .con = con),
    sep = ", "
  )

  region_starts <- paste0(regions, "_start_value")
  region_ends <- paste0(regions, "_end_value")
  region_values <- glue::glue_sql_collapse(
    glue::glue_sql(
      "({regions}, cast({`region_starts`} AS BIGINT),
        cast({`region_ends`} AS BIGINT))",
      .con = con
    ),
    sep = ", "
  )

  regional_fields <- lapply(regions, function(region) {
    output_cols <- paste0(
      c("mu_count_", "mu_count_", "mu_freq_", "mu_freq_"),
      region,
      c("_r", "_s", "_r", "_s")
    )
    values <- DBI::SQL(c("r_count", "s_count", "r_freq", "s_freq"))

    glue::glue_sql(
      "{`output_cols`} := max({values}) FILTER (region = {region})",
      .con = con
    )
  })
  regional_fields <- DBI::SQL(unlist(regional_fields, use.names = FALSE))

  result_fields <- c(
    glue::glue_sql(
      "mu_count := CASE
        WHEN count(denominator) = 0 THEN NULL::BIGINT
        ELSE sum(r_count + s_count)::BIGINT
      END",
      .con = con
    ),
    glue::glue_sql(
      "mu_freq := CASE
        WHEN coalesce(sum(denominator), 0) = 0 THEN NULL::DOUBLE
        ELSE sum(r_count + s_count)::DOUBLE / sum(denominator)
      END",
      .con = con
    ),
    ig_shm_aggregate_fields("", con = con),
    ig_shm_aggregate_fields("_cdr", paste0("cdr", 1:3), con),
    ig_shm_aggregate_fields("_fwr", paste0("fwr", 1:4), con),
    regional_fields
  )
  result_fields <- glue::glue_sql_collapse(result_fields, sep = ", ")

  macro_parameters <- c(
    "observed_alignment", "germline_alignment",
    "v_sequence_start_value", "v_alignment_start_value",
    unlist(lapply(regions, function(region) {
      paste0(region, c("_start_value", "_end_value"))
    }), use.names = FALSE)
  )
  macro_parameters <- glue::glue_sql("{`macro_parameters`*}", .con = con)
  macro_sql <- glue::glue_sql(
    "
    CREATE OR REPLACE TEMP MACRO immunarch_annotate_ig_shm(
      {{macro_parameters}}
    ) AS (
      CASE
        WHEN observed_alignment IS NULL OR germline_alignment IS NULL
          THEN NULL
        ELSE (
          WITH bases AS (
            SELECT
              pos,
              substr(observed_alignment, pos, 1) AS obs_base,
              substr(germline_alignment, pos, 1) AS germ_base
            FROM range(
              cast(v_alignment_start_value AS BIGINT),
              length(observed_alignment) + 1
            ) AS p(pos)
          ), positions AS (
            SELECT
              *,
              CASE WHEN obs_base NOT IN ('-', '.') THEN
                cast(v_sequence_start_value AS BIGINT) +
                sum(
                  CASE WHEN obs_base NOT IN ('-', '.') THEN 1 ELSE 0 END
                ) OVER (
                  ORDER BY pos
                  ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
                ) - 1
              END AS query_pos
            FROM bases
          ), regions(region, region_start, region_end) AS (
            VALUES {{region_values}}
          ), genetic_code(codon, amino_acid) AS (
            VALUES {{genetic_code_values}}
          ), mapped AS (
            SELECT
              r.*, p.pos, p.obs_base, p.germ_base, p.query_pos,
              ((p.query_pos - r.region_start) % 3 + 1)::INTEGER AS codon_offset
            FROM regions r
            LEFT JOIN positions p
              ON p.query_pos BETWEEN r.region_start AND r.region_end
          ), codon_strings AS (
            SELECT
              region,
              floor((query_pos - region_start) / 3) AS codon_index,
              string_agg(obs_base, '' ORDER BY pos) AS observed_codon,
              string_agg(germ_base, '' ORDER BY pos) AS germ_codon,
              NOT regexp_matches(
                substr(
                  observed_alignment, min(pos), max(pos) - min(pos) + 1
                ),
                '[-.]'
              ) AND NOT regexp_matches(
                substr(
                  germline_alignment, min(pos), max(pos) - min(pos) + 1
                ),
                '[-.]'
              ) AS gap_free_codon
            FROM mapped
            WHERE query_pos IS NOT NULL
            GROUP BY region, codon_index
          ), codons AS (
            SELECT
              m.*, c.observed_codon, c.germ_codon, c.gap_free_codon
            FROM mapped m
            LEFT JOIN codon_strings c
              ON m.region = c.region
              AND floor((m.query_pos - m.region_start) / 3) = c.codon_index
          ), substitutions AS (
            SELECT
              *,
              regexp_full_match(obs_base, '[ACGT]')
                AND regexp_full_match(germ_base, '[ACGT]') AS comparable,
              substr(germ_codon, 1, codon_offset - 1)
                || obs_base
                || substr(germ_codon, codon_offset + 1) AS mutant_codon
            FROM codons
          ), translated AS (
            SELECT
              s.*,
              germ_code.amino_acid AS germ_amino_acid,
              mutant_code.amino_acid AS mutant_amino_acid
            FROM substitutions s
            LEFT JOIN genetic_code germ_code
              ON s.germ_codon = germ_code.codon
            LEFT JOIN genetic_code mutant_code
              ON s.mutant_codon = mutant_code.codon
          ), classified AS (
            SELECT
              *,
              CASE
                WHEN obs_base <> germ_base
                  AND gap_free_codon
                  AND length(observed_codon) = 3
                  AND regexp_full_match(observed_codon, '[ACGT]{3}')
                  AND length(germ_codon) = 3
                  AND regexp_full_match(germ_codon, '[ACGT]{3}')
                  AND mutant_amino_acid <> '*'
                THEN CASE
                  WHEN mutant_amino_acid = germ_amino_acid THEN 'S'
                  ELSE 'R'
                END
              END AS mutation_class
            FROM translated
          ), region_stats AS (
            SELECT
              region, region_start, region_end,
              CASE
                WHEN region_start IS NULL OR region_end IS NULL THEN NULL
                ELSE count(*) FILTER (WHERE comparable)::BIGINT
              END AS denominator,
              CASE
                WHEN region_start IS NULL OR region_end IS NULL THEN NULL
                ELSE count(*) FILTER (WHERE mutation_class = 'R')::BIGINT
              END AS r_count,
              CASE
                WHEN region_start IS NULL OR region_end IS NULL THEN NULL
                ELSE count(*) FILTER (WHERE mutation_class = 'S')::BIGINT
              END AS s_count
            FROM classified
            GROUP BY region, region_start, region_end
          ), frequencies AS (
            SELECT
              *,
              CASE WHEN denominator = 0 THEN NULL::DOUBLE
                ELSE r_count::DOUBLE / denominator
              END AS r_freq,
              CASE WHEN denominator = 0 THEN NULL::DOUBLE
                ELSE s_count::DOUBLE / denominator
              END AS s_freq
            FROM region_stats
          )
          SELECT struct_pack({{result_fields}})
          FROM frequencies
        )
      END
    )",
    .open = "{{",
    .close = "}}",
    .con = con
  )

  shm_output_cols <- c(
    "mu_count", "mu_freq",
    paste0(
      c("mu_count", "mu_count", "mu_freq", "mu_freq"),
      rep(c("", "_cdr", "_fwr"), each = 4L),
      rep(c("_r", "_s", "_r", "_s"), 3L)
    ),
    unlist(lapply(regions, function(region) {
      paste0(
        c("mu_count_", "mu_count_", "mu_freq_", "mu_freq_"), region,
        c("_r", "_s", "_r", "_s")
      )
    }), use.names = FALSE)
  )

  group_output_cols <- c("mu_freq_median", "mu_freq_mean", "mu_freq_n")

  output_cols <- c(
    shm_output_cols,
    if (length(resolved_by)) group_output_cols
  )

  DBI::dbExecute(con, macro_sql)

  chain_col <- immundata::imd_schema("chain")
  input_cols <- unique(c(
    chain_col, sequence_col, germline_col, coordinate_cols, resolved_by
  ))

  annotations <- idata$annotations |>
    dplyr::select(dplyr::all_of(input_cols)) |>
    dplyr::mutate(.ig_shm = dd$immunarch_annotate_ig_shm(
      !!rlang::sym(sequence_col),
      !!rlang::sym(germline_col),
      v_sequence_start,
      v_alignment_start,
      fwr1_start, fwr1_end,
      cdr1_start, cdr1_end,
      fwr2_start, fwr2_end,
      cdr2_start, cdr2_end,
      fwr3_start, fwr3_end,
      cdr3_start, cdr3_end,
      fwr4_start, fwr4_end
    ))

  for (column in shm_output_cols) {
    annotations <- annotations |>
      dplyr::mutate(
        !!column := dd$struct_extract(!!rlang::sym(".ig_shm"), !!column)
      )
  }

  annotations <- annotations |>
    dplyr::select(-.ig_shm)

  if (length(resolved_by)) {
    complete_group <- Reduce(function(condition, column) {
      rlang::expr(!!condition & !is.na(!!rlang::sym(column)))
    }, resolved_by, init = TRUE)

    annotations <- annotations |>
      dplyr::mutate(.ig_complete_group = !!complete_group) |>
      dplyr::mutate(
        mu_freq_median = dplyr::if_else(
          .ig_complete_group,
          stats::median(mu_freq, na.rm = TRUE),
          NA_real_
        ),
        mu_freq_mean = dplyr::if_else(
          .ig_complete_group,
          mean(mu_freq, na.rm = TRUE),
          NA_real_
        ),
        mu_freq_n = dplyr::if_else(
          .ig_complete_group,
          sum(!is.na(mu_freq), na.rm = TRUE),
          0L
        ),
        .by = dplyr::all_of(resolved_by)
      ) |>
      dplyr::select(-.ig_complete_group)
  }

  shm_annotations <- annotations |>
    dplyr::select(dplyr::all_of(c(chain_col, output_cols)))

  immundata::annotate_chains(
    idata,
    annotations = shm_annotations,
    annot_col = chain_col,
    conflicts = "replace"
  )
}


ig_shm_aggregate_fields <- function(label, selected_regions = NULL, con) {
  filter <- if (is.null(selected_regions)) {
    DBI::SQL("")
  } else {
    glue::glue_sql(
      " FILTER (WHERE region IN ({selected_regions*}))",
      .con = con
    )
  }
  available <- glue::glue_sql("count(denominator){filter}", .con = con)
  denominator <- glue::glue_sql("sum(denominator){filter}", .con = con)
  r_count <- glue::glue_sql("sum(r_count){filter}", .con = con)
  s_count <- glue::glue_sql("sum(s_count){filter}", .con = con)
  output_cols <- paste0(
    c("mu_count", "mu_count", "mu_freq", "mu_freq"),
    label,
    c("_r", "_s", "_r", "_s")
  )
  r_count_col <- output_cols[[1L]]
  s_count_col <- output_cols[[2L]]
  r_freq_col <- output_cols[[3L]]
  s_freq_col <- output_cols[[4L]]

  fields <- list(
    glue::glue_sql(
      "{`r_count_col`} := CASE
        WHEN {available} = 0 THEN NULL::BIGINT
        ELSE ({r_count})::BIGINT
      END",
      .con = con
    ),
    glue::glue_sql(
      "{`s_count_col`} := CASE
        WHEN {available} = 0 THEN NULL::BIGINT
        ELSE ({s_count})::BIGINT
      END",
      .con = con
    ),
    glue::glue_sql(
      "{`r_freq_col`} := CASE
        WHEN coalesce({denominator}, 0) = 0 THEN NULL::DOUBLE
        ELSE ({r_count})::DOUBLE / {denominator}
      END",
      .con = con
    ),
    glue::glue_sql(
      "{`s_freq_col`} := CASE
        WHEN coalesce({denominator}, 0) = 0 THEN NULL::DOUBLE
        ELSE ({s_count})::DOUBLE / {denominator}
      END",
      .con = con
    )
  )

  DBI::SQL(unlist(fields, use.names = FALSE))
}


#' @export
annotate_ig_shm <- register_method_annotate(
  core = annotate_ig_shm_impl,
  family = "annotate_ig",
  name = "shm",
  required_cols = c("sequence_col", "germline_col"),
  need_repertoires = FALSE
)
