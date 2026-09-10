#' @title Annotate immunoglobulin somatic hypermutation
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' Calculate overall and AIRR-region somatic-hypermutation counts and
#' frequencies from aligned observed and reconstructed germline sequences.
#'
#' The input must contain the columns selected by `sequence_col` and
#' `germline_col`, `v_sequence_start`, `v_alignment_start`, and the canonical
#' AIRR `*_start` and `*_end` columns for FWR1, CDR1, FWR2, CDR2, FWR3, CDR3,
#' and FWR4. The default germline column is normally created with
#' [annotate_ig_germlines()].
#'
#' @param idata [immundata::ImmunData]. IgBLAST AIRR data to annotate.
#' @param sequence_col Character scalar. Aligned observed-sequence column.
#' @param germline_col Character scalar. Aligned reconstructed germline column.
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
  germline_col = "germline_alignment",
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

  chain_col <- immundata::imd_schema("chain")
  input_cols <- unique(c(
    chain_col, sequence_col, germline_col, coordinate_cols, resolved_by
  ))

  source_annotations <- idata$annotations |>
    dplyr::select(dplyr::all_of(input_cols))
  source_tbl <- duckplyr::as_tbl(source_annotations)
  source_sql <- dbplyr::sql_render(source_tbl)
  quoted_coordinates <- DBI::dbQuoteIdentifier(con, coordinate_cols)
  sql_input <- paste0(
    "SELECT ", DBI::dbQuoteIdentifier(con, chain_col), " AS imd_chain_id, ",
    DBI::dbQuoteIdentifier(con, sequence_col), " AS observed_alignment, ",
    DBI::dbQuoteIdentifier(con, germline_col), " AS germline_alignment, ",
    paste(quoted_coordinates, collapse = ", "),
    " FROM (", source_sql, ")"
  )

  codon_lookup <- ig_shm_codon_lookup(genetic_code)
  calculated <- dplyr::tbl(
    con,
    dbplyr::sql(ig_shm_per_row_sql(
      sql_input, codon_lookup, con
    ))
  ) |>
    duckplyr::as_duckdb_tibble() |>
    dplyr::compute(
      name = basename(tempfile(pattern = "immunarch_ig_shm_")),
      prudence = "stingy"
    )
  if (chain_col != "imd_chain_id") {
    calculated <- calculated |>
      dplyr::rename(!!chain_col := imd_chain_id)
  }
  annotations <- source_annotations |>
    dplyr::left_join(calculated, by = chain_col)

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


ig_shm_codon_lookup <- function(genetic_code) {
  codon_pairs <- expand.grid(
    germ_codon = names(genetic_code),
    observed_codon = names(genetic_code),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  counts <- Map(function(germ_codon, observed_codon) {
    germ_bases <- strsplit(germ_codon, "", fixed = TRUE)[[1L]]
    observed_bases <- strsplit(observed_codon, "", fixed = TRUE)[[1L]]
    mutant_codons <- vapply(seq_len(3L), function(position) {
      mutant <- germ_bases
      mutant[[position]] <- observed_bases[[position]]
      paste0(mutant, collapse = "")
    }, character(1))
    changed <- germ_bases != observed_bases
    mutant_amino_acids <- unname(genetic_code[mutant_codons])
    classifiable <- changed & mutant_amino_acids != "*"
    c(
      r_count = sum(
        classifiable & mutant_amino_acids != genetic_code[[germ_codon]]
      ),
      s_count = sum(
        classifiable & mutant_amino_acids == genetic_code[[germ_codon]]
      )
    )
  }, codon_pairs$germ_codon, codon_pairs$observed_codon)
  counts <- do.call(rbind, counts)
  codon_pairs$r_count <- as.integer(counts[, "r_count"])
  codon_pairs$s_count <- as.integer(counts[, "s_count"])
  codon_pairs
}


ig_shm_per_row_sql <- function(source, codon_lookup, con) {
  regions <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
  region_values <- paste0(
    "struct_pack(region := '", regions,
    "', region_start := cast(", regions,
    "_start AS BIGINT), region_end := cast(", regions, "_end AS BIGINT))",
    collapse = ", "
  )

  direct_codons <- paste0("WITH region_bounds AS (
    SELECT *, list_transform([", region_values, "], r -> struct_pack(
      region := r.region,
      region_start := r.region_start,
      query_start := CASE
        WHEN r.region_start IS NULL OR r.region_end IS NULL THEN NULL
        ELSE greatest(r.region_start, cast(v_sequence_start AS BIGINT))
      END,
      query_end := CASE
        WHEN r.region_start IS NULL OR r.region_end IS NULL THEN NULL
        ELSE least(
          r.region_end,
          cast(v_sequence_start AS BIGINT) + length(observed_alignment)
            - cast(v_alignment_start AS BIGINT)
        )
      END
    )) AS region_values
    FROM input
    WHERE observed_alignment IS NOT NULL AND germline_alignment IS NOT NULL
      AND v_sequence_start IS NOT NULL AND v_alignment_start IS NOT NULL
      AND NOT regexp_matches(observed_alignment, '[-.]')
  ), codon_bounds AS (
    SELECT * EXCLUDE (region_values), list_transform(region_values, r -> struct_pack(
      region := r.region,
      region_start := r.region_start,
      query_start := r.query_start,
      query_end := r.query_end,
      first_codon := CASE
        WHEN r.query_start IS NULL OR r.query_start > r.query_end THEN 0
        ELSE floor((r.query_start - r.region_start) / 3)::BIGINT
      END,
      codon_end := CASE
        WHEN r.query_start IS NULL OR r.query_start > r.query_end THEN 0
        ELSE floor((r.query_end - r.region_start) / 3)::BIGINT + 1
      END
    )) AS region_values
    FROM region_bounds
  ), nested_codons AS (
    SELECT imd_chain_id, list_transform(
      region_values,
      r -> list_transform(range(r.first_codon, r.codon_end), codon_index ->
        struct_pack(
          region := r.region,
          observed_codon := substr(
            observed_alignment,
            cast(v_alignment_start AS BIGINT)
              + greatest(r.query_start, r.region_start + 3 * codon_index)
              - cast(v_sequence_start AS BIGINT),
            least(r.query_end, r.region_start + 3 * codon_index + 2)
              - greatest(r.query_start, r.region_start + 3 * codon_index) + 1
          ),
          germ_codon := substr(
            germline_alignment,
            cast(v_alignment_start AS BIGINT)
              + greatest(r.query_start, r.region_start + 3 * codon_index)
              - cast(v_sequence_start AS BIGINT),
            least(r.query_end, r.region_start + 3 * codon_index + 2)
              - greatest(r.query_start, r.region_start + 3 * codon_index) + 1
          ),
          gap_free_codon := TRUE
        )
      )
    ) AS region_codons
    FROM codon_bounds
  )
  SELECT imd_chain_id, flatten(region_codons) AS codon_values
  FROM nested_codons")

  list_codons <- paste0("WITH position_lists AS (
    SELECT *, list_transform(
      list_filter(
        list_transform(
          range(cast(v_alignment_start AS BIGINT), length(observed_alignment) + 1),
          pos -> struct_pack(
            pos := pos,
            obs_base := substr(observed_alignment, pos, 1),
            germ_base := substr(germline_alignment, pos, 1)
          )
        ),
        x -> x.obs_base NOT IN ('-', '.')
      ),
      (x, index) -> struct_pack(
        pos := x.pos,
        obs_base := x.obs_base,
        germ_base := x.germ_base,
        query_pos := cast(v_sequence_start AS BIGINT) + index - 1
      )
    ) AS positions
    FROM input
    WHERE observed_alignment IS NOT NULL AND germline_alignment IS NOT NULL
      AND (v_sequence_start IS NULL OR v_alignment_start IS NULL
        OR regexp_matches(observed_alignment, '[-.]'))
  ), region_lists AS (
    SELECT imd_chain_id, list_transform([", region_values, "], r -> struct_pack(
      region := r.region,
      region_start := r.region_start,
      positions := list_filter(
        positions,
        p -> p.query_pos BETWEEN r.region_start AND r.region_end
      )
    )) AS region_values
    FROM position_lists
  ), codon_lists AS (
    SELECT imd_chain_id, list_transform(region_values, r -> struct_pack(
      region := r.region,
      region_start := r.region_start,
      positions := r.positions,
      codon_indices := list_sort(list_distinct(list_transform(
        r.positions,
        p -> floor((p.query_pos - r.region_start) / 3)
      )))
    )) AS region_values
    FROM region_lists
  ), nested_codons AS (
    SELECT imd_chain_id, list_transform(region_values, r -> list_transform(
      r.codon_indices,
      codon_index -> struct_pack(
        region := r.region,
        observed_codon := list_reduce(list_transform(list_filter(
          r.positions,
          p -> floor((p.query_pos - r.region_start) / 3) = codon_index
        ), p -> p.obs_base),
          (left_value, right_value) -> left_value || right_value, ''),
        germ_codon := list_reduce(list_transform(list_filter(
          r.positions,
          p -> floor((p.query_pos - r.region_start) / 3) = codon_index
        ), p -> p.germ_base),
          (left_value, right_value) -> left_value || right_value, ''),
        gap_free_codon :=
          list_max(list_transform(list_filter(
            r.positions,
            p -> floor((p.query_pos - r.region_start) / 3) = codon_index
          ), p -> p.pos))
            - list_min(list_transform(list_filter(
              r.positions,
              p -> floor((p.query_pos - r.region_start) / 3) = codon_index
            ), p -> p.pos)) + 1
            = list_count(list_filter(
              r.positions,
              p -> floor((p.query_pos - r.region_start) / 3) = codon_index
            ))
          AND list_count(list_filter(
            r.positions,
            p -> floor((p.query_pos - r.region_start) / 3) = codon_index
              AND p.germ_base IN ('-', '.')
          )) = 0
      )
    )) AS region_codons
    FROM codon_lists
  )
  SELECT imd_chain_id, flatten(region_codons) AS codon_values
  FROM nested_codons")

  aggregate_fields <- function(label, selected_regions = NULL) {
    filter <- if (is.null(selected_regions)) {
      ""
    } else {
      paste0(
        " FILTER (WHERE region IN (",
        paste0("'", selected_regions, "'", collapse = ", "), "))"
      )
    }
    available <- paste0("count(denominator)", filter)
    denominator <- paste0("sum(denominator)", filter)
    r_count <- paste0("sum(r_count)", filter)
    s_count <- paste0("sum(s_count)", filter)
    paste0(
      "CASE WHEN ", available, " = 0 THEN NULL::BIGINT ELSE ", r_count,
      "::BIGINT END AS mu_count", label, "_r, ",
      "CASE WHEN ", available, " = 0 THEN NULL::BIGINT ELSE ", s_count,
      "::BIGINT END AS mu_count", label, "_s, ",
      "CASE WHEN coalesce(", denominator, ", 0) = 0 THEN NULL::DOUBLE ELSE ",
      r_count, "::DOUBLE / ", denominator, " END AS mu_freq", label, "_r, ",
      "CASE WHEN coalesce(", denominator, ", 0) = 0 THEN NULL::DOUBLE ELSE ",
      s_count, "::DOUBLE / ", denominator, " END AS mu_freq", label, "_s"
    )
  }
  regional_fields <- unlist(lapply(regions, function(region) {
    c(
      paste0(
        "max(r_count) FILTER (WHERE region = '", region,
        "') AS mu_count_", region, "_r"
      ),
      paste0(
        "max(s_count) FILTER (WHERE region = '", region,
        "') AS mu_count_", region, "_s"
      ),
      paste0(
        "max(r_freq) FILTER (WHERE region = '", region,
        "') AS mu_freq_", region, "_r"
      ),
      paste0(
        "max(s_freq) FILTER (WHERE region = '", region,
        "') AS mu_freq_", region, "_s"
      )
    )
  }), use.names = FALSE)
  result_fields <- paste(c(
    "CASE WHEN count(denominator) = 0 THEN NULL::BIGINT
      ELSE sum(r_count + s_count)::BIGINT END AS mu_count",
    "CASE WHEN coalesce(sum(denominator), 0) = 0 THEN NULL::DOUBLE
      ELSE sum(r_count + s_count)::DOUBLE / sum(denominator) END AS mu_freq",
    aggregate_fields(""),
    aggregate_fields("_cdr", paste0("cdr", 1:3)),
    aggregate_fields("_fwr", paste0("fwr", 1:4)),
    regional_fields
  ), collapse = ", ")
  lookup_values <- paste0(
    "(", DBI::dbQuoteString(con, codon_lookup$germ_codon), ", ",
    DBI::dbQuoteString(con, codon_lookup$observed_codon), ", ",
    codon_lookup$r_count, ", ", codon_lookup$s_count, ")",
    collapse = ", "
  )

  paste0("WITH input AS MATERIALIZED (", source, "),
    codon_lookup(germ_codon, observed_codon, r_count, s_count) AS (
      VALUES ", lookup_values, "
    ),
    codon_input AS (
      SELECT * FROM (", direct_codons, ")
      UNION ALL SELECT * FROM (", list_codons, ")
    ), codon_rows AS (
      SELECT c.imd_chain_id, codon
      FROM codon_input c
      CROSS JOIN UNNEST(c.codon_values) u(codon)
    ), classified_codons AS (
      SELECT c.imd_chain_id, c.codon.region,
        list_count(list_filter(list_zip(
          string_split(c.codon.observed_codon, ''),
          string_split(c.codon.germ_codon, '')
        ), bases ->
          bases[1] IN ('A', 'C', 'G', 'T')
            AND bases[2] IN ('A', 'C', 'G', 'T')
        )) AS denominator,
        CASE WHEN c.codon.gap_free_codon
          THEN coalesce(l.r_count, 0) ELSE 0 END AS r_count,
        CASE WHEN c.codon.gap_free_codon
          THEN coalesce(l.s_count, 0) ELSE 0 END AS s_count
      FROM codon_rows c
      LEFT JOIN codon_lookup l
        ON c.codon.germ_codon = l.germ_codon
        AND c.codon.observed_codon = l.observed_codon
      WHERE c.codon.region IS NOT NULL
    ), calculated_region_stats AS (
      SELECT imd_chain_id, region,
        sum(denominator)::BIGINT AS denominator,
        sum(r_count)::BIGINT AS r_count,
        sum(s_count)::BIGINT AS s_count
      FROM classified_codons
      GROUP BY imd_chain_id, region
    ), regions AS (
      SELECT i.imd_chain_id,
        i.observed_alignment IS NOT NULL AND i.germline_alignment IS NOT NULL
          AS alignments_available,
        r.region, r.region_start, r.region_end
      FROM input i,
      LATERAL (VALUES
        ('fwr1', cast(i.fwr1_start AS BIGINT), cast(i.fwr1_end AS BIGINT)),
        ('cdr1', cast(i.cdr1_start AS BIGINT), cast(i.cdr1_end AS BIGINT)),
        ('fwr2', cast(i.fwr2_start AS BIGINT), cast(i.fwr2_end AS BIGINT)),
        ('cdr2', cast(i.cdr2_start AS BIGINT), cast(i.cdr2_end AS BIGINT)),
        ('fwr3', cast(i.fwr3_start AS BIGINT), cast(i.fwr3_end AS BIGINT)),
        ('cdr3', cast(i.cdr3_start AS BIGINT), cast(i.cdr3_end AS BIGINT)),
        ('fwr4', cast(i.fwr4_start AS BIGINT), cast(i.fwr4_end AS BIGINT))
      ) r(region, region_start, region_end)
    ), region_stats AS (
      SELECT r.imd_chain_id, r.region,
        CASE WHEN NOT r.alignments_available OR r.region_start IS NULL
            OR r.region_end IS NULL THEN NULL
          ELSE coalesce(c.denominator, 0)::BIGINT END AS denominator,
        CASE WHEN NOT r.alignments_available OR r.region_start IS NULL
            OR r.region_end IS NULL THEN NULL
          ELSE coalesce(c.r_count, 0)::BIGINT END AS r_count,
        CASE WHEN NOT r.alignments_available OR r.region_start IS NULL
            OR r.region_end IS NULL THEN NULL
          ELSE coalesce(c.s_count, 0)::BIGINT END AS s_count
      FROM regions r
      LEFT JOIN calculated_region_stats c
        ON r.imd_chain_id = c.imd_chain_id AND r.region = c.region
    ), frequencies AS (
      SELECT *,
        CASE WHEN denominator = 0 THEN NULL::DOUBLE
          ELSE r_count::DOUBLE / denominator END AS r_freq,
        CASE WHEN denominator = 0 THEN NULL::DOUBLE
          ELSE s_count::DOUBLE / denominator END AS s_freq
      FROM region_stats
    )
    SELECT imd_chain_id, ", result_fields, "
    FROM frequencies GROUP BY imd_chain_id")
}


#' @export
annotate_ig_shm <- register_method_annotate(
  core = annotate_ig_shm_impl,
  family = "annotate_ig",
  name = "shm",
  required_cols = c("sequence_col", "germline_col"),
  need_repertoires = FALSE
)
