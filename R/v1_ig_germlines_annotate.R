#' @title Reconstruct clone-level immunoglobulin germlines
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' Select one representative sequence per clone and locus, then reconstruct its
#' germline from V and J references with the intervening sequence masked by `N`.
#'
#' @section Reference compatibility:
#'
#' **Warning:** `references` must contain the same IMGT-gapped reference
#' sequences, with the same numbering scheme, that were used to produce the
#' input AIRR alignments and germline coordinates. Mixing reference sets or
#' using ungapped references can silently shift germline bases and junction
#' boundaries.
#'
#' @param idata [immundata::ImmunData]. AIRR data to annotate.
#' @param references A table returned by [read_germline_references()].
#' @param species Character scalar, or `NULL` to use the sole reference species.
#' @param by Character vector or `"<auto>"`. Clone grouping columns. Automatic
#'   grouping uses `family_id`, otherwise `subject_id` and `clone_id`. Locus is
#'   always added.
#' @param sequence_id_col,sequence_col,v_col,j_col Character scalars selecting
#'   the sequence ID, aligned sequence, V-call, and J-call columns.
#' @param v_start_col,v_end_col,j_start_col,j_end_col Character scalars selecting
#'   the closed germline-coordinate columns.
#' @param germline_col Character scalar naming the reconstructed output column.
#'
#' @return [immundata::ImmunData] with one clone/locus germline assigned to each
#'   input row. Groups which cannot be reconstructed receive `NA`.
#' @name annotate_ig_germlines
NULL


#' @keywords internal
annotate_ig_germlines_impl <- function(
  idata,
  references,
  species = NULL,
  by = "<auto>",
  sequence_id_col = "sequence_id",
  sequence_col = "sequence_alignment",
  v_col = "v_call",
  j_col = "j_call",
  v_start_col = "v_germline_start",
  v_end_col = "v_germline_end",
  j_start_col = "j_germline_start",
  j_end_col = "j_germline_end",
  germline_col = "germline_alignment"
) {
  reference_cols <- c("species", "locus", "segment", "allele_id", "sequence")
  if (!is.data.frame(references) || !all(reference_cols %in% names(references))) {
    cli::cli_abort(
      "{.arg references} must contain {.field {reference_cols}} columns."
    )
  }

  reference_species <- unique(references$species)
  reference_species <- reference_species[!is.na(reference_species)]
  if (is.null(species)) {
    if (length(reference_species) != 1L) {
      cli::cli_abort(
        "{.arg species} must be supplied when references contain multiple species."
      )
    }
    species <- reference_species[[1L]]
  } else {
    checkmate::assert_string(species)
    if (!species %in% reference_species) {
      cli::cli_abort("Unknown {.arg species}: {.val {species}}.")
    }
  }

  selected_references <- references |>
    dplyr::filter(species == !!species, segment %in% c("V", "J")) |>
    dplyr::select(dplyr::all_of(reference_cols))

  annotation_cols <- names(idata$annotations)
  locus_col <- immundata::imd_schema("locus")
  if (identical(by, "<auto>")) {
    if ("family_id" %in% annotation_cols) {
      by <- "family_id"
    } else if (all(c("subject_id", "clone_id") %in% annotation_cols)) {
      by <- c("subject_id", "clone_id")
    } else {
      cli::cli_abort(
        "{.arg by = \"<auto>\"} requires {.field family_id} or both {.field subject_id} and {.field clone_id}."
      )
    }
  } else {
    checkmate::assert_character(
      by, min.len = 1L, any.missing = FALSE, unique = TRUE
    )
    if ("<auto>" %in% by) {
      cli::cli_abort("{.code \"<auto>\"} must be supplied alone through {.arg by}.")
    }
  }
  by <- unique(c(by, locus_col))

  input_cols <- c(
    by, sequence_id_col, sequence_col, v_col, j_col,
    v_start_col, v_end_col, j_start_col, j_end_col
  )
  missing_cols <- setdiff(input_cols, annotation_cols)
  if (length(missing_cols)) {
    cli::cli_abort(
      "Passed column name(s) [{.code {missing_cols}}] is not in the input ImmunData."
    )
  }

  first_allele <- function(column) {
    rlang::expr(dd$trim(dd$regexp_extract(
      !!rlang::sym(column), "[^,|[:space:]][^,|]*", 0L
    )))
  }
  annotations <- idata$annotations |>
    dplyr::select(dplyr::all_of(unique(input_cols))) |>
    dplyr::mutate(
      .ig_v = !!first_allele(v_col),
      .ig_j = !!first_allele(j_col),
      .ig_id = dd$concat("", !!rlang::sym(sequence_id_col)),
      .ig_length = dd$length(!!rlang::sym(sequence_col)),
      .ig_max_length = max(.ig_length, na.rm = TRUE),
      .by = dplyr::all_of(by)
    )

  modal_calls <- function(segment_col, output_col) {
    annotations |>
      dplyr::filter(
        !is.na(!!rlang::sym(segment_col)),
        !!rlang::sym(segment_col) != ""
      ) |>
      dplyr::summarise(
        .ig_n = dplyr::n(),
        .by = dplyr::all_of(c(by, segment_col))
      ) |>
      dplyr::mutate(
        .ig_max_n = max(.ig_n, na.rm = TRUE),
        .by = dplyr::all_of(by)
      ) |>
      dplyr::filter(.ig_n == .ig_max_n) |>
      dplyr::select(
        dplyr::all_of(by),
        !!output_col := !!rlang::sym(segment_col)
      )
  }
  modal_v <- modal_calls(".ig_v", ".ig_modal_v")
  modal_j <- modal_calls(".ig_j", ".ig_modal_j")

  representatives <- annotations |>
    dplyr::inner_join(
      modal_v,
      by = c(stats::setNames(by, by), ".ig_v" = ".ig_modal_v")
    ) |>
    dplyr::inner_join(
      modal_j,
      by = c(stats::setNames(by, by), ".ig_j" = ".ig_modal_j")
    ) |>
    dplyr::mutate(
      .ig_at_max = .ig_length == .ig_max_length,
      .ig_has_max = max(.ig_at_max, na.rm = TRUE),
      .by = dplyr::all_of(by)
    ) |>
    dplyr::filter(!.ig_has_max | .ig_at_max) |>
    dplyr::arrange(!!!rlang::syms(by), .ig_id) |>
    dplyr::mutate(
      .ig_rank = dplyr::row_number(),
      .by = dplyr::all_of(by)
    ) |>
    dplyr::filter(.ig_rank == 1L)

  groups <- annotations |>
    dplyr::summarise(
      .ig_any_v = any(!is.na(.ig_v) & .ig_v != ""),
      .ig_any_j = any(!is.na(.ig_j) & .ig_j != ""),
      .ig_max_length = max(.ig_length, na.rm = TRUE),
      .by = dplyr::all_of(by)
    )

  reference_table <- duckplyr::as_duckdb_tibble(
    selected_references,
    prudence = "stingy"
  )
  v_references <- reference_table |>
    dplyr::filter(segment == "V") |>
    dplyr::select(
      !!locus_col := locus,
      .ig_v = allele_id,
      .ig_v_reference = sequence
    )
  j_references <- reference_table |>
    dplyr::filter(segment == "J") |>
    dplyr::select(
      !!locus_col := locus,
      .ig_j = allele_id,
      .ig_j_reference = sequence
    )

  reconstructed <- groups |>
    dplyr::left_join(representatives, by = by, suffix = c("", ".rep")) |>
    dplyr::left_join(
      v_references,
      by = c(locus_col, ".ig_v"),
      copy = TRUE
    ) |>
    dplyr::left_join(
      j_references,
      by = c(locus_col, ".ig_j"),
      copy = TRUE
    ) |>
    dplyr::mutate(
      .ig_v_length = !!rlang::sym(v_end_col) -
        !!rlang::sym(v_start_col) + 1L,
      .ig_j_length = !!rlang::sym(j_end_col) -
        !!rlang::sym(j_start_col) + 1L,
      .ig_extension = .ig_max_length - .ig_length,
      .ig_middle_length = .ig_length - .ig_v_length - .ig_j_length,
      .ig_status = dplyr::if_else(
        !.ig_any_v | !.ig_any_j,
        "missing calls",
        dplyr::if_else(
          is.na(.ig_v) | is.na(.ig_j),
          "no modal V/J co-occurrence",
          dplyr::if_else(
            is.na(.ig_v_reference) | is.na(.ig_j_reference),
            "absent references",
            dplyr::if_else(
              is.na(!!rlang::sym(v_start_col)) |
                is.na(!!rlang::sym(v_end_col)) |
                is.na(!!rlang::sym(j_start_col)) |
                is.na(!!rlang::sym(j_end_col)) |
                !!rlang::sym(v_start_col) < 1L |
                !!rlang::sym(j_start_col) < 1L |
                .ig_v_length < 1L | .ig_j_length < 1L,
              "invalid coordinates",
              dplyr::if_else(
                is.na(.ig_middle_length) | .ig_middle_length < 0L |
                  is.na(.ig_extension) | .ig_extension < 0L,
                "impossible lengths",
                "ok"
              )
            )
          )
        )
      ),
      !!germline_col := dplyr::if_else(
        .ig_status == "ok",
        dd$concat(
          dd$rpad(
            dd$substr(
              .ig_v_reference,
              as.integer(!!rlang::sym(v_start_col)),
              as.integer(.ig_v_length)
            ),
            as.integer(.ig_v_length),
            "N"
          ),
          dd$rpad("", as.integer(.ig_middle_length), "N"),
          dd$rpad(
            dd$substr(
              .ig_j_reference,
              as.integer(!!rlang::sym(j_start_col)),
              as.integer(.ig_j_length + .ig_extension)
            ),
            as.integer(.ig_j_length + .ig_extension),
            "N"
          )
        ),
        NA_character_
      )
    ) |>
    dplyr::select(dplyr::all_of(c(by, germline_col)))

  join_by <- stats::setNames(by, by)
  immundata::annotate(
    idata,
    annotations = reconstructed,
    by = join_by,
    conflicts = "replace"
  )
}


#' @rdname annotate_ig_germlines
#' @export
annotate_ig_germlines <- register_method_annotate(
  core = annotate_ig_germlines_impl,
  family = "annotate_ig",
  name = "germlines",
  need_repertoires = FALSE
)
