#' @title Read IMGT immunoglobulin germline references
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' V and J genes provide parts of an antibody sequence. To measure mutations
#' in a B-cell clone, we need to compare its sequences with reference versions of
#' those genes before somatic hypermutation.
#'
#' This function reads nucleotide V and J genes from an IMGT FASTA file for
#' use by [annotate_ig_germlines()]. It keeps the IMGT gaps (`.`) in the
#' reference sequences.
#'
#' @section Choosing files:
#' Use IMGT-gapped V and J FASTA files from the same IMGT release used for
#' your IgBlast run. Read each needed file and combine the
#' resulting tables. Do not use the ungapped FASTA from the searchable
#' IgBlast database: this function does not add gaps or convert positions.
#'
#' @param path Character scalar. Path to an IMGT FASTA file.
#'
#' @return A tibble with `species`, `locus`, `segment`, `allele_id`, and
#'   `sequence` columns.
#' @export
read_germline_references <- function(path) {
  checkmate::assert_string(path)
  checkmate::assert_file_exists(path, access = "r")
  if (dir.exists(path)) {
    cli::cli_abort("{.arg path} must be a single FASTA file.")
  }

  lines <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (!length(headers)) {
    cli::cli_abort("No usable nucleotide IG V/J records were found in {.arg path}.")
  }

  ends <- c(headers[-1L] - 1L, length(lines))
  records <- lapply(seq_along(headers), function(i) {
    header <- substring(lines[[headers[[i]]]], 2L)
    fields <- strsplit(header, "|", fixed = TRUE)[[1L]]
    allele <- if (length(fields) >= 2L) trimws(fields[[2L]]) else ""
    region <- if (length(fields) >= 5L) trimws(fields[[5L]]) else ""
    ig_vj <- grepl("^IG[HKL][VJ]", allele)

    if (ig_vj && length(fields) < 5L) {
      cli::cli_abort("Malformed relevant IMGT header: {.val {header}}.")
    }
    relevant <- ig_vj && region %in% c("V-REGION", "J-REGION")

    if (!relevant) {
      return(NULL)
    }
    if (length(fields) < 7L) {
      cli::cli_abort("Malformed relevant IMGT header: {.val {header}}.")
    }
    if (!grepl("(^|[[:space:]])nt$", fields[[7L]])) {
      return(NULL)
    }

    species <- trimws(fields[[3L]])
    sequence_lines <- if (headers[[i]] < ends[[i]]) {
      lines[seq.int(headers[[i]] + 1L, ends[[i]])]
    } else {
      character()
    }
    sequence <- toupper(paste0(trimws(sequence_lines), collapse = ""))
    if (!nzchar(species) || !grepl("^IG[HKL][VJ].+[*].+$", allele) ||
      !nzchar(sequence) || !grepl("^[ACGTRYSWKMBDHVN.-]+$", sequence)) {
      cli::cli_abort("Malformed relevant IMGT record: {.val {header}}.")
    }

    tibble::tibble(
      species = species,
      locus = sub("^IG([HKL]).*$", "IG\\1", allele),
      segment = substring(region, 1L, 1L),
      allele_id = allele,
      sequence = sequence
    )
  })

  references <- dplyr::bind_rows(records)
  if (!nrow(references)) {
    cli::cli_abort("No usable nucleotide IG V/J records were found in {.arg path}.")
  }

  keys <- c("species", "locus", "segment", "allele_id")
  duplicates <- references |>
    dplyr::count(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::filter(n > 1L)
  if (nrow(duplicates)) {
    cli::cli_abort(
      "Duplicate IMGT records share {nrow(duplicates)} germline reference key{?s}."
    )
  }

  references |>
    dplyr::select(dplyr::all_of(c(keys, "sequence"))) |>
    tibble::as_tibble()
}
