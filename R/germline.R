#' This function creates germlines for clonal lineages. B cell clonal lineage represents a set of B cells
#' that presumably have a common origin (arising from the same VDJ rearrangement event) and a common ancestor.
#' Each clonal lineage has its own germline sequence that represents the ancestral sequence
#' for each BCR in clonal lineage. In other words, germline sequence is a sequence of B-cells immediately
#' after VDJ recombination, before B-cell maturation and hypermutation process. Germline sequence is useful
#' for assessing the degree of mutation and maturity of the repertoire.
#'
#' @concept germline
#'
#' @aliases repGermline germline_single_df generate_germline_sequence merge_reference_sequences validate_genes_edges validate_chains_length
#'
#' @importFrom stringr str_sub str_length str_replace fixed
#' @importFrom purrr imap
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom tidyr drop_na
#' @importFrom dplyr rowwise

#' @description Creates germlines for clonal lineages
#'
#' @usage
#'
#' repGermline(.data, species)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}.
#'
#' @param species Species from which the data was acquired. Available options:
#' "HomoSapiens" (default), "MusMusculus", "BosTaurus", "CamelusDromedarius",
#' "CanisLupusFamiliaris", "DanioRerio", "MacacaMulatta", "MusMusculusDomesticus",
#' "MusMusculusCastaneus", "MusMusculusMolossinus", "MusMusculusMusculus", "MusSpretus",
#' "OncorhynchusMykiss", "OrnithorhynchusAnatinus", "OryctolagusCuniculus", "RattusNorvegicus",
#' "SusScrofa".
#'
#' @param min_nuc_outside_cdr3 This parameter sets how many nucleotides should have V or J chain
#' outside of CDR3 to be considered good for further alignment.
#'
#' @return
#'
#' Data with added columns V.first.allele, J.first.allele (with first alleles of V and J genes),
#' V.sequence, J.sequence (with V and J reference sequences),
#' Germline.sequence (with combined germline sequence)
#'
#' @examples
#'
#' data(bcrdata)
#'
#' bcrdata$data %>%
#'   top(2000) %>% # reduce the dataset to save time on examples
#'   repGermline()
#' @export repGermline
repGermline <- function(.data, species = "HomoSapiens", min_nuc_outside_cdr3 = 5) {
  .data %<>%
    apply_to_sample_or_list(
      germline_single_df,
      .with_names = TRUE,
      species = species,
      min_nuc_outside_cdr3 = min_nuc_outside_cdr3
    )
  return(.data)
}

germline_single_df <- function(data, species, min_nuc_outside_cdr3, sample_name = NA) {
  data %<>%
    validate_genes_edges(sample_name) %>%
    add_column_with_first_gene(
      "V.name",
      "V.first.allele",
      .with_allele = TRUE
    ) %>%
    merge_reference_sequences("V", species, sample_name) %>%
    add_column_with_first_gene(
      "J.name",
      "J.first.allele",
      .with_allele = TRUE
    ) %>%
    merge_reference_sequences("J", species, sample_name) %>%
    validate_chains_length(min_nuc_outside_cdr3, sample_name) %>%
    rowwise() %>%
    mutate(Germline.sequence = generate_germline_sequence(
      Sequence, V.sequence, J.sequence, V.end, CDR3.start, CDR3.end, J.start, sample_name
    )) %>%
    drop_na(Germline.sequence)
  return(data)
}

generate_germline_sequence <- function(seq, v_ref, j_ref, v_end, cdr3_start, cdr3_end, j_start, sample_name) {
  if (any(is.na(c(seq, v_ref, j_ref, v_end, cdr3_start, cdr3_end, j_start))) || (seq == "")) {
    warning(
      "Some of mandatory fields in a row ",
      optional_sample("from sample ", sample_name, " "),
      "contain unexpected NA or empty strings! Found values:\n",
      "Sequence = \"",
      seq,
      "\",\nV.sequence = \"",
      v_ref,
      "\",\nJ.sequence = \"",
      j_ref,
      "\",\nV.end = ",
      v_end,
      ", CDR3.start = ",
      cdr3_start,
      ", CDR3.end = ",
      cdr3_end,
      ", J.start = ",
      j_start,
      ".\nThe row will be dropped!"
    )
    return(NA)
  } else {
    cdr3_start %<>% as.numeric()
    cdr3_end %<>% as.numeric()

    if (v_end <= cdr3_start) {
      v_part <- v_ref
    } else {
      # trim intersection of V and CDR3 from reference V gene
      v_part <- stringr::str_sub(
        v_ref, 1,
        max(0, stringr::str_length(v_ref) - (v_end - cdr3_start))
      )
    }

    cdr3_part <- paste(rep("n", cdr3_end - cdr3_start), collapse = "")

    if (j_start >= cdr3_end) {
      j_part <- j_ref
    } else {
      # trim intersection of J and CDR3 from reference J gene
      j_part <- stringr::str_sub(j_ref, cdr3_end - j_start + 1)
    }

    germline <- paste0(v_part, cdr3_part, j_part) %>%
      toupper()
    return(germline)
  }
}

merge_reference_sequences <- function(data, chain_letter, species, sample_name) {
  data(genesegments)
  reference_df <- GENE_SEGMENTS %>% filter(species == species)
  reference_df <- reference_df[c("sequence", "allele_id")]
  chain_seq_colname <- paste0(chain_letter, ".sequence")
  chain_allele_colname <- paste0(chain_letter, ".first.allele")
  colnames(reference_df) <- c(chain_seq_colname, chain_allele_colname)

  # check for alleles in data that don't exist in the reference
  alleles_in_data <- unique(data[[chain_allele_colname]])
  alleles_in_ref <- unique(reference_df[[chain_allele_colname]])
  missing_alleles <- alleles_in_data[!(alleles_in_data %in% alleles_in_ref)]
  if (length(missing_alleles) > 0) {
    warning(
      "Alleles ",
      paste(missing_alleles, collapse = ", "),
      " ",
      optional_sample("from sample ", sample_name, " "),
      "not found in the reference and will be dropped!\n",
      "Probably, species argument is wrong (current value: ",
      species,
      ") or the data contains non-BCR genes."
    )
  }

  data %<>% merge(reference_df, by = chain_allele_colname)
  if (nrow(data) == 0) {
    stop(
      "After merging with reference, the data ",
      optional_sample("from sample ", sample_name, " "),
      "is empty.\n",
      "There were no valid alleles in the data!"
    )
  }
  return(data)
}

validate_genes_edges <- function(data, sample_name) {
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty!"
    )
  }
  for (column in c("V.end", "J.start")) {
    if (!(column %in% colnames(data))) {
      stop(
        "Missing mandatory ",
        column,
        " column",
        optional_sample(" in sample ", sample_name, ""),
        "!"
      )
    }
    if (all(is.na(data[, column]))) {
      stop(
        "No data in mandatory ",
        column,
        " column",
        optional_sample(" in sample ", sample_name, ""),
        "!"
      )
    }
  }
  old_length <- nrow(data)
  data %<>% drop_na(V.end, J.start)
  dropped_num <- old_length - nrow(data)
  if (dropped_num > 0) {
    warning(
      dropped_num,
      " rows from ",
      old_length,
      optional_sample(" in sample ", sample_name, ""),
      " were dropped because of missing values in mandatory columns V.end and J.start!"
    )
  }
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty after dropping missing values!"
    )
  }
  return(data)
}

validate_chains_length <- function(data, min_nuc_outside_cdr3, sample_name) {
  old_length_v <- nrow(data)
  data %<>% filter(v_len_outside_cdr3(V.end, CDR3.start) >= min_nuc_outside_cdr3)
  dropped_v <- old_length_v - nrow(data)
  old_length_j <- nrow(data)
  if (nrow(data) > 0) {
    data %<>% filter(j_len_outside_cdr3(Sequence, J.start, CDR3.end) >= min_nuc_outside_cdr3)
  }
  dropped_j <- old_length_j - nrow(data)

  warning_prefix <- paste0(
    " clonotype(s) ",
    optional_sample("in sample ", sample_name, " "),
    "were dropped because they have too short part of "
  )
  warning_suffix <- paste0(
    " chain that doesn't intersect with CDR3!\n",
    "Lengths less than ",
    min_nuc_outside_cdr3,
    " are considered too short."
  )
  if (dropped_v > 0) {
    warning(dropped_v, warning_prefix, "V", warning_suffix)
  }
  if (dropped_j > 0) {
    warning(dropped_j, warning_prefix, "J", warning_suffix)
  }
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty after dropping sequences with too short V and J chains!"
    )
  }
  return(data)
}
