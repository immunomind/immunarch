#' This function creates germlines for clonal lineages. B cell clonal lineage represents a set of B cells
#' that presumably have a common origin (arising from the same VDJ rearrangement event) and a common ancestor.
#' Each clonal lineage has its own germline sequence that represents the ancestral sequence
#' for each BCR in clonal lineage. In other words, germline sequence is a sequence of B-cells immediately
#' after VDJ recombination, before B-cell maturation and hypermutation process. Germline sequence is useful
#' for assessing the degree of mutation and maturity of the repertoire.
#'
#' @concept germline
#'
#' @aliases repGermline germline_single_df take_first_allele generate_germline_sequence merge_reference_sequences auto_detect_positions auto_detect_position validate_chains_length
#'
#' @importFrom stringr str_sub str_length str_replace fixed
#' @importFrom purrr imap
#' @importFrom magrittr %>% %<>%

#' @description Creates germlines for clonal lineages
#'
#' @usage
#'
#' repGermline(.data)
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
#' @param substring_length If V.end and/or J.start values are missing in the data, this function
#' will try to auto-detect them by fuzzy searching part of reference sequence in clonotype sequence.
#' This parameter is length of reference subsequence (end part for V, start part for J) that
#' will be searched in clonotype sequence.
#'
#' @param max_mismatches If V.end and/or J.start values are missing in the data, this function
#' will try to auto-detect them by fuzzy searching part of reference sequence in clonotype sequence.
#' This parameter sets how many mismatches are allowed when performing the search.
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
repGermline <- function(.data, species = "HomoSapiens", substring_length = 8, max_mismatches = 1) {
  if (inherits(.data, "list")) {
    .validate_repertoires_data(.data)
    .data %<>%
      purrr::imap(function(sample_data, sample_name) {
        sample_data %>%
          as_tibble() %>%
          germline_single_df(species, substring_length, max_mismatches, sample_name)
      })
    return(.data)
  } else {
    .data %<>%
      as_tibble() %>%
      germline_single_df(species, substring_length, max_mismatches)
    return(.data)
  }
}

germline_single_df <- function(data, species, substring_length, max_mismatches, sample_name = NA) {
  data %<>%
    rowwise() %>%
    mutate(V.first.allele = take_first_allele(V.name)) %>%
    merge_reference_sequences("V", species, sample_name) %>%
    rowwise() %>%
    mutate(J.first.allele = take_first_allele(J.name)) %>%
    merge_reference_sequences("J", species, sample_name) %>%
    auto_detect_positions(substring_length, max_mismatches, sample_name) %>%
    validate_chains_length(sample_name) %>%
    rowwise() %>%
    mutate(Germline.sequence = generate_germline_sequence(
      Sequence, V.sequence, J.sequence, V.end, CDR3.start, CDR3.end, J.start, sample_name
    )) %>%
    drop_na(Germline.sequence)
  return(data)
}

take_first_allele <- function(string) {
  string %<>%
    # first allele is substring until first ',' or '(' in string taken from column with gene names
    strsplit(",|\\(") %>%
    unlist() %>%
    extract2(1) %>%
    # MiXCR uses *00 for unknown alleles; replace *00 to *01 to find them in reference
    stringr::str_replace(stringr::fixed("*00"), "*01")
  return(string)
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
      "There were no valid alleles, or the data was initially empty."
    )
  }
  return(data)
}

# if V.end and/or J.start are missing somewhere, try to auto-detect them
auto_detect_positions <- function(data, substring_length, max_mismatches, sample_name) {
  if (sum(is.na(data$V.end)) > 0) {
    old_size <- nrow(data)

    data %<>%
      rowwise() %>%
      mutate(V.end = auto_detect_position(
        "V", V.end, Sequence, V.sequence, substring_length, max_mismatches
      )) %>%
      drop_na(V.end)

    new_size <- nrow(data)
    size_diff <- old_size - new_size
    if (new_size == 0) {
      stop(
        "Auto-detection of V.end values failed",
        optional_sample(" in sample ", sample_name, ""),
        ":\nafter removing clonotypes with not detected V.end, data size is 0!"
      )
    } else if (size_diff > 0) {
      warning(
        "Auto-detection of V.end values failed in some clonotypes",
        optional_sample(" in sample ", sample_name, ""),
        ":\ndropped ",
        size_diff,
        " from ",
        old_size,
        " clonotypes because of not detected V.end values!"
      )
    }
  }
  if (sum(is.na(data$J.start)) > 0) {
    old_size <- nrow(data)

    data %<>%
      rowwise() %>%
      mutate(J.start = auto_detect_position(
        "J", J.start, Sequence, J.sequence, substring_length, max_mismatches
      )) %>%
      drop_na(J.start)

    new_size <- nrow(data)
    size_diff <- old_size - new_size
    if (new_size == 0) {
      stop(
        "Auto-detection of J.start values failed",
        optional_sample(" in sample ", sample_name, ""),
        "\nafter removing clonotypes with not detected J.start, data size is 0!"
      )
    } else if (size_diff > 0) {
      warning(
        "Auto-detection of J.start values failed in some clonotypes",
        optional_sample(" in sample ", sample_name, ""),
        ":\ndropped ",
        size_diff,
        " from ",
        old_size,
        " clonotypes because of not detected J.start values!"
      )
    }
  }
  return(data)
}

# if the coordinate (V.end or J.start) is missing, try to find it by fuzzy searching substring from reference
auto_detect_position <- function(chain, coordinate, seq, ref, substring_length, max_mismatches) {
  if (is.na(coordinate)) {
    # for V we search end location and for J - start location
    if (chain == "V") {
      ref_substring <- stringr::str_sub(ref, -substring_length)
    } else if (chain == "J") {
      ref_substring <- stringr::str_sub(ref, 1, substring_length)
    } else {
      stop("Wrong chain argument: ", chain)
    }

    found_coordinate <- aregexec(ref_substring, seq, max_mismatches, fixed = TRUE)[[1]][1]

    if (found_coordinate == -1) {
      return(NA)
    } else {
      if (chain == "V") {
        return(as.integer(found_coordinate + substring_length))
      } else {
        return(as.integer(found_coordinate))
      }
    }
  } else {
    return(coordinate)
  }
}

validate_chains_length <- function(data, sample_name, min_good_length = 5) {
  too_short_v_chains_num <- data %>%
    rowwise() %>%
    mutate(Too.short.V = min(V.end, as.numeric(CDR3.start)) < min_good_length) %>%
    pull(Too.short.V) %>%
    sum()
  too_short_j_chains_num <- data %>%
    rowwise() %>%
    mutate(
      Too.short.J =
        stringr::str_length(Sequence) + 1 - max(J.start, as.numeric(CDR3.end)) < min_good_length
    ) %>%
    pull(Too.short.J) %>%
    sum()
  warning_prefix <- paste0(
    too_short_v_chains_num,
    " clonotype(s) ",
    optional_sample("in sample ", sample_name, " "),
    "have too short part of "
  )
  warning_suffix <- paste0(
    " chain that doesn't intersect with CDR3!\n",
    "Lengths less than ",
    min_good_length,
    " are considered too short."
  )
  if (too_short_v_chains_num > 0) {
    warning(warning_prefix, "V", warning_suffix)
  }
  if (too_short_j_chains_num > 0) {
    warning(warning_prefix, "J", warning_suffix)
  }
  return(data)
}
