#' This function creates germlines for clonal lineages. B cell clonal lineage represents a set of B cells
#' that presumably have a common origin (arising from the same VDJ rearrangement event) and a common ancestor.
#' Each clonal lineage has its own germline sequence that represents the ancestral sequence
#' for each BCR in clonal lineage. In other words, germline sequence is a sequence of B-cells immediately
#' after VDJ recombination, before B-cell maturation and hypermutation process. Germline sequence is useful
#' for assessing the degree of mutation and maturity of the repertoire.
#'
#' @concept germline
#'
#' @aliases repGermline germline_single_df take_first_allele extract_gene_name generate_germline_sequence load_reference_sequences
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
#' @param species Species from which the data was acquired. Available options:
#' "HomoSapiens" (default), "MusMusculus", "BosTaurus", "CamelusDromedarius",
#' "CanisLupusFamiliaris", "DanioRerio", "MacacaMulatta", "MusMusculusDomesticus",
#' "MusMusculusCastaneus", "MusMusculusMolossinus", "MusMusculusMusculus", "MusSpretus",
#' "OncorhynchusMykiss", "OrnithorhynchusAnatinus", "OryctolagusCuniculus", "RattusNorvegicus",
#' "SusScrofa".
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}.
#'
#' @return
#'
#' Data with added columns V.first.allele (with first allele of V gene),
#' V.sequence (with V reference sequence) and Germline.sequence (with combined germline sequence)
#'
#' @examples
#'
#' data(immdata)
#'
#' immdata$data %>%
#'   top(2000) %>% # reduce the dataset to save time on examples
#'   repGermline()
#' @export repGermline
repGermline <- function(.data, species = "HomoSapiens") {
  if (inherits(.data, "list")) {
    .validate_repertoires_data(.data)
    .data %<>%
      purrr::imap(function(sample_data, sample_name) {
        sample_data %>%
          as_tibble() %>%
          germline_single_df(species, sample_name)
      })
    return(.data)
  } else {
    .data %<>%
      as_tibble() %>%
      germline_single_df(species)
    return(.data)
  }
}

germline_single_df <- function(data, species, sample_name = NA) {
  data %>%
    rowwise() %>%
    mutate(V.first.allele = take_first_allele(V.name)) %>%
    merge_reference_sequences("V", species, sample_name) %>%
    rowwise() %>%
    mutate(Germline.sequence = generate_germline_sequence(
      Sequence, V.sequence, V.end, CDR3.start, CDR3.end, J.start, sample_name
    )) %>%
    drop_na(Germline.sequence)
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

generate_germline_sequence <- function(seq, v_ref, v_end, cdr3_start, cdr3_end, j_start, sample_name) {
  if (any(is.na(c(seq, v_ref, v_end, cdr3_start, cdr3_end, j_start))) || (seq == "")) {
    warning(
      "Some of mandatory fields in a row ",
      optional_from_sample(sample_name),
      "contain unexpected NA or empty strings! Found values:\n",
      "Sequence = \"",
      seq,
      "\",\nV.sequence = \"",
      v_ref,
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

    j_part <- stringr::str_sub(seq, max(cdr3_end, j_start))

    paste0(v_part, cdr3_part, j_part) %>%
      toupper()
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
      optional_from_sample(sample_name),
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
      optional_from_sample(sample_name),
      "is empty.\n",
      "There were no valid alleles, or the data was initially empty."
    )
  }
  return(data)
}
