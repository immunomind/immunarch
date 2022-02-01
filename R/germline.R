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
#' @importFrom stringr str_sub
#' @importFrom purrr modify
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
#' Data with added columns V.first.allele and J.first.allele (with first alleles of the genes);
#' V.sequence and J.sequence (with V/J germline sequences)
#' and Germline.Sequence (with combined germline sequence)
#'
#' @examples
#'
#' data(immdata)
#' repGermline(immdata$data)
#' @export repGermline
repGermline <- function(.data, species = "HomoSapiens") {
  if (inherits(.data, "list")) {
    .validate_repertoires_data(.data)
    .data %>%
      seq_along() %>%
      lapply(function(i) {
        .data[[i]] %>%
          as_tibble() %>%
          germline_single_df(species, sample_name = names(.data)[i]) %>%
          return()
      }) %>%
      return()
  } else {
    .data %>%
      as_tibble() %>%
      germline_single_df(species) %>%
      return()
  }
}

germline_single_df <- function(data, species, sample_name = NA) {
  data["V.first.allele"] <- data %>%
    select(V.name) %>%
    apply(MARGIN = 1, FUN = take_first_allele)
  data["J.first.allele"] <- data %>%
    select(J.name) %>%
    apply(MARGIN = 1, FUN = take_first_allele)

  data %<>% merge_reference_sequences("IGHV", species, sample_name)
  data %<>% merge_reference_sequences("IGHJ", species, sample_name)

  data %>%
    mutate(Germline.Sequence = purrr::map2(
      V.sequence, J.sequence, generate_germline_sequence
    )) %>%
    return()
}

take_first_allele <- function(string) {
  unlist(strsplit(string, ","))[1]
}

generate_germline_sequence <- function(v_seq, j_seq) {
  if (is.na(v_seq) || is.na(j_seq)) {
    return(NA)
  } else {
    return(paste(v_seq, j_seq, sep = "..."))
  }
}

merge_reference_sequences <- function(data, chain, species, sample_name) {
  data(genesegments)
  reference_df <- GENE_SEGMENTS %>% filter(species == species)
  reference_df <- reference_df[c("sequence", "allele_id")]
  chain_letter <- stringr::str_sub(chain, -1)
  chain_seq_colname <- paste0(chain_letter, ".sequence")
  chain_allele_colname <- paste0(chain_letter, ".first.allele")
  colnames(reference_df) <- c(chain_seq_colname, chain_allele_colname)

  alleles_in_data <- unique(data[[chain_allele_colname]])
  alleles_in_ref <- unique(reference_df[[chain_allele_colname]])
  missing_alleles <- alleles_in_data[!(alleles_in_data %in% alleles_in_ref)]
  if (length(missing_alleles) > 0) {
    warning(
      "Alleles ", paste(missing_alleles, collapse = ", "),
      if (is.na(sample_name)) "" else paste0(" from sample ", sample_name),
      " not found in the reference and will be dropped!\n",
      "Probably, species argument is wrong (current value: ", species,
      ") or the data contains non-BCR genes."
    )
  }

  data %>%
    merge(reference_df, by = chain_allele_colname) %>%
    return()
}
