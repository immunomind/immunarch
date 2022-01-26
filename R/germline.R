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
#' @importFrom VDJgermlines extractSequencesR

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
repGermline <- function(.data) {
  if (inherits(.data, "list")) {
    .validate_repertoires_data(.data)
    .data %>%
      lapply(function(sample_data) {
        sample_data %>%
          as_tibble() %>%
          germline_single_df() %>%
          return()
      }) %>%
      return()
  } else {
    .data %>%
      as_tibble() %>%
      germline_single_df() %>%
      return()
  }
}

germline_single_df <- function(data) {
  # add first allele of V and J genes
  data["V.first.allele"] <- data %>%
    select(V.name) %>%
    apply(MARGIN = 1, FUN = take_first_allele)
  data["J.first.allele"] <- data %>%
    select(J.name) %>%
    apply(MARGIN = 1, FUN = take_first_allele)

  # load V genes
  v_genes <- load_reference_sequences("IGHV")

  # add V genes
  data <- merge(x = v_genes, y = data, by = "V.first.allele", all.y = TRUE)

  # load J genes
  j_genes <- load_reference_sequences("IGHJ")

  # add J genes
  data <- merge(x = j_genes, y = data, by = "J.first.allele", all.y = TRUE)

  data %>%
    mutate(Germline.Sequence = purrr::map2(
      V.sequence, J.sequence, generate_germline_sequence
    )) %>%
    return()
}

take_first_allele <- function(string) {
  unlist(strsplit(string, ","))[1]
}

# example input: "J00256|IGHJ1*01|Homo"; example output: "IGHJ1*01"
extract_gene_name <- function(full_name) {
  unlist(strsplit(full_name, "\\|"))[2]
}

generate_germline_sequence <- function(v_seq, j_seq) {
  if (is.na(v_seq) || is.na(j_seq)) {
    return(NA)
  } else {
    return(paste(v_seq, j_seq, sep = "..."))
  }
}

load_reference_sequences <- function(chain) {
  sequences_df <- VDJgermlines::extractSequencesR("human", chain, "IMGT", FALSE)
  sequences_df <- sequences_df[c("sequence", "names")]
  sequences_df[["names"]] <- purrr::modify(sequences_df[["names"]], extract_gene_name)
  chain_letter <- stringr::str_sub(chain, -1)
  colnames(sequences_df) <- c(paste0(chain_letter, ".sequence"), paste0(chain_letter, ".first.allele"))
  return(sequences_df)
}
