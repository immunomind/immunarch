#' This function creates germlines for clonal lineages. B cell clonal lineage represents a set of B cells
#' that presumably have a common origin (arising from the same VDJ rearrangement event) and a common ancestor.
#' Each clonal lineage has its own germline sequence that represents the ancestral sequence
#' for each BCR in clonal lineage. In other words, germline sequence is a sequence of B-cells immediately
#' after VDJ recombination, before B-cell maturation and hypermutation process. Germline sequence is useful
#' for assessing the degree of mutation and maturity of the repertoire.
#'
#' @concept germline
#'
#' @aliases repGermline take_first_allele extract_gene_name load_reference_sequences
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
#' @param .data The data to be processed. Can be \link{data.frame} or \link{data.table}.
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}
#'
#' @return
#'
#' Data with added columns V.allele and J.allele (with first alleles of the genes),
#' and V.sequence and J.sequence (with germline sequences).
#'
#' @examples
#'
#' data(immdata)
#' repGermline(immdata$data)
#' @export repGermline

repGermline <- function(.data) {
  data <- as_tibble(.data)
  # add first allele of V and J genes
  data["V.allele"] <- data %>%
    select(V.name) %>%
    apply(MARGIN = 1, FUN = take_first_allele)
  data["J.allele"] <- data %>%
    select(J.name) %>%
    apply(MARGIN = 1, FUN = take_first_allele)

  # load V genes
  v_genes <- load_reference_sequences("IGHV")

  # add V genes
  data <- merge(v_genes, data, by = "V.allele")

  # load J genes
  j_genes <- load_reference_sequences("IGHJ")

  # add J genes
  data <- merge(j_genes, data, by = "J.allele")

  data %<>% mutate(germline_sequence = paste(V_sequence, J_sequence, sep = "..."))
  return(data)
}

take_first_allele <- function(string) {
  unlist(strsplit(string, ","))[1]
}

# example input: "J00256|IGHJ1*01|Homo"; example output: "IGHJ1*01"
extract_gene_name <- function(full_name) {
  unlist(strsplit(full_name, '\\|'))[2]
}

load_reference_sequences <- function(chain) {
  sequences_df <- VDJgermlines::extractSequencesR("human", chain, "IMGT", FALSE)
  sequences_df <- sequences_df[c("sequence", "names")]
  sequences_df[["names"]] <- purrr::modify(sequences_df[["names"]], extract_gene_name)
  chain_letter <- stringr::str_sub(chain, -1)
  colnames(sequences_df) <- c(paste0(chain_letter, ".sequence"), paste0(chain_letter, ".allele"))
  return(sequences_df)
}
