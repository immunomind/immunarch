seqGermline <- function(.data) {
  data <- .data
  # add first allele of V and J genes
  data["V_allele"] <- data %>%
    select(V_CALL) %>%
    apply(MARGIN = 1, FUN = take_first_allele)
  data["J_allele"] <- data %>%
    select(J_CALL) %>%
    apply(MARGIN = 1, FUN = take_first_allele)

  # load V genes
  v_genes <- load_reference_sequences("IGHV")

  # add V genes
  data <- merge(v_genes, data, by = "V_allele")

  # load J genes
  j_genes <- load_reference_sequences("IGHJ")

  # add J genes
  data <- merge(j_genes, data, by = "J_allele")

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
  colnames(sequences_df) <- c(paste0(chain_letter, "_sequence"), paste0(chain_letter, "_allele"))
  return(sequences_df)
}
