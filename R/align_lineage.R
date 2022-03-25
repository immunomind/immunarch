#' This function aligns all sequences incliding germline that belong to one clonal lineage and one cluster.
#' After clustering, building clonal lineage and germline, the next step is to analyze the degree of mutation
#' and maturity of each clonal lineage. This allows you to find high mature cells and cells with a large
#' number of offspring. The phylogenetic analysis will find mutations that increase the affinity of BCR.
#' Making alignment of the sequence is the first step towards sequence analysis including BCR.
#'
#' @concept align_lineage
#'
#' @aliases repAlignLineage
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom stringr str_extract_all str_sub str_length
#' @importFrom ape as.DNAbin muscle
#' @importFrom parallel mclapply detectCores

#' @description Aligns all sequences incliding germline within each clonal lineage within each cluster
#'
#' @usage
#'
#' repAlignLineage(.data, .min.lineage.sequences)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' @param .min.lineage.sequences If number of sequences in the same clonal lineage and the same
#' cluster (not including germline) is lower than this threshold, this group of sequences
#' will not be aligned and will not be used in next steps of BCR pipeline
#' (but still will be saved in the output table).
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}, and also
#' must contain 'Cluster' column, which is added by seqCluster() function, and 'Sequence.germline'
#' column, which is added by repGermline() function.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has number of rows equal to number of combinations of unique clusters
#' and unique germlines, and has these columns:
#' * Cluster: cluster name
#' * Germline: germline sequence
#' * Aligned: FALSE if this group of sequences was not aligned with lineage (.min.lineage.sequences
#'   is below the threshold); TRUE if it was aligned
#' * Alignment: DNAbin object with alignment (if Aligned=TRUE) or DNAbin object with unaligned
#'   sequences (if Aligned=FALSE)
#' * V.length: shortest length of V gene part outside of CDR3 in this group of sequences;
#'   longer V genes (including germline) are trimmed to this length before alignment
#' * J.length: shortest length of J gene part outside of CDR3 in this group of sequences;
#'   longer J genes (including germline) are trimmed to this length before alignment
#' * Sequences: nested dataframe containing all sequences for this combination
#'   of cluster and germline; it has columns
#'   Sequence, V.end, J.start, CDR3.start, CDR3.end; all values taken from the input dataframe
#'
#' @examples
#'
#' data(bcrdata)
#' bcr_data <- bcrdata$data %>% top(100) # reduce the dataset to save time on examples
#'
#' bcr_data %>%
#'   seqCluster(seqDist(bcr_data)) %>%
#'   repGermline() %>%
#'   repAlignLineage()
#' @export repAlignLineage
repAlignLineage <- function(.data, .min.lineage.sequences = 3) {
  require_system_package("muscle", error_message = paste0(
    "repAlignLineage requires MUSCLE app to be installed!\n",
    "Please download it from here: https://github.com/rcedgar/muscle/releases/latest\n",
    "or install it with your system package manager (such as apt or dnf)."
  ))
  .data %<>%
    apply_to_sample_or_list(
      align_single_df,
      .min.lineage.sequences = .min.lineage.sequences
    )
  return(.data)
}

align_single_df <- function(data, .min.lineage.sequences) {
  for (required_column in c("Cluster", "Germline.sequence")) {
    if (!(required_column %in% colnames(data))) {
      stop(
        "Found dataframe without required column ",
        required_column,
        ";\nexisting columns: ",
        str(colnames(data))
      )
    }
  }

  groups <- unique(data[, c("Cluster", "Germline.sequence")])

  # results is a dataframe, see repAlignLineage @return description
  results <- groups %>%
    apply(1, get_lineage_subset, data = data) %>%
    lapply(align_sequences,
      .min.lineage.sequences = .min.lineage.sequences
#      mc.preschedule = FALSE, mc.cores = parallel::detectCores()
    ) %>%
    map_dfr(~.)
  return(results)
}

# return data subset containing only rows with specified combination of cluster and germline
get_lineage_subset <- function(group, data) {
  cluster <- group[["Cluster"]]
  germline <- group[["Germline.sequence"]]
  data %<>% subset(Cluster == cluster & Germline.sequence == germline)
  return(data)
}

# this function returns named list containing 1 row for results dataframe
align_sequences <- function(lineage_subset, .min.lineage.sequences) {
  cluster_name <- lineage_subset[[1, "Cluster"]]
  germline_seq <- lineage_subset[[1, "Germline.sequence"]]
  aligned <- nrow(lineage_subset) >= .min.lineage.sequences
  lineage_subset[["V.lengths"]] <- v_len_outside_cdr3(
    lineage_subset[["V.end"]], lineage_subset[["CDR3.start"]]
  )
  lineage_subset[["J.lengths"]] <- j_len_outside_cdr3(
    lineage_subset[["Sequence"]], lineage_subset[["J.start"]], lineage_subset[["CDR3.end"]]
  )
  v_min_len <- min(lineage_subset[["V.lengths"]])
  j_min_len <- min(lineage_subset[["J.lengths"]])
  sequences <- lineage_subset[c("Sequence", "V.end", "J.start", "CDR3.start", "CDR3.end")]

  germline_parts <- strsplit(germline_seq, "N")[[1]]
  germline_v_len <- stringr::str_length(germline_parts[1])
  germline_j_len <- stringr::str_length(tail(germline_parts, 1))
  germline_trimmed <- trim_seq(germline_seq, germline_v_len, v_min_len, germline_j_len, j_min_len)
  clonotypes_trimmed <- trim_seq(
    lineage_subset[["Sequence"]],
    lineage_subset[["V.lengths"]],
    v_min_len,
    lineage_subset[["J.lengths"]],
    j_min_len
  )
  alignment <- convert_to_dnabin(germline_trimmed, clonotypes_trimmed)
  if (aligned) {
    alignment %<>% ape::muscle()
  }

  return(list(
    Cluster = cluster,
    Germline = germline_seq,
    Aligned = aligned,
    Alignment = alignment,
    V.length = v_min_len,
    J.length = j_min_len,
    Sequences = sequences
  ))
}

convert_to_dnabin <- function(germline_seq, clonotypes) {
  all_sequences_list <- c(list(germline = germline_seq), as.list(clonotypes))
  dnabin <- lapply(
    function(sequence) {
      sequence %>%
        stringr::str_extract_all(boundary("character")) %>%
        unlist()
    }
  ) %>%
    ape::as.DNAbin()
  return(dnabin)
}

# trim V/J tails in sequence to the specified lenghts v_min, j_min
trim_seq <- function(seq, v_len, v_min, j_len, j_min) {
  stringr::str_sub(seq, v_len - v_min + 1, -(j_len - j_min + 1))
}
