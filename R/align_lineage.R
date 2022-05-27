#' This function aligns all sequences (incliding germline) that belong to one clonal lineage and one cluster.
#' After clustering and building the clonal lineage and germline, the next step is to analyze the degree of mutation
#' and maturity of each clonal lineage. This allows for finding high mature cells and cells with a large
#' number of offspring. The phylogenetic analysis will find mutations that increase the affinity of BCR.
#' Making alignment of the sequence is the first step towards sequence analysis including BCR.
#'
#' @concept align_lineage
#'
#' @aliases repAlignLineage
#'
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom stringr str_extract_all str_sub str_length boundary
#' @importFrom plyr dlply .
#' @importFrom purrr map_dfr
#' @importFrom rlist list.remove
#' @importFrom utils str
#' @importFrom ape as.DNAbin muscle
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel mclapply

#' @description Aligns all sequences incliding germline within each clonal lineage within each cluster
#'
#' @usage
#'
#' repAlignLineage(.data,
#' .min_lineage_sequences, .prepare_threads, .align_threads, .verbose_output, .nofail)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' @param .min_lineage_sequences If number of sequences in the same clonal lineage and the same
#' cluster (not including germline) is lower than this threshold, this group of sequences
#' will not be aligned and will not be used in next steps of BCR pipeline
#' (will be saved in output table only if .verbose_output parameter is set to TRUE).
#'
#' @param .prepare_threads Number of threads to prepare results table.
#' Please note that high number can cause heavy memory usage!
#'
#' @param .align_threads Number of threads for lineage alignment.
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}, and also
#' must contain 'Cluster' column, which is added by seqCluster() function, and 'Sequence.germline'
#' column, which is added by repGermline() function.
#'
#' @param .verbose_output If TRUE, all output dataframe columns will be included (see documentation about this
#' function return), and unaligned clusters will be included in the output. Setting this to TRUE significantly
#' increases memory usage. If FALSE, only aligned clusters and columns required for repClonalFamily() calculation
#' will be included in the output.
#'
#' @param .nofail Will return NA instead of stopping if Clustal W is not installed.
#' Used to avoid raising errors in examples on computers where Clustal W is not installed.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has these columns:
#' * Cluster: cluster name
#' * Germline: germline sequence
#' * Aligned (included if .verbose_output=TRUE): FALSE if this group of sequences was not aligned with lineage
#'   (.min_lineage_sequences is below the threshold); TRUE if it was aligned
#' * Alignment: DNAbin object with alignment or DNAbin object with unaligned sequences (if Aligned=FALSE)
#' * V.length (included if .verbose_output=TRUE): shortest length of V gene part outside of CDR3 region in this
#'   group of sequences; longer V genes (including germline) are trimmed to this length before alignment
#' * J.length (included if .verbose_output=TRUE): shortest length of J gene part outside of CDR3 region in this
#'   group of sequences; longer J genes (including germline) are trimmed to this length before alignment
#' * Sequences (included if .verbose_output=TRUE): nested dataframe containing all sequences for this combination
#'   of cluster and germline; it has columns
#'   Sequence, V.end, J.start, CDR3.start, CDR3.end; all values taken from the input dataframe
#'
#' @examples
#'
#' data(bcrdata)
#' bcr_data <- bcrdata$data
#'
#' bcr_data %>%
#'   seqCluster(seqDist(bcr_data), .fixed_threshold = 3) %>%
#'   repGermline() %>%
#'   repAlignLineage(.min_lineage_sequences = 2, .align_threads = 2, .nofail = TRUE)
#' @export repAlignLineage
repAlignLineage <- function(.data,
                            .min_lineage_sequences = 3,
                            .prepare_threads = 2,
                            .align_threads = 4,
                            .verbose_output = FALSE,
                            .nofail = FALSE) {
  if (!require_system_package("clustalw", error_message = paste0(
    "repAlignLineage requires Clustal W app to be installed!\n",
    "Please download it from here: http://www.clustal.org/download/current/\n",
    "or install it with your system package manager (such as apt or dnf)."
  ), .nofail)) {
    return(NA)
  }

  doParallel::registerDoParallel(cores = .prepare_threads)
  .data %<>%
    apply_to_sample_or_list(
      align_single_df,
      .min_lineage_sequences = .min_lineage_sequences,
      .align_threads = .align_threads,
      .verbose_output = .verbose_output
    )
  return(.data)
}

align_single_df <- function(data, .min_lineage_sequences, .align_threads, .verbose_output) {
  for (required_column in c("Cluster", "Germline.sequence")) {
    if (!(required_column %in% colnames(data))) {
      stop(
        "Found dataframe without required column ",
        required_column,
        ";\nexisting columns: ",
        utils::str(colnames(data))
      )
    }
  }

  results <- data %>%
    plyr::dlply(
      .variables = .(get("Cluster"), get("Germline.sequence")),
      .fun = prepare_results_row,
      .min_lineage_sequences = .min_lineage_sequences,
      .verbose_output = .verbose_output,
      .parallel = TRUE
    ) %>%
    `[`(!is.na(.)) %>%
    unname()

  if (length(results) == 0) {
    stop("There are no lineages containing at least ", .min_lineage_sequences, " sequences!")
  }

  # only required columns are passed to alignment function to reduce consumed memory
  if (.verbose_output) {
    alignments <- lapply(results, "[", c("Aligned", "Alignment"))
  } else {
    alignments <- lapply(results, "[", "Alignment")
  }
  alignments %<>% parallel::mclapply(
    align_sequences,
    .verbose_output = .verbose_output,
    mc.preschedule = TRUE,
    mc.cores = .align_threads
  )

  return(convert_results_to_df(results, alignments, .verbose_output))
}

# this function accepts dataframe subset containing rows only for current lineage
# and returns named list containing 1 row for results dataframe
prepare_results_row <- function(lineage_subset, .min_lineage_sequences, .verbose_output) {
  cluster_name <- lineage_subset[[1, "Cluster"]]
  germline_seq <- lineage_subset[[1, "Germline.sequence"]]
  aligned <- nrow(lineage_subset) >= .min_lineage_sequences

  if (!aligned & !.verbose_output) {
    return(NA)
  }

  lineage_subset[["V.lengths"]] <- v_len_outside_cdr3(
    lineage_subset[["V.end"]], lineage_subset[["CDR3.start"]]
  )
  lineage_subset[["J.lengths"]] <- j_len_outside_cdr3(
    lineage_subset[["Sequence"]], lineage_subset[["J.start"]], lineage_subset[["CDR3.end"]]
  )

  if (.verbose_output) {
    sequences <- lineage_subset[c("Sequence", "V.end", "J.start", "CDR3.start", "CDR3.end")]
  }

  germline_parts <- strsplit(germline_seq, "N")[[1]]
  germline_v_len <- stringr::str_length(germline_parts[1])
  germline_j_len <- stringr::str_length(tail(germline_parts, 1))
  v_min_len <- min(lineage_subset[["V.lengths"]], germline_v_len)
  j_min_len <- min(lineage_subset[["J.lengths"]], germline_j_len)

  germline_trimmed <- trim_seq(germline_seq, germline_v_len, v_min_len, germline_j_len, j_min_len)
  clonotypes_trimmed <- trim_seq(
    lineage_subset[["Sequence"]],
    lineage_subset[["V.lengths"]],
    v_min_len,
    lineage_subset[["J.lengths"]],
    j_min_len
  )
  alignment <- convert_to_dnabin(germline_trimmed, clonotypes_trimmed)

  if (.verbose_output) {
    return(list(
      Cluster = cluster_name,
      Germline = germline_seq,
      Aligned = aligned,
      Alignment = alignment,
      V.length = v_min_len,
      J.length = j_min_len,
      Sequences = sequences
    ))
  } else {
    return(list(
      Cluster = cluster_name,
      Germline = germline_seq,
      Alignment = alignment
    ))
  }
}

convert_to_dnabin <- function(germline_seq, clonotypes) {
  all_sequences_list <- c(list(germline = germline_seq), as.list(clonotypes))
  dnabin <- all_sequences_list %>%
    lapply(
      function(sequence) {
        sequence %>%
          stringr::str_extract_all(stringr::boundary("character")) %>%
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

convert_results_to_df <- function(nested_results_list, nested_alignments_list, .verbose_output) {
  alignments <- nested_alignments_list %>%
    lapply(magrittr::extract2, "Alignment") %>%
    tibble(Alignment = .)
  df <- nested_results_list %>%
    lapply(rlist::list.remove, c("Alignment", "Sequences")) %>%
    purrr::map_dfr(~.) %>%
    cbind(alignments)
  if (.verbose_output) {
    sequences <- nested_results_list %>%
      lapply(magrittr::extract2, "Sequences") %>%
      tibble(Sequences = .)
    df %<>% cbind(sequences)
  }
  return(df)
}

align_sequences <- function(df_row, .verbose_output) {
  if (.verbose_output) {
    aligned <- df_row[["Aligned"]]
  } else {
    aligned <- TRUE
  }
  if (aligned) {
    df_row[["Alignment"]] %<>% ape::clustal()
  }
  return(df_row)
}
