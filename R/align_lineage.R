#' Aligns all sequences incliding germline within each clonal lineage within each cluster
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
#' @importFrom ape as.DNAbin clustal
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel mclapply

#' @description This function aligns all sequences (incliding germline) that belong to one clonal
#' lineage and one cluster. After clustering and building the clonal lineage and germline, the next
#' step is to analyze the degree of mutation and maturity of each clonal lineage. This allows for
#' finding high mature cells and cells with a large number of offspring. The phylogenetic analysis
#' will find mutations that increase the affinity of BCR. Making alignment of the sequence
#' is the first step towards sequence analysis including BCR.
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
#' must contain 'Cluster' column, which is added by seqCluster() function, and 'Germline.sequence'
#' column, which is added by repGermline() function.
#'
#' @param .verbose_output If TRUE, all output dataframe columns will be included (see documentation about this
#' function return), and unaligned clusters will be included in the output. Setting this to TRUE significantly
#' increases memory usage. If FALSE, only aligned clusters and columns required for repClonalFamily() and
#' repSomaticHypermutation() calculation will be included in the output.
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
#' * V.germline.nt: germline V gene sequence
#' * J.germline.nt: germline J gene sequence
#' * CDR3.germline.length: length of CDR3 in germline
#' * Aligned (included if .verbose_output=TRUE): FALSE if this group of sequences was not aligned with lineage
#'   (.min_lineage_sequences is below the threshold); TRUE if it was aligned
#' * Alignment: DNAbin object with alignment or DNAbin object with unaligned sequences (if Aligned=FALSE)
#' * V.length: shortest length of V gene part outside of CDR3 region in this
#'   group of sequences; longer V genes (including germline) are trimmed to this length before alignment
#' * J.length: shortest length of J gene part outside of CDR3 region in this
#'   group of sequences; longer J genes (including germline) are trimmed to this length before alignment
#' * Sequences: nested dataframe containing all sequences for this combination
#'   of cluster and germline; it has columns
#'   Sequence, Clone.ID, Clones, CDR1.nt, CDR2.nt, CDR3.nt, FR1.nt, FR2.nt, FR3.nt, FR4.nt
#'   and, if .verbose_output=TRUE, also V.end, J.start, CDR3.start, CDR3.end;
#'   all values taken from the input dataframe
#' * AA.frame.starts: start positions for amino acid translation for germline and all sequences
#'   after trimming (possible values: 1, 2 and 3)
#'
#' @examples
#'
#' data(bcrdata)
#' bcr_data <- bcrdata$data
#'
#' bcr_data %>%
#'   seqCluster(seqDist(bcr_data), .fixed_threshold = 3) %>%
#'   repGermline(.threads = 1) %>%
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
    return(get_empty_object_with_class("step_failure_ignored"))
  }
  if (.min_lineage_sequences < 2) {
    warning(
      ".min_lineage_sequences is set to less than 2; ",
      "results will not be valid to build trees with repClonalLineage()!"
    )
  }

  parallel_prepare <- .prepare_threads > 1
  if (parallel_prepare) {
    doParallel::registerDoParallel(cores = .prepare_threads)
  }
  .data %<>%
    apply_to_sample_or_list(
      align_single_df,
      .min_lineage_sequences = .min_lineage_sequences,
      .parallel_prepare = parallel_prepare,
      .align_threads = .align_threads,
      .verbose_output = .verbose_output
    )
  if (parallel_prepare) {
    doParallel::stopImplicitCluster()
  }
  return(.data)
}

align_single_df <- function(data,
                            .min_lineage_sequences,
                            .parallel_prepare,
                            .align_threads,
                            .verbose_output) {
  for (required_column in c(
    "Cluster", "Germline.sequence", "V.germline.nt", "J.germline.nt", "CDR3.germline.length"
  )) {
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
      .parallel = .parallel_prepare
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
  alignments %<>% par_or_normal_lapply(
    align_sequences,
    .verbose_output = .verbose_output,
    mc.preschedule = TRUE,
    mc.cores = .align_threads
  )

  return(convert_results_to_df(results, alignments))
}

# this function accepts dataframe subset containing rows only for current lineage
# and returns named list containing 1 row for results dataframe
prepare_results_row <- function(lineage_subset, .min_lineage_sequences, .verbose_output) {
  cluster_name <- lineage_subset[[1, "Cluster"]]
  germline_seq <- lineage_subset[[1, "Germline.sequence"]]
  germline_v <- lineage_subset[[1, "V.germline.nt"]]
  germline_j <- lineage_subset[[1, "J.germline.nt"]]
  germline_cdr3_len <- lineage_subset[[1, "CDR3.germline.length"]]
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

  sequences_columns <- c(
    "Sequence", "Clone.ID", "Clones",
    "CDR1.nt", "CDR2.nt", "CDR3.nt", "FR1.nt", "FR2.nt", "FR3.nt", "FR4.nt"
  )
  if (.verbose_output) {
    sequences_columns %<>% c("V.end", "J.start", "CDR3.start", "CDR3.end")
  }
  sequences <- lineage_subset[sequences_columns]
  sequences[["Clone.ID"]] %<>% as.integer()
  sequences[["Clones"]] %<>% as.integer()

  germline_v_len <- str_length(germline_v)
  germline_j_len <- str_length(germline_j)
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

  clonotypes_names <- sapply(lineage_subset[["Clone.ID"]], function(id) {
    paste0("ID_", id)
  })
  all_sequences_list <- c(list(germline_trimmed), as.list(clonotypes_trimmed))
  names(all_sequences_list) <- c("Germline", clonotypes_names)
  alignment <- convert_seq_list_to_dnabin(all_sequences_list)

  # calculate amino acid frame starts for all trimmed sequences including germline
  germline_aa_start <- (germline_v_len - v_min_len) %% 3 + 1
  clonotypes_aa_starts <- (lineage_subset[["V.lengths"]] - v_min_len) %% 3 + 1
  all_sequences_aa_starts <- c(list(germline_aa_start), as.list(clonotypes_aa_starts))
  names(all_sequences_aa_starts) <- names(all_sequences_list)

  if (.verbose_output) {
    return(list(
      Cluster = cluster_name,
      Germline = germline_seq,
      V.germline.nt = germline_v,
      J.germline.nt = germline_j,
      CDR3.germline.length = germline_cdr3_len,
      Aligned = aligned,
      Alignment = alignment,
      V.length = v_min_len,
      J.length = j_min_len,
      Sequences = sequences,
      AA.frame.starts = all_sequences_aa_starts
    ))
  } else {
    return(list(
      Cluster = cluster_name,
      Germline = germline_seq,
      V.germline.nt = germline_v,
      J.germline.nt = germline_j,
      CDR3.germline.length = germline_cdr3_len,
      Alignment = alignment,
      V.length = v_min_len,
      J.length = j_min_len,
      Sequences = sequences,
      AA.frame.starts = all_sequences_aa_starts
    ))
  }
}

# trim V/J tails in sequence to the specified lenghts v_min, j_min
trim_seq <- function(seq, v_len, v_min, j_len, j_min) {
  str_sub(seq, v_len - v_min + 1, -(j_len - j_min + 1))
}

convert_results_to_df <- function(nested_results_list, nested_alignments_list) {
  alignments <- nested_alignments_list %>%
    lapply(magrittr::extract2, "Alignment") %>%
    tibble(Alignment = .)
  sequences <- nested_results_list %>%
    lapply(magrittr::extract2, "Sequences") %>%
    tibble(Sequences = .)
  frame_starts <- nested_results_list %>%
    lapply(magrittr::extract2, "AA.frame.starts") %>%
    tibble(AA.frame.starts = .)
  df <- nested_results_list %>%
    lapply(rlist::list.remove, c("Alignment", "Sequences", "AA.frame.starts")) %>%
    purrr::map_dfr(~.) %>%
    cbind(alignments, sequences, frame_starts)
  # fix column types after dataframe rebuilding
  for (column in c("CDR3.germline.length", "V.length", "J.length")) {
    df[[column]] %<>% as.integer()
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
