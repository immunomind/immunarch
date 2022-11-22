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
#' repAlignLineage(.data, .min_lineage_sequences, .prepare_threads, .align_threads, .nofail)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' @param .min_lineage_sequences If number of sequences in the same clonal lineage and the same
#' cluster (not including germline) is lower than this threshold, this group of sequences
#' will be filtered out from the dataframe; so only large enough lineages will be included.
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
#' @param .nofail Will return NA instead of stopping if Clustal W is not installed.
#' Used to avoid raising errors in examples on computers where Clustal W is not installed.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has these columns:
#' * Cluster: cluster name
#' * Germline: germline sequence
#' * Alignment: DNAbin object with alignment
#' * Sequences: nested dataframe containing all sequences for this combination
#'   of cluster and germline; it has columns
#'   * Sequence, CDR1.nt, CDR2.nt, CDR3.nt, FR1.nt, FR2.nt, FR3.nt, FR4.nt, V.allele, J.allele,
#'     V.aa, J.aa: all values taken from the input dataframe
#'   * Clone.ID: taken from the input dataframe, or created (filled with row numbers) if missing
#'   * Clones: taken from the input dataframe, or created (filled with '1' values) if missing
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
                            .nofail = FALSE) {
  if (!require_system_package(c("clustalw", "clustalw2"), error_message = paste0(
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
      .align_threads = .align_threads
    )
  if (parallel_prepare) {
    doParallel::stopImplicitCluster()
  }
  return(.data)
}

align_single_df <- function(data,
                            .min_lineage_sequences,
                            .parallel_prepare,
                            .align_threads) {
  for (required_column in c(
    "Cluster", "Germline.sequence", "V.allele", "J.allele",
    "FR1.nt", "CDR1.nt", "FR2.nt", "CDR2.nt", "FR3.nt", "CDR3.nt", "FR4.nt", "V.aa", "J.aa"
  )) {
    if (!(required_column %in% colnames(data))) {
      stop(
        "Found dataframe without required column ",
        required_column,
        ";\nexisting columns: ",
        toString(colnames(data))
      )
    }
  }

  results <- data %>%
    fill_missing_columns() %>%
    plyr::dlply(
      .variables = .(get("Cluster"), get("Germline.sequence")),
      .fun = prepare_results_row,
      .min_lineage_sequences = .min_lineage_sequences,
      .parallel = .parallel_prepare
    ) %>%
    `[`(!is.na(.)) %>%
    unname()

  if (length(results) == 0) {
    stop("There are no lineages containing at least ", .min_lineage_sequences, " sequences!")
  }

  # only Alignment column are passed to alignment function to reduce consumed memory
  alignments <- lapply(results, "[", "Alignment") %>%
    par_or_normal_lapply(mc.preschedule = TRUE, mc.cores = .align_threads, function(df_row) {
      df_row[["Alignment"]] %<>% ape::clustal()
    })

  return(convert_results_to_df(results, alignments))
}

# fill Clone.ID and Clones columns if they are missing
fill_missing_columns <- function(data) {
  if (!("Clone.ID" %in% colnames(data))) {
    data[["Clone.ID"]] <- seq.int(nrow(data))
  }
  if (!("Clones" %in% colnames(data))) {
    data[["Clones"]] <- as.integer(1)
  }
  return(data)
}

# this function accepts dataframe subset containing rows only for current lineage
# and returns named list containing 1 row for results dataframe
prepare_results_row <- function(lineage_subset, .min_lineage_sequences) {
  if (nrow(lineage_subset) < .min_lineage_sequences) {
    # NA rows will be filtered out
    return(NA)
  }

  cluster_name <- lineage_subset[[1, "Cluster"]]
  germline_seq <- lineage_subset[[1, "Germline.sequence"]]

  sequences_columns <- c(
    "Sequence", "Clone.ID", "Clones", "V.allele", "J.allele",
    "CDR1.nt", "CDR2.nt", "CDR3.nt", "FR1.nt", "FR2.nt", "FR3.nt", "FR4.nt", "V.aa", "J.aa"
  )

  sequences <- lineage_subset[sequences_columns]
  sequences[["Clone.ID"]] %<>% as.integer()
  sequences[["Clones"]] %<>% as.integer()

  clonotypes_names <- sapply(lineage_subset[["Clone.ID"]], function(id) {
    paste0("ID_", id)
  })
  all_sequences_list <- c(list(germline_seq), as.list(lineage_subset[["Sequence"]]))
  names(all_sequences_list) <- c("Germline", clonotypes_names)
  alignment <- convert_seq_list_to_dnabin(all_sequences_list)

  return(list(
    Cluster = cluster_name,
    Germline = germline_seq,
    Alignment = alignment,
    Sequences = sequences
  ))
}

convert_results_to_df <- function(nested_results_list, alignments_list) {
  alignments <- tibble(Alignment = alignments_list)
  sequences <- nested_results_list %>%
    lapply(magrittr::extract2, "Sequences") %>%
    tibble(Sequences = .)
  df <- nested_results_list %>%
    lapply(rlist::list.remove, c("Alignment", "Sequences")) %>%
    purrr::map_dfr(~.) %>%
    cbind(alignments, sequences)
  return(df)
}
