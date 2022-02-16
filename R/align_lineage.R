#' This function aligns all sequences incliding germline within one clonal lineage. After building
#' clonal lineage and germline, the next step is to analyze the degree of mutation and maturity
#' of each clonal lineage. This allows you to find high mature cells and cells with a large number of offspring.
#' The phylogenetic analysis will find mutations that increase the affinity of BCR.
#' Making alignment of the sequence is the first step towards sequence analysis including BCR.
#'
#' @concept align_lineage
#'
#' @aliases repAlignLineage
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom stringr str_extract_all
#' @importFrom ape as.DNAbin clustal
#' @importFrom parallel mclapply detectCores

#' @description Aligns all sequences incliding germline within one clonal lineage
#'
#' @usage
#'
#' repAlignLineage(.data)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' @param .min.lineage.sequences If number of sequences in the clonal lineage
#' (not including germline) is lower than this threshold, this lineage will not be aligned and
#' will not be used in next steps of BCR pipeline (but still will be saved in the output table).
#'
#' @param .germline.tails Extra germline V/J length relative to the longest V/J gene
#' in it's lineage (counted separately for V and J genes).
#' Value 0 means that gene in the germline after trimming will have the same length as
#' the longest gene in it's lineage; default value 0.1 means that gene in the germline will be
#' 10% longer.
#'
#' @param .pw.gapopen Gap opening penalty used by Clustal during pairwise alignments.
#'
#' @param .pw.gapext Gap extension penalty used by Clustal during pairwise alignments.
#'
#' @param .gapopen Gap opening penalty used by Clustal during global alignments.
#'
#' @param .gapext Gap extension penalty used by Clustal during global alignments.
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}, and also
#' must contain Sequence.germline column, which is added by repGermline() function.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has number of rows equal to number of unique germlines, and has these columns:
#' * Germline: germline sequence
#' * Aligned: FALSE if it was not aligned with lineage (.min.lineage.sequences is below
#'   the threshold); TRUE if it was aligned
#' * Alignment: DNAbin object with alignment (if Aligned=TRUE) or DNAbin object with unaligned
#'   sequences (if Aligned=FALSE)
#' * V.trim: start (starting from 1, inclusive) of the germline part that is included
#'   in the alignment; all positions on the left of V.trim are present in Germline column,
#'   but not present in Alignment object
#' * J.trim: end (starting from 1, exclusive) of the germline part that is included
#'   in the alignment; all positions on J.trim position and on the right are present
#'   in Germline column, but not present in Alignment object
#' * Sequences: nested dataframe containing all sequences for this germline; has these columns:
#'   * Sequence: original sequence; V/J genes padded with nucleotides taken from
#'     the germline if they are shorter than the longest V/J genes in this lineage
#'     (only in sequences, not in the germline)
#'   * V.start: start (starting from 1, inclusive) of the original V gene of this sequence:
#'     this number shows where it's located in the padded sequence (that is in Sequence column)
#'   * V.end: end (starting from 1, exclusive) of the original V gene of this sequence:
#'     this number shows where it's located in the padded sequence (that is in Sequence column)
#'   * J.start: start (starting from 1, inclusive) of the original J gene of this sequence:
#'     this number shows where it's located in the padded sequence (that is in Sequence column)
#'   * J.end: end (starting from 1, exclusive) of the original J gene of this sequence:
#'     this number shows where it's located in the padded sequence (that is in Sequence column)
#' @examples
#'
#' data(bcrdata)
#'
#' bcrdata$data %>%
#'   top(100) %>% # reduce the dataset to save time on examples
#'   repGermline() %>%
#'   repAlignLineage()
#' @export repAlignLineage
repAlignLineage <- function(.data, .min.lineage.sequences = 3, .germline.tails = 0.1, .pw.gapopen = NA, .pw.gapext = NA, .gapopen = NA, .gapext = NA) {
  require_system_package("clustalw", error_message = paste0(
    "repAlignLineage requires Clustal W app to be installed!\n",
    "Please download it from here: http://www.clustal.org/download/current/\n",
    "or install it with your system package manager (such as apt or dnf)."
  ))
  if (inherits(.data, "list")) {
    .validate_repertoires_data(.data)
    .data %<>%
      lapply(function(sample_data) {
        sample_data %>%
          as_tibble() %>%
          align_single_df(
            .min.lineage.sequences, .germline.tails, .pw.gapopen, .pw.gapext, .gapopen, .gapext
          )
      })
    return(.data)
  } else {
    .data %<>%
      as_tibble() %>%
      align_single_df(
        .min.lineage.sequences, .germline.tails, .pw.gapopen, .pw.gapext, .gapopen, .gapext
      )
    return(.data)
  }
}

align_single_df <- function(data, .min.lineage.sequences, .germline.tails, .pw.gapopen, .pw.gapext, .gapopen, .gapext) {
  if (!("Germline.sequence" %in% colnames(data))) {
    stop(
      "Found dataframe without required column Germline.sequence;\n",
      "existing columns: ", str(colnames(data))
    )
  }

  germlines <- unique(data[["Germline.sequence"]])

  germlines %>%
    lapply(get_germline_with_lineage, data = data) %>%
    parallel::mclapply(align_sequences,
      mc.preschedule = FALSE, mc.cores = parallel::detectCores()
    )
}

# return a list containing the germline and all sequences with this germline
get_germline_with_lineage <- function(germline, data) {
  c(list(germline = germline), as.list(
    subset(data, Germline.sequence == germline)[["Sequence"]]
  ))
}

align_sequences <- function(list_of_sequences) {
  list_of_sequences %>%
    lapply(function(sequence) {
      sequence %>%
        stringr::str_extract_all(boundary("character")) %>%
        unlist()
    }) %>%
    ape::as.DNAbin() %>%
    ape::clustal()
}
