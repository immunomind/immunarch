#' This function aligns all sequences incliding germline within one clonal lineage. After building
#' clonal lineage and germline, the next step isto analyze the degree of mutation and maturity
#' of each clonal lineage. This allows you to find high mature cells and cells with a large number of offspring.
#' The phylogenetic analysis will find mutations that increase the affinity of BCR.
#' Making alignment of the sequence is the first step towards sequence analysis including BCR.
#'
#' @concept align_lineage
#'
#' @aliases repAlignLineage
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom ape as.DNAbin muscle

#' @description Aligns all sequences incliding germline within one clonal lineage
#'
#' @usage
#'
#' repAlignLineage(.data)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}, and also
#' must contain Sequence.Germline column, which is added by repGermline() function.
#'
#' @return
#'
#' DNAbin object with alignment if there was single input dataframe,
#' or list of DNAbin objects if the input was a list of dataframes
#'
#' @examples
#'
#' data(immdata)
#'
#' immdata$data %>%
#'   top(2000) %>% # reduce the dataset to save time on examples
#'   repGermline() %>%
#'   repAlignLineage()
#' @export repAlignLineage
repAlignLineage <- function(.data) {
  require_system_package("muscle", error_message = paste0(
    "repAlignLineage requires MUSCLE app to be installed!\n",
    "Please download it from here: https://github.com/rcedgar/muscle/releases/latest\n",
    "or install it with your system package manager (such as apt or dnf)."
  ))
  if (inherits(.data, "list")) {
    .validate_repertoires_data(.data)
    .data %>%
      lapply(function(sample_data) {
        sample_data %>%
          as_tibble() %>%
          align_single_df() %>%
          return()
      }) %>%
      return()
  } else {
    .data %>%
      as_tibble() %>%
      align_single_df() %>%
      return()
  }
}

align_single_df <- function(data) {
  if (!("Germline.Sequence" %in% colnames(data))) {
    stop(
      "Found dataframe without required column Germline.Sequence;\n",
      "existing columns: ", str(colnames(data))
    )
  }

  # fastaLines = c(as.character(paste(">", str_pad('germline', 10, 'right', ' '), sep = '')))
  # fastaLines = c(fastaLines, as.character(data[1, 'germline_sequence']))
  # for (rowNum in 1:nrow(data)){
  #   fastaLines = c(fastaLines, as.character(paste(">", str_pad(paste(paste(data[rowNum, 'group_id'], data[rowNum, 'CLONE_LINEAGE_ID'], sep = '_'), rowNum, sep = '_'), 10, 'right', ' '), sep = "")))
  #   fastaLines = c(fastaLines,as.character(data[rowNum,"SEQUENCE_INPUT"]))
  # }

  seq_list <- list()
  seq_list["germline"] <- data[1, "Germline.Sequence"]

  seq_list %>%
    ape::as.DNAbin() %>%
    ape::muscle() %>%
    return()
}
