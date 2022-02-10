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
#' @importFrom ape as.DNAbin muscle
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
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}, and also
#' must contain Sequence.germline column, which is added by repGermline() function.
#'
#' @return
#'
#' DNAbin object with alignment if there was single input dataframe,
#' or list of DNAbin objects if the input was a list of dataframes
#'
#' @examples
#'
#' data(bcrdata)
#'
#' bcrdata$data %>%
#'   top(100) %>% # reduce the dataset to save time on examples
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
    .data %<>%
      lapply(function(sample_data) {
        sample_data %>%
          as_tibble() %>%
          align_single_df()
      })
    return(.data)
  } else {
    .data %<>%
      as_tibble() %>%
      align_single_df()
    return(.data)
  }
}

align_single_df <- function(data) {
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
    ape::muscle()
}
