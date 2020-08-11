#' Track clonotypes across time and data points
#'
#' @concept dynamics
#'
#' @importFrom dplyr arrange desc distinct
#' @importFrom tibble tibble
#'
#' @aliases trackClonotypes
#'
#' @description track top abundant clonotypes across timepoints
#'
#' @param .data The data to process. It can be a \link{data.frame}, a
#' \link{data.table}, or a list of these objects.
#'
#' Every object must have columns in the immunarch compatible format.
#' \link{immunarch_data_format}
#'
#' Competent users may provide advanced data representations:
#' DBI database connections, Apache Spark DataFrame from \link{copy_to} or a list
#' of these objects. They are supported with the same limitations as basic objects.
#'
#' Note: each connection must represent a separate repertoire.
#'
#' @param .timepoint An argument that regulates which clonotypes to choose for tracking. There are three options for this argument:
#'
#' 1) pass a list with two elements \code{list(X, Y)}, where \code{X} is the name or the index of a target repertoire from ".data", and
#' \code{Y} is the number of the most abundant clonotypes to take from \code{X}.
#'
#' 2) pass a character vector of sequences to take from all data frames;
#'
#' 3) pass a data frame (data table, database) with one or more columns - first for sequences, and other for gene segments (if applicable).
#'
#' See the "Examples" below with examples for each option.
#'
#' @param .col A character vector of length 1. Specifies an identifier for a column, from which the function
#' chooses clonotype sequences. Specify "nt" for nucleotide sequences, "aa" for amino acid sequences,
#' "aa+v" for amino acid sequences and Variable genes, "nt+j" for nucleotide
#' sequences with Joining genes, or any combination of the above.
#' Used only if ".timepoint" has option 1) or option 2).
#'
#' @param .norm Logical. If TRUE then use Proportion instead of the number of Clones per clonotype to store
#' in the function output.
#'
#' @description
#' Note: duplicated clonotypes are merged and their counts are summed up.
#'
#' @return Data frame with input sequences and counts or proportions for each of the input repertoire.
#' @export trackClonotypes
#' 

trackClonotypes2 <- function(.data, .meta, .timepoints = list(), .num, .col = "aa", .norm = TRUE) {
  # subsets data using timepoint field
  # top n most abundant clonotypes
  # run pca
  if (!hasArg(.timepoints)) {
    stop("Error: please pass a list of timepoints to the function.")
  }

  if (!has_class(.data, "list")) {
    stop("Error: please pass a list with immune repertoires to track clonotypes.")
  }
  if (length(.data) < 2) {
    stop("Error: please pass a list with 2 or more immune repertoires to track clonotypes.")
  }

  .col <- unlist(strsplit(.col, split = "\\+"))
  .col <- sapply(.col, switch_type, USE.NAMES = FALSE)

  if (length(.timepoints) != 2) {
    stop("Error: please pass a list with two elements for the .timepoint argument. Run ?trackClonotypes for more details in the documentation.")
  }

  # n_clonotypes assigned by .num parameter
  n_clonotypes <- .num

  subset_df <- .data[.timepoints]

  # for each sample in the list in .timepoint, get top .num clonotypes
  topdf <- lapply(subset_df, function(w) { w$top <- trackClonotypes(w, list(1, .num), .col, .norm) })

  temp <- data.frame()

  for (sample in topdf) {
    temp <- temp %>% full_join(sample$top)
  }

 # @TODO fix .var so it works for all columns
  temp <- temp %>% column_to_rownames(., var = "CDR3.aa")

  m_raw <- data.matrix(temp)

  #### NORMALISATION
  # adding column with max proportion for each clonotype across time points
  m <- m_raw %>% rowwise() %>% mutate(max_prop = max(across(test)))

  # dividing proportions by the max proportion across all time points
  m <- m %>%
  ungroup() %>%
  mutate(across(test), ~ . / max_prop))

  # removing max prop column
  m$max_prop <- NULL

  m <- as.data.frame(m)

  rownames(m) <- rownames(m.raw)

  result_df$pca <- immunr_pca(m)

  ### GENERATE TRAJECTORIES
  
  library(factoextra)

  res <- hcut(as.matrix(m), k = 3, stand = FALSE) 
  
  m$clust <- res$cluster

  result_df$data <- m

  m <- m %>% reshape2::melt(id.vars="clust")

  result_df$traj <- m %>% group_by(clust, day) %>% summarise(mean_traj = mean(value), n_traj = n(), SE = sd(value)/sqrt(n()), SE_scaled = 2.96*sd(value)/sqrt(n()))

  add_class(result_df, "immunr_trajectories")


}