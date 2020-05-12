#' Track clonotypes across time and data points
#'
#' @concept dynamics
#'
#' @importFrom dplyr arrange desc distinct
#' @importFrom tibble tibble
#'
#' @aliases trackClonotypes
#'
#' @description
#' Track the temporal dynamics of clonotypes in repertoires. For example, tracking across multiple
#' time points after vaccination.
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
#' @param .which An argument that regulates which clonotypes to choose for tracking. There are three options for this argument:
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
#' Used only if ".which" has option 1) or option 2).
#'
#' @param .norm Logical. If TRUE then use Proportion instead of the number of Clones per clonotype to store
#' in the function output.
#'
#' @description
#' Note: duplicated clonotypes are merged and their counts are summed up.
#'
#' @return Data frame with input sequences and counts or proportions for each of the input repertoire.
#'
#' @examples
#' # Load an example data that comes with immunarch
#' data(immdata)
#'
#' # Make the data smaller in order to speed up the examples
#' immdata$data <- immdata$data[c(1, 2, 3, 7, 8, 9)]
#' immdata$meta <- immdata$meta[c(1, 2, 3, 7, 8, 9), ]
#'
#' # Option 1
#' # Choose the first 10 amino acid clonotype sequences
#' # from the first repertoire to track
#' tc <- trackClonotypes(immdata$data, list(1, 10), .col = "aa")
#' # Choose the first 20 nucleotide clonotype sequences
#' # and their V genes from the "MS1" repertoire to track
#' tc <- trackClonotypes(immdata$data, list("MS1", 20), .col = "nt+v")
#'
#' # Option 2
#' # Choose clonotypes with amino acid sequences "CASRGLITDTQYF" or "CSASRGSPNEQYF"
#' tc <- trackClonotypes(immdata$data, c("CASRGLITDTQYF", "CSASRGSPNEQYF"), .col = "aa")
#'
#' # Option 3
#' # Choose the first 10 clonotypes from the first repertoire
#' # with amino acid sequences and V segments
#' target <- immdata$data[[1]] %>%
#'   select(CDR3.aa, V.name) %>%
#'   head(10)
#' tc <- trackClonotypes(immdata$data, target)
#'
#' # Visualise the output regardless of the chosen option
#' # Therea are three way to visualise it, regulated by the .plot argument
#' vis(tc, .plot = "smooth")
#' vis(tc, .plot = "area")
#' vis(tc, .plot = "line")
#'
#' # Visualising timepoints
#' # First, we create an additional column in the metadata with randomly choosen timepoints:
#' immdata$meta$Timepoint <- sample(1:length(immdata$data))
#' immdata$meta
#' # Next, we create a vector with samples in the right order,
#' # according to the "Timepoint" column (from smallest to greatest):
#' sample_order <- order(immdata$meta$Timepoint)
#' # Sanity check: timepoints are following the right order:
#' immdata$meta$Timepoint[sample_order]
#' # Samples, sorted by the timepoints:
#' immdata$meta$Sample[sample_order]
#' # And finally, we visualise the data:
#' vis(tc, .order = sample_order)
#' @export trackClonotypes
trackClonotypes <- function(.data, .which = list(1, 15), .col = "aa", .norm = TRUE) {
  if (!has_class(.data, "list")) {
    stop("Error: please pass a list with immune repertoires to track clonotypes.")
  }
  if (length(.data) < 2) {
    stop("Error: please pass a list with 2 or more immune repertoires to track clonotypes.")
  }

  .col <- unlist(strsplit(.col, split = "\\+"))
  .col <- sapply(.col, switch_type, USE.NAMES = FALSE)

  if (has_class(.which, "list")) {
    # Option 1
    if (length(.which) != 2) {
      stop("Error: please pass a list with two elements for the .which argument. Run ?trackClonotypes for more details in the documentation.")
    }

    target_df <- .which[[1]]
    n_clonotypes <- .which[[2]]

    # ToDo: replace Clones with IMMCOL$count
    target_seq <- .data[[target_df]] %>%
      arrange(desc(Clones)) %>%
      select(.col) %>%
      head(n_clonotypes) %>%
      distinct() %>%
      collect(n = Inf)
  } else if (has_class(.which, "character")) {
    # Option 2
    target_seq <- tibble(.which) %>% distinct()

    # ToDo: extract sequence column names
    names(target_seq)[1] <- .col[1]
  } else {
    # Option 3
    target_seq <- .which %>% distinct()
    .col <- colnames(target_seq)

    for (col_name in .col) {
      if (!(col_name %in% names(.data[[1]]))) {
        stop("Error: can't find a column with the name \"", col_name, "\" in the input data. Please check column names for the .which argument.")
      }
    }
  }

  setDT(target_seq)
  count_col <- IMMCOL$count

  result_df <- NULL
  for (i_df in 1:length(.data)) {
    temp_df <- .data[[i_df]] %>%
      select(.col, Count = count_col)
    setDT(temp_df)

    if (.norm) {
      temp_df$Count <- temp_df$Count / sum(temp_df$Count)
    }

    temp_df <- temp_df[, .(sum(Count)), .col]

    if (is.null(result_df)) {
      result_df <- merge(target_seq, temp_df, all.x = TRUE)
    } else {
      result_df <- merge(result_df, temp_df, all.x = TRUE)
    }

    setnames(result_df, ncol(result_df), names(.data)[i_df])
  }

  for (j in seq_len(ncol(result_df))) {
    set(result_df, which(is.na(result_df[[j]])), j, 0)
  }

  add_class(result_df, "immunr_dynamics")
}


# poisson_model <- function (.data, .q = c(.025, .975)) {
#   col = collect(select(.data, Count))[[1]]
#   lo = qpois(.q[1], col)
#   hi = qpois(.q[2], col)
#   res = cbind(Data = col, Lo = lo, Hi = hi)
#   add_class(res, "immunr_dynamics")
# }
#
# norm_model <- function (.data, .q = c(.025, .975), .sd = 1.05) {
#   col = collect(select(.data, Count))[[1]]
#   lo = qnorm(.q[1], col, .sd)
#   hi = qnorm(.q[2], col, .sd)
#   res = cbind(Data = col, Lo = lo, Hi = hi)
#   add_class(res, "immunr_dynamics")
# }
