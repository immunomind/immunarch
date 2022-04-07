#' This function uses PHYLIP package to make phylogenetic analysis. For making trees it uses
#' maximum parsimony methods.
#'
#' @concept phylip
#'
#' @aliases repClonalFamily
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr discard imap
#' @importFrom utils capture.output
#' @importFrom parallel mclapply detectCores
#' @importFrom phangorn write.phyDat
#' @importFrom uuid UUIDgenerate

#' @description This function builds a phylogenetic tree using the sequences of a clonal lineage
#'
#' @usage
#'
#' repClonalFamily(.data, .threads)
#'
#' @param .data The data to be processed. Can be output of repAlignLineage, Alignment column
#' from repAlignLineage output (filtered by Aligned column or not), or list of samples with
#' Alignment column for each sample.
#'
#' @param .threads Number of threads to use.
#'
#' @return
#'
#' list of trees in "phylo" format and output files of PHYLIP dnapars function;
#' or list of these lists if the input was a list of samples.
#'
#' @examples
#'
#' data(bcrdata)
#' bcr_data <- bcrdata$data %>% top(100) # reduce the dataset to save time on examples
#'
#' bcr_data %>%
#'   seqCluster(seqDist(bcr_data), .fixed_threshold = 3) %>%
#'   repGermline() %>%
#'   repAlignLineage() %>%
#'   repLineagePhylogeny()
#' @export repClonalFamily
repClonalFamily <- function(.data, .threads = parallel::detectCores()) {
  require_system_package("phylip", error_message = paste0(
    "repLineagePhylogeny requires PHYLIP app to be installed!\n",
    "Please install it as described here:\n",
    "https://evolution.genetics.washington.edu/phylip/install.html"
  ))

  # verify input data, detect format and run process function that fits the format
  if (length(.data) == 0) {
    stop("repClonalFamily: input data is empty!")
  } else if (inherits(.data, "list")) {
    if (inherits(.data[[1]], "DNAbin")) {
      # if we got only Alignment column, filter out unaligned rows manually;
      # unaligned rows are DNAbin lists, so for them nrow() returns NULL
      results <- .data %>%
        purrr::discard(~ is.null(nrow(.x))) %>%
        process_list_of_alignments(.threads)
    } else if (inherits(.data[[1]], "data.frame")) {
      results <- .data %>%
        lapply(process_dataframe, .threads = .threads, sample_name = sample_name)
    } else if (inherits(.data[[1]], "list")) {
      # this is list of samples, and each has only Alignment column in it
      results <- .data %>%
        purrr::imap(function(sample_data, sample_name) {
          if (length(sample_data) == 0) {
            stop("repClonalFamily: input contains empty sample ", sample_name)
          }
          if (inherits(sample_data[[1]], "DNAbin")) {
            sample_results <- sample_data %>%
              purrr::discard(~ is.null(nrow(.x))) %>%
              process_list_of_alignments(.threads, sample_name = sample_name)
          } else {
            stop("Unrecognized format of input data for repClonalFamily in sample ", sample_name)
          }
          return(sample_results)
        })
    } else {
      stop("Unrecognized format of input data for repClonalFamily!")
    }
  } else if (inherits(.data, "data.frame")) {
    results <- .data %>%
      process_dataframe(.threads)
  } else {
    stop("Unrecognized format of input data for repClonalFamily!")
  }
  return(results)
}

process_dataframe <- function(df, .threads, sample_name = NA) {
  for (column in c("Aligned", "Alignment")) {
    if (!(column %in% colnames(df))) {
      stop(
        "Unrecognized input dataframe format for repClonalFamily: missing \"",
        column,
        "\" column",
        optional_sample(" in sample ", sample_name, ""),
        "!"
      )
    }
  }

  results <- df %>%
    filter(Aligned) %>%
    magrittr::extract2("Alignment") %>%
    process_list_of_alignments(.threads, sample_name)
  return(results)
}

process_list_of_alignments <- function(alignments, .threads, sample_name = NA) {
  if (length(alignments) == 0) {
    stop(
      "repClonalFamily: aligned sequences not found in input data",
      optional_sample(" in sample ", sample_name, ""),
      "!"
    )
  }

  results <- parallel::mclapply(
    alignments,
    process_alignment,
    mc.preschedule = FALSE,
    mc.cores = .threads
  )
  return(results)
}

process_alignment <- function(alignment) {
  temp_dir <- file.path(tempdir(check = TRUE), uuid::UUIDgenerate(use.time = FALSE))
  dir.create(temp_dir)
  # bugfix for phylip: it shows "Unexpected end-of-file" for too short sequence labels
  rownames(alignment) %<>% paste0("\t")
  phangorn::write.phyDat(alignment, file.path(temp_dir, "infile"))
  system(
    paste0("sh -c \"cd ", temp_dir, "; phylip dnapars infile\""),
    input = (c("V", 1, 5, "Y"))
  ) %>%
    quiet()
  outfile <- read_file(file.path(temp_dir, "outfile"))
  outtree <- read_file(file.path(temp_dir, "outtree"))
  unlink(temp_dir, recursive = TRUE)
  return(list(OUTFILE = outfile, OUTTREE = outtree))
}
