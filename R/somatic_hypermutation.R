#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr unnest
#' @importFrom plyr adply
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom ape as.DNAbin clustal
#'
#' @export repSomaticHypermutation
repSomaticHypermutation <- function(.data, .threads = parallel::detectCores(), .nofail = FALSE) {
  if (!require_system_package("clustalw", error_message = paste0(
    "repSomaticHypermutation requires Clustal W app to be installed!\n",
    "Please download it from here: http://www.clustal.org/download/current/\n",
    "or install it with your system package manager (such as apt or dnf)."
  ), .nofail, identical(.data, NA))) {
    return(NA)
  }

  doParallel::registerDoParallel(cores = .threads)
  results <- .data %>% apply_to_sample_or_list(
    shm_process_dataframe,
    .with_names = TRUE,
    .validate = FALSE
  )
  return(results)
}

shm_process_dataframe <- function(df, sample_name = NA) {
  # convert dataframe from cluster-per-row to clonotype-per-row
  # and add columns with alignment and number of mutations
  df %<>% unnest("Sequences") %>% plyr::adply(
    .fun = shm_process_clonotype_row,
    .margins = 1,
    .parallel = TRUE,
    sample_name = sample_name
  )
  return(df)
}

shm_process_clonotype_row <- function(row, sample_name) {
  v_gene_germline <- NA
  j_gene_germline <- NA
  v_gene_clonotype <- paste0(
    row[["FR1.nt"]], row[["CDR1.nt"]], row[["FR2.nt"]], row[["CDR2.nt"]], row[["FR3.nt"]]
  )
  j_gene_clonotype <- row[["FR4.nt"]]
  return(row)
}
