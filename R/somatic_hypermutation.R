#' Calculates number of mutations against the germline for each clonotype
#'
#' @concept somatic_hypermutation
#'
#' @aliases repSomaticHypermutation
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr unnest
#' @importFrom plyr adply
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom ape as.DNAbin clustal
#'
#' @description This function aligns V and J genes from the germline in each cluster
#' with corresponding genes in each clonotype, saves the alignments for purpose of visualization,
#' and calculates number of mutations for each clonotype.
#'
#' @usage
#'
#' repSomaticHypermutation(.data, .threads, .nofail)
#'
#' @param .data The data to be processed: an output of repClonalFamily();
#' variants with one sample and list of samples are both supported.
#'
#' @param .threads Number of threads to use.
#'
#' @param .nofail Will return NA instead of stopping if Clustal W is not installed.
#' Used to avoid raising errors in examples on computers where Clustal W is not installed.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has all the columns from repClonalFamily() output dataframe, with Sequence
#' column unnested: the resulting dataframe has one line per clonotype. Clone.ID column
#' contains original IDs for clonotypes, and can be used as dataframe key.
#' New columns are added:
#' * Germline.Alignment.V: contains V gene alignment of current clonotype with the germline
#' * Germline.Alignment.J: contains J gene alignment of current clonotype with the germline
#' * Substitutions: contains number of substitutions in the alignment (summary for V and J)
#' * Insertions: contains number of insertions in the clonotype relative to germline
#'   (summary for V and J)
#' * Deletions: contains number of deletions in the clonotype relative to germline
#'   (summary for V and J)
#' * Mutations: contains total number of mutations in the alignment (summary for V and J)
#'
#' @examples
#'
#' data(bcrdata)
#' bcr_data <- bcrdata$data
#'
#' bcr_data %>%
#'   seqCluster(seqDist(bcr_data), .fixed_threshold = 3) %>%
#'   repGermline(.threads = 2) %>%
#'   repAlignLineage(.min_lineage_sequences = 2, .align_threads = 2, .nofail = TRUE) %>%
#'   repClonalFamily(.threads = 2, .nofail = TRUE) %>%
#'   repSomaticHypermutation(.threads = 2, .nofail = TRUE)
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
    .parallel = .threads > 1,
    .validate = FALSE
  )
  doParallel::stopImplicitCluster()
  return(results)
}

shm_process_dataframe <- function(df, .parallel) {
  # convert dataframe from cluster-per-row to clonotype-per-row
  # and add columns with alignment and number of mutations
  df %<>% unnest("Sequences") %>% plyr::adply(
    .fun = shm_process_clonotype_row,
    .margins = 1,
    .parallel = .parallel
  )
  # fix column types after dataframe rebuilding
  for (column in c(
    "Clone.ID", "Clones", "CDR3.germline.length", "Trunk.Length",
    "Substitutions", "Insertions", "Deletions", "Mutations"
  )) {
    df[[column]] %<>% as.integer()
  }
  return(df)
}

shm_process_clonotype_row <- function(row) {
  v_gene_germline <- row[["V.germline.nt"]]
  j_gene_germline <- row[["J.germline.nt"]]
  v_gene_clonotype <- paste0(
    row[["FR1.nt"]], row[["CDR1.nt"]], row[["FR2.nt"]], row[["CDR2.nt"]], row[["FR3.nt"]]
  )
  j_gene_clonotype <- row[["FR4.nt"]]

  alignments <- list(
    V = list(v_gene_germline, v_gene_clonotype),
    J = list(j_gene_germline, j_gene_clonotype)
  )

  substitutions <- 0
  insertions <- 0
  deletions <- 0
  for (gene in c("V", "J")) {
    names(alignments[[gene]]) <- c("Germline", paste0("ID_", row[["Clone.ID"]]))
    dnabin <- convert_seq_list_to_dnabin(alignments[[gene]])

    # align gene from germline and gene from clonotype
    alignments[[gene]] <- ape::clustal(dnabin)

    # calculate numbers of mutations for current gene, add them to the sums
    alignment <- as.character(alignments[[gene]])
    germline <- alignment[1, ]
    clonotype <- alignment[2, ]
    substitutions <- substitutions +
      sum((germline != clonotype) & (germline != "-") & (clonotype != "-"))
    insertions <- insertions + sum(germline == "-")
    deletions <- deletions + sum(clonotype == "-")
  }

  # using list(NA) as workaround for "Assigned data must be compatible with existing data":
  # writing list(NA), which is considered compatible by adply, then storing alignment in list
  row[["Germline.Alignment.V"]] <- list(NA)
  row[["Germline.Alignment.V"]][[1]] <- alignments[["V"]]
  row[["Germline.Alignment.J"]] <- list(NA)
  row[["Germline.Alignment.J"]][[1]] <- alignments[["J"]]
  row[["Substitutions"]] <- substitutions
  row[["Insertions"]] <- insertions
  row[["Deletions"]] <- deletions
  row[["Mutations"]] <- substitutions + insertions + deletions

  return(row)
}
