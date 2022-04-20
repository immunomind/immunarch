#' This function uses PHYLIP package to make phylogenetic analysis. For making trees it uses
#' maximum parsimony methods.
#'
#' @concept phylip
#'
#' @aliases repClonalFamily
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map_dfr
#' @importFrom utils capture.output
#' @importFrom parallel mclapply detectCores
#' @importFrom phangorn write.phyDat
#' @importFrom ape read.tree
#' @importFrom uuid UUIDgenerate
#' @importFrom data.table fread

#' @description This function builds a phylogenetic tree using the sequences of a clonal lineage
#'
#' @usage
#'
#' repClonalFamily(.data, .threads)
#'
#' @param .data The data to be processed. Can be output of repAlignLineage() with normal
#' or verbose output; variants with one sample and list of samples are both supported.
#'
#' @param .threads Number of threads to use.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has these columns:
#' * Cluster: cluster name
#' * Germline.Input: germline sequence, like it was in the input
#' * Germline.Output: germline sequence, parsed from PHYLIP dnapars function output
#' * Common.Ancestor: common ancestor sequence, parsed from PHYLIP dnapars function output
#' * Trunk.Length: mean trunk length, representing the distance between the most recent
#'   common ancestor and germline sequence as a measure of the maturity of a lineage
#' * Tree: output tree in "phylo" format, loaded from by PHYLIP dnapars function output
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
#'   repClonalFamily()
#' @export repClonalFamily
repClonalFamily <- function(.data, .threads = parallel::detectCores()) {
  require_system_package("phylip", error_message = paste0(
    "repLineagePhylogeny requires PHYLIP app to be installed!\n",
    "Please install it as described here:\n",
    "https://evolution.genetics.washington.edu/phylip/install.html"
  ))
  results <- .data %>% apply_to_sample_or_list(
    process_dataframe,
    .with_names = TRUE,
    .validate = FALSE,
    .threads = .threads
  )
  return(results)
}

process_dataframe <- function(df, .threads, sample_name = NA) {
  if (nrow(df) == 0) {
    stop(
      "repClonalFamily: input dataframe is empty",
      optional_sample(" in sample ", sample_name, ""),
      "!"
    )
  }
  if (!("Alignment" %in% colnames(df))) {
    stop(
      "Unrecognized input dataframe format for repClonalFamily: missing \"Alignment\" column",
      optional_sample(" in sample ", sample_name, ""),
      "!"
    )
  }
  if ("Aligned" %in% colnames(df)) {
    df %<>% filter(Aligned)
  }
  df <- df[c("Cluster", "Germline", "Alignment")]
  clusters_list <- split(df, seq(nrow(df)))

  results <- parallel::mclapply(
    clusters_list,
    process_cluster,
    mc.preschedule = FALSE,
    mc.cores = .threads
  ) %>%
    purrr::map_dfr(~.)
  return(results)
}

process_cluster <- function(cluster_row) {
  cluster_name <- cluster_row[["Cluster"]]
  cluster_germline <- cluster_row[["Germline"]]
  alignment <- cluster_row[["Alignment"]]

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

  tree <- ape::read.tree(file.path(temp_dir, "outtree"))
  outfile_path <- file.path(temp_dir, "outfile")
  common_ancestor_path <- file.path(temp_dir, "common_ancestor.csv")
  germline_path <- file.path(temp_dir, "germline.csv")

  # parse common ancestor
  system(paste0(
    "grep  '         1   ' ",
    outfile_path,
    " > ",
    common_ancestor_path
  ))
  system(paste0("cat ", common_ancestor_path))
  common_ancestor <- data.table::fread(common_ancestor_path, drop = c(1)) %>%
    t() %>%
    paste(collapse = "")

  # parse germline
  system(paste0(
    "grep  '1   germline' ",
    outfile_path,
    " > ",
    germline_path
  ))
  system(paste0("cat ", germline_path))
  germline <- read.delim(germline_path, delim)
  print(germline)
  germline <- data.table::fread(germline_path, drop = c(1:3)) %>%
    t() %>%
    paste(collapse = "") %>%
    as.character()

  unlink(temp_dir, recursive = TRUE)

  # return row of output dataframe as named list
  return(list(
    Cluster = cluster,
    Germline.Input = cluster_germline,
    Germline.Output = germline,
    Common.Ancestor = common_ancestor,
    Trunk.Length = trunk_length,
    Tree = tree
  ))
}
