#' This function uses PHYLIP package to make phylogenetic analysis. For making trees it uses
#' maximum parsimony methods.
#'
#' @concept phylip
#'
#' @aliases repClonalFamily
#'
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom purrr map_dfr
#' @importFrom rlist list.remove
#' @importFrom stringr str_match str_count fixed
#' @importFrom stringi stri_replace_all_fixed
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
#' repClonalFamily(.data, .threads, .nofail)
#'
#' @param .data The data to be processed. Can be output of repAlignLineage() with normal
#' or verbose output; variants with one sample and list of samples are both supported.
#'
#' @param .threads Number of threads to use.
#'
#' @param .nofail Return NA instead of stopping if PHYLIP is not installed.
#' Used to avoid raising errors in examples on computers where PHYLIP is not installed.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has these columns:
#' * Cluster: cluster name
#' * Germline.Input: germline sequence, like it was in the input; not trimmed and not aligned
#' * Germline.Output: germline sequence, parsed from PHYLIP dnapars function output;
#'   it contains difference of germline from the common ancestor; "." characters mean
#'   matching letters. It's usually shorter than Germline.Input, because germline and
#'   clonotype sequences were trimmed to the same length before alignment.
#' * Common.Ancestor: common ancestor sequence, parsed from PHYLIP dnapars function output
#' * Trunk.Length: mean trunk length, representing the distance between the most recent
#'   common ancestor and germline sequence as a measure of the maturity of a lineage
#' * Tree: output tree in "phylo" format, loaded from by PHYLIP dnapars function output
#'
#' @examples
#'
#' data(bcrdata)
#' bcr_data <- bcrdata$data
#'
#' bcr_data %>%
#'   seqCluster(seqDist(bcr_data), .fixed_threshold = 3) %>%
#'   repGermline() %>%
#'   repAlignLineage(.min_lineage_sequences = 2, .align_threads = 2, .nofail = TRUE) %>%
#'   repClonalFamily(.threads = 2, .nofail = TRUE)
#' @export repClonalFamily
repClonalFamily <- function(.data, .threads = parallel::detectCores(), .nofail = FALSE) {
  if (!require_system_package("phylip", error_message = paste0(
    "repLineagePhylogeny requires PHYLIP app to be installed!\n",
    "Please install it as described here:\n",
    "https://evolution.genetics.washington.edu/phylip/install.html"
  ), .nofail, is.na(.data))) {
    return(NA)
  }

  results <- .data %>% apply_to_sample_or_list(
    process_dataframe,
    .with_names = TRUE,
    .validate = FALSE,
    .threads = .threads
  )
  return(results)
}

process_dataframe <- function(df, .threads, sample_name = NA) {
  for (column in c("Cluster", "Germline", "Alignment")) {
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

  if ("Aligned" %in% colnames(df)) {
    df %<>% filter(get("Aligned"))
  }

  if (nrow(df) == 0) {
    stop(
      "repClonalFamily: input dataframe ",
      optional_sample("in sample ", sample_name, " "),
      "doesn't contain any aligned clusters!"
    )
  }

  df <- df[c("Cluster", "Germline", "Alignment")]
  clusters_list <- split(df, seq(nrow(df)))

  results <- parallel::mclapply(
    clusters_list,
    process_cluster,
    mc.preschedule = FALSE,
    mc.cores = .threads
  ) %>%
    convert_nested_to_df()
  return(results)
}

process_cluster <- function(cluster_row) {
  cluster_name <- cluster_row[["Cluster"]]
  cluster_germline <- cluster_row[["Germline"]]
  # alignment should be extracted from 1-element list because of Alignment column format
  alignment <- cluster_row[["Alignment"]][[1]]

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
  outfile_lines <- scan(outfile_path, what = "", sep = "\n", quiet = TRUE)
  common_ancestor <- outfile_lines[grepl("         1   ", outfile_lines)]
  germline <- outfile_lines[grepl("1   germline", outfile_lines)]

  common_ancestor %<>%
    stringr::str_match(" +1 +(.+) *") %>%
    join_to_string()
  germline %<>%
    stringr::str_match(" +1 +germline\\s+[a-z]* *(.+) *") %>%
    join_to_string()
  trunk_length <- nchar(germline) - stringr::str_count(germline, stringr::fixed("."))

  unlink(temp_dir, recursive = TRUE)

  # return row of output dataframe as named list
  return(list(
    Cluster = cluster_name,
    Germline.Input = cluster_germline,
    Germline.Output = germline,
    Common.Ancestor = common_ancestor,
    Trunk.Length = trunk_length,
    Tree = tree
  ))
}

join_to_string <- function(matching_results) {
  matching_results[, 2] %>%
    paste(collapse = "") %>%
    stringi::stri_replace_all_fixed(" ", "")
}

convert_nested_to_df <- function(nested_results_list) {
  tree <- nested_results_list %>%
    lapply(magrittr::extract2, "Tree") %>%
    tibble(Tree = .)
  df <- nested_results_list %>%
    lapply(rlist::list.remove, "Tree") %>%
    purrr::map_dfr(~.) %>%
    cbind(tree)
  return(df)
}
