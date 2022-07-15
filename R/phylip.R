#' This function uses the PHYLIP package to make phylogenetic analysis. For making trees it uses
#' maximum parsimony methods.
#'
#' @concept phylip
#'
#' @aliases repClonalFamily
#'
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom purrr map_dfr
#' @importFrom rlist list.remove
#' @importFrom stringr str_match str_count fixed str_extract_all
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
#' @param .nofail Returns NA instead of stopping if PHYLIP is not installed.
#' Used to avoid raising errors in examples on computers where PHYLIP is not installed.
#'
#' @return
#'
#' Dataframe or list of dataframes (if input is a list with multiple samples).
#' The dataframe has these columns:
#' * Cluster: cluster name
#' * Germline.Input: germline sequence, like it was in the input; not trimmed and not aligned
#' * V.germline.nt: input germline V gene sequence
#' * J.germline.nt: input germline J gene sequence
#' * CDR3.germline.length: length of CDR3 in input germline
#' * Germline.Output: germline sequence, parsed from PHYLIP dnapars function output;
#'   it contains difference of germline from the common ancestor; "." characters mean
#'   matching letters. It's usually shorter than Germline.Input, because germline and
#'   clonotype sequences were trimmed to the same length before alignment.
#' * Common.Ancestor: common ancestor sequence, parsed from PHYLIP dnapars function output
#' * Trunk.Length: mean trunk length, representing the distance between the most recent
#'   common ancestor and germline sequence as a measure of the maturity of a lineage
#' * Tree: output tree in "phylo" format, loaded from by PHYLIP dnapars function output
#' * TreeStats: nested dataframe containing data about tree nodes, needed for visualization
#' * Sequences: nested dataframe containing all sequences for this combination of cluster
#'   and germline; it contains regions from original sequences, saved for
#'   repSomaticHypermutation() calculation, and also data needed for visualizations
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
#'   repClonalFamily(.threads = 2, .nofail = TRUE)
#' @export repClonalFamily
repClonalFamily <- function(.data, .threads = parallel::detectCores(), .nofail = FALSE) {
  if (!require_system_package("phylip", error_message = paste0(
    "repLineagePhylogeny requires PHYLIP app to be installed!\n",
    "Please install it as described here:\n",
    "https://evolution.genetics.washington.edu/phylip/install.html"
  ), .nofail, identical(.data, NA))) {
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
  required_columns <- c(
    "Cluster", "Germline", "V.germline.nt", "J.germline.nt", "CDR3.germline.length",
    "Alignment", "Sequences"
  )
  for (column in required_columns) {
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

  df <- df[required_columns]
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
  # alignment and sequences should be extracted from 1-element lists because of these columns format
  alignment <- cluster_row[["Alignment"]][[1]]
  sequences <- cluster_row[["Sequences"]][[1]]
  cdr3_germline_length <- cluster_row[["CDR3.germline.length"]]

  temp_dir <- file.path(tempdir(check = TRUE), uuid::UUIDgenerate(use.time = FALSE))
  dir.create(temp_dir)
  # workaround for phylip: it shows "Unexpected end-of-file" for too short sequence labels;
  # these \t are also used to read outfile as table
  rownames(alignment) %<>% paste0("\t")
  phangorn::write.phyDat(alignment, file.path(temp_dir, "infile"))
  system(
    paste0("sh -c \"cd ", temp_dir, "; phylip dnapars infile\""),
    input = (c("V", 1, 5, "Y"))
  ) %>%
    quiet()

  tree <- ape::read.tree(file.path(temp_dir, "outtree"))

  outfile_path <- file.path(temp_dir, "outfile")
  outfile_table <- read.table(outfile_path,
    sep = "\t", header = FALSE, na.strings = "", stringsAsFactors = FALSE,
    fill = TRUE, blank.lines.skip = TRUE
  )
  outfile_rows <- split(outfile_table, seq(nrow(outfile_table)))
  rm(outfile_table)
  header_line_1_idx <- which(grepl("From    To", outfile_rows, fixed = TRUE))[1]
  # remove start of file, keep only lines with sequences
  outfile_rows <- outfile_rows[-seq(1, header_line_1_idx + 1)]

  # iterate over rows, assemble sequences, fill table values
  tree_stats <- data.frame(matrix(
    ncol = 5, nrow = 0,
    dimnames = list(NULL, c("Name", "Type", "Clones", "Ancestor", "Sequence"))
  ))
  for (row in outfile_rows) {
    clones <- 1 # for all sequences except clonotypes
    col1_strings <- str_extract_all(row[[1]], "[[^\\s]]+")[[1]]
    if (is.na(row[[2]])) {
      ancestor <- NA
      seq_name <- col1_strings[1]
      if (seq_name == "Germline") {
        seq_type <- "Germline"
      } else {
        seq_type <- "CommonAncestor"
      }
      seq <- paste(col1_strings[-1], collapse = "")
    } else {
      col2_strings <- str_extract_all(row[[2]], "[[^\\s]]+")[[1]]
      ancestor <- col1_strings[1]
      seq_name <- col1_strings[2]
      seq <- paste(col2_strings[-1], collapse = "")
      if (seq_name == "Germline") {
        seq_type <- "Germline"
      } else if (startsWith(seq_name, "ID_")) {
        seq_type <- "Clonotype"
        # find clones value by sequence ID
        clones <- sequences[
          which(sequences["Clone.ID"] == as.integer(substring(seq_name, 4))),
        ][["Clones"]]
      } else if (ancestor == "Germline") {
        seq_type <- "CommonAncestor"
      } else {
        seq_type <- "Presumable"
      }
    }

    # add sequence to table if it's new, otherwise append nucleotides to the end
    existing_row <- tree_stats[which(seq_stats["Name"] == seq_name)]
    if (nrow(existing_row == 0)) {
      tree_stats[nrow(seq_stats) + 1, ] <- c(seq_name, seq_type, clones, ancestor, seq)
    } else {
      existing_row[["Sequence"]] %<>% paste0(seq)
    }
  }

  common_ancestor <- tree_stats[which(seq_stats["Type"] == "CommonAncestor")][1][["Sequence"]]
  germline <- tree_stats[which(seq_stats["Type"] == "Germline")][1][["Sequence"]]

  trunk_length <- nchar(germline) - cdr3_germline_length - str_count(germline, fixed("."))

  unlink(temp_dir, recursive = TRUE)

  # return row of output dataframe as named list
  return(list(
    Cluster = cluster_row[["Cluster"]],
    Germline.Input = cluster_row[["Germline"]],
    V.germline.nt = cluster_row[["V.germline.nt"]],
    J.germline.nt = cluster_row[["J.germline.nt"]],
    CDR3.germline.length = cdr3_germline_length,
    Germline.Output = germline,
    Common.Ancestor = common_ancestor,
    Trunk.Length = trunk_length,
    Tree = tree,
    TreeStats = tree_stats,
    Sequences = sequences,
  ))
}

convert_nested_to_df <- function(nested_results_list) {
  tree <- nested_results_list %>%
    lapply(magrittr::extract2, "Tree") %>%
    tibble(Tree = .)
  tree_stats <- nested_results_list %>%
    lapply(magrittr::extract2, "TreeStats") %>%
    tibble(TreeStats = .)
  sequences <- nested_results_list %>%
    lapply(magrittr::extract2, "Sequences") %>%
    tibble(Sequences = .)
  df <- nested_results_list %>%
    lapply(rlist::list.remove, c("Tree", "TreeStats", "Sequences")) %>%
    purrr::map_dfr(~.) %>%
    cbind(tree, tree_stats, sequences)
  df[["CDR3.germline.length"]] %<>% as.integer()
  df[["Trunk.Length"]] %<>% as.integer()
  return(df)
}
