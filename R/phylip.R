#' Builds a phylogenetic tree using the sequences of a clonal lineage
#'
#' @concept phylip
#'
#' @aliases repClonalFamily
#'
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom purrr map_dfr
#' @importFrom rlist list.remove
#' @importFrom stringr str_match str_count fixed str_extract_all str_length str_sub
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom utils capture.output
#' @importFrom parallel mclapply detectCores
#' @importFrom phangorn write.phyDat
#' @importFrom ape read.tree
#' @importFrom uuid UUIDgenerate
#' @importFrom data.table fread

#' @description This function uses the PHYLIP package to make phylogenetic analysis.
#' For making trees it uses maximum parsimony methods.
#'
#' @usage
#'
#' repClonalFamily(.data, .threads, .nofail)
#'
#' @param .data The data to be processed, output of repAlignLineage() function.
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
#' * Germline.Input: germline sequence, like it was in the input; not aligned
#' * Germline.Output: germline sequence, parsed from PHYLIP dnapars function output;
#'   it contains difference of germline from the common ancestor; "." characters mean
#'   matching letters
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
#'   repGermline(.threads = 1) %>%
#'   repAlignLineage(.min_lineage_sequences = 2, .align_threads = 2, .nofail = TRUE) %>%
#'   repClonalFamily(.threads = 1, .nofail = TRUE)
#' @export repClonalFamily
repClonalFamily <- function(.data, .threads = parallel::detectCores(), .nofail = FALSE) {
  if (!require_system_package(c("phylip", "dnapars"), error_message = paste0(
    "repLineagePhylogeny requires PHYLIP app to be installed!\n",
    "Please install it as described here:\n",
    "https://evolution.genetics.washington.edu/phylip/install.html"
  ), .nofail, has_class(.data, "step_failure_ignored"))) {
    return(get_empty_object_with_class("step_failure_ignored"))
  }

  results <- .data %>%
    apply_to_sample_or_list(
      process_dataframe,
      .with_names = TRUE,
      .validate = FALSE,
      .threads = .threads
    ) %>%
    add_class("clonal_family")
  return(results)
}

process_dataframe <- function(df, .threads, sample_name = NA) {
  required_columns <- c("Cluster", "Germline", "Alignment", "Sequences")
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

  if (nrow(df) == 0) {
    stop(
      "repClonalFamily: input dataframe ",
      optional_sample("in sample ", sample_name, " "),
      "doesn't contain any aligned clusters!"
    )
  }

  df <- df[required_columns]
  clusters_list <- split(df, seq(nrow(df)))

  results <- par_or_normal_lapply(
    clusters_list,
    process_cluster,
    mc.preschedule = FALSE,
    mc.cores = .threads
  ) %>%
    convert_nested_to_df()
  return(results)
}

process_cluster <- function(cluster_row) {
  # alignment, sequences and aa_frame_starts should be extracted from 1-element lists
  # because of these columns format
  alignment <- cluster_row[["Alignment"]][[1]]
  sequences <- cluster_row[["Sequences"]][[1]]
  cluster_name <- cluster_row[["Cluster"]]

  fsep <- if (.Platform$OS.type == "windows") "\\" else "/"
  shell <- if (.Platform$OS.type == "windows") "powershell /c " else "sh -c "
  temp_dir <- file.path(tempdir(check = TRUE), uuid::UUIDgenerate(use.time = FALSE), fsep = fsep)
  dir.create(temp_dir)
  # workaround for phylip: it shows "Unexpected end-of-file" for too short sequence labels;
  # these \t are also used to read outfile as table
  rownames(alignment) %<>% paste0("\t")
  phangorn::write.phyDat(alignment, file.path(temp_dir, "infile", fsep = fsep))
  dnapars <- if (Sys.which("phylip") == "") "dnapars" else "phylip dnapars"
  system(
    paste0(shell, "\"cd ", temp_dir, "; ", dnapars, " infile\""),
    input = (c("V", 1, 5, ".", "Y"))
  ) %>%
    quiet(capture_output = TRUE)

  tree <- ape::read.tree(file.path(temp_dir, "outtree", fsep = fsep))

  outfile_path <- file.path(temp_dir, "outfile", fsep = fsep)
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
    ncol = 7, nrow = 0,
    dimnames = list(NULL, c(
      "Name", "Type", "Clones", "Ancestor", "DistanceNT", "DistanceAA", "Sequence"
    ))
  )) %>%
    add_class("clonal_family_tree")
  for (row in outfile_rows) {
    clones <- 1 # for all sequences except clonotypes
    col1_strings <- str_extract_all(row[[1]], "[[^\\s]]+")[[1]]
    if (is.na(row[[2]])) {
      # if second element is not a number, then ancestor is NA
      if (is.na(quiet(as.numeric(col1_strings[2])))) {
        ancestor <- NA
        seq_name <- col1_strings[1]
        seq <- paste(col1_strings[-1], collapse = "")
      } else {
        ancestor <- col1_strings[1]
        seq_name <- col1_strings[2]
        col1_strings <- col1_strings[-c(1, 2)]
        if (col1_strings[1] %in% c("yes", "no", "maybe")) {
          col1_strings <- col1_strings[-1]
        }
        seq <- paste(col1_strings, collapse = "")
      }
    } else {
      col2_strings <- str_extract_all(row[[2]], "[[^\\s]]+")[[1]]
      ancestor <- col1_strings[1]
      seq_name <- col1_strings[2]
      seq <- paste(col2_strings[-1], collapse = "")
    }

    if (seq_name == "Germline") {
      seq_type <- "Germline"
    } else if (startsWith(seq_name, "ID_")) {
      seq_type <- "Clonotype"
      # find clones value by sequence ID
      clones <- sequences[
        which(sequences["Clone.ID"] == as.integer(substring(seq_name, 4))),
        "Clones"
      ]
    } else if (identical(ancestor, "Germline")) {
      seq_type <- "CommonAncestor"
    } else {
      seq_type <- "Presumable"
    }

    # add sequence to table if it's new, otherwise append nucleotides to the end
    if (nrow(tree_stats[which(tree_stats["Name"] == seq_name), ]) == 0) {
      tree_stats[nrow(tree_stats) + 1, ] <- c(seq_name, seq_type, clones, ancestor, 0, 0, seq)
    } else {
      tree_stats[which(tree_stats["Name"] == seq_name), "Sequence"] %<>% paste0(seq)
    }
  }

  # force germline to be root if phylip had set "1" as its ancestor
  germline_ancestor <- tree_stats[which(tree_stats["Type"] == "Germline"), ][1, "Ancestor"]
  if (!is.na(germline_ancestor)) {
    ancestor_of_ancestor <-
      tree_stats[which(tree_stats["Name"] == germline_ancestor), ][1, "Ancestor"]
    if (is.na(ancestor_of_ancestor)) {
      tree_stats[which(tree_stats["Name"] == "Germline"), ][1, "Ancestor"] <- NA
      tree_stats[which(tree_stats["Name"] == germline_ancestor), ][1, "Type"] <- "CommonAncestor"
      tree_stats[which(tree_stats["Name"] == germline_ancestor), ][1, "Ancestor"] <- "Germline"
    } else {
      warning(
        "Germline's ancestor is ",
        germline_ancestor,
        ", and its' ancestor is ",
        ancestor_of_ancestor,
        ": tree of cluster ",
        cluster_name,
        " can't be corrected!"
      )
    }
  }

  # remove ancestors that are not in names
  names <- tree_stats[["Name"]]
  ancestors <- tree_stats[["Ancestor"]]
  for (missing_name in ancestors[!ancestors %in% names]) {
    tree_stats["Ancestor"][tree_stats["Ancestor"] == missing_name] <- NA
  }
  # rename CommonAncestor and Presumable nodes, both in Name and Ancestor columns
  renamed_rows <- tree_stats[which(tree_stats[["Type"]] %in% c("CommonAncestor", "Presumable")), ]
  old_names <- renamed_rows[["Name"]]
  new_names <- as.list(paste0(renamed_rows[["Type"]], "_", renamed_rows[["Name"]]))
  names(new_names) <- old_names
  for (column in c("Name", "Ancestor")) {
    for (row in which(tree_stats[[column]] %in% old_names)) {
      tree_stats[[row, column]] <- new_names[[tree_stats[[row, column]]]]
    }
  }

  germline <- tree_stats[which(tree_stats["Type"] == "Germline"), ][1, "Sequence"]
  common_ancestor <- tree_stats[which(tree_stats["Type"] == "CommonAncestor"), ][1, "Sequence"]

  # calculate distances of all sequences from their ancestors
  v_nt_length <- str_length(paste0(
    sequences[[1, "FR1.nt"]], sequences[[1, "CDR1.nt"]], sequences[[1, "FR2.nt"]],
    sequences[[1, "CDR2.nt"]], sequences[[1, "FR3.nt"]]
  ))
  j_nt_length <- str_length(sequences[[1, "FR4.nt"]])
  v_aa_length <- str_length(sequences[[1, "V.aa"]])
  j_aa_length <- str_length(sequences[[1, "J.aa"]])

  for (row in seq_len(nrow(tree_stats))) {
    if (!has_no_data(tree_stats[row, "Ancestor"])) {
      seq <- tree_stats[row, "Sequence"]
      ancestor_name <- tree_stats[row, "Ancestor"]
      ancestor <- tree_stats[which(tree_stats["Name"] == ancestor_name), ][1, "Sequence"]
      seq_v <- str_sub(seq, 1, v_nt_length)
      seq_j <- str_sub(seq, -j_nt_length)
      ancestor_v <- str_sub(ancestor, 1, v_nt_length)
      ancestor_j <- str_sub(ancestor, -j_nt_length)
      seq_nt_chars <- strsplit(paste0(seq_v, seq_j), "")[[1]]
      ancestor_nt_chars <- strsplit(paste0(ancestor_v, ancestor_j), "")[[1]]
      if (length(seq_nt_chars) != length(ancestor_nt_chars)) {
        warning(
          "Sequence and ancestor lengths are different; ",
          "something is wrong in repClonalFamily calculation!\n",
          "seq_nt_chars = ", seq_nt_chars, "\nancestor_nt_chars = ", ancestor_nt_chars
        )
      }
      tree_stats[row, "DistanceNT"] <- sum(seq_nt_chars != ancestor_nt_chars)

      seq_aa <- bunch_translate(seq, .two.way = FALSE, .ignore.n = TRUE)
      ancestor_aa <- bunch_translate(ancestor, .two.way = FALSE, .ignore.n = TRUE)
      seq_v_aa <- str_sub(seq_aa, 1, v_aa_length)
      seq_j_aa <- str_sub(seq_aa, -j_aa_length)
      ancestor_v_aa <- str_sub(ancestor_aa, 1, v_aa_length)
      ancestor_j_aa <- str_sub(ancestor_aa, -j_aa_length)
      seq_aa_chars <- strsplit(paste0(seq_v_aa, seq_j_aa), "")[[1]]
      ancestor_aa_chars <- strsplit(paste0(ancestor_v_aa, ancestor_j_aa), "")[[1]]
      if (length(seq_aa_chars) != length(ancestor_aa_chars)) {
        warning(
          "Sequence and ancestor lengths are different; ",
          "something is wrong in repClonalFamily calculation!\n",
          "seq_aa_chars = ", seq_aa_chars, "\nancestor_aa_chars = ", ancestor_aa_chars
        )
      }
      tree_stats[row, "DistanceAA"] <- sum(seq_aa_chars != ancestor_aa_chars)
    }
  }

  for (column in c("Clones", "DistanceNT", "DistanceAA")) {
    tree_stats[[column]] %<>% as.integer()
  }

  trunk_length <- tree_stats[which(tree_stats["Ancestor"] == "Germline"), ][1, "DistanceNT"]

  unlink(temp_dir, recursive = TRUE)

  # return row of output dataframe as named list
  return(list(
    Cluster = cluster_name,
    Germline.Input = cluster_row[["Germline"]],
    Germline.Output = germline,
    Common.Ancestor = common_ancestor,
    Trunk.Length = trunk_length,
    Tree = tree,
    TreeStats = tree_stats,
    Sequences = sequences
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
  # fix column types after dataframe rebuilding
  df[["Trunk.Length"]] %<>% as.integer()
  df %<>% add_class("clonal_family_df")
  return(df)
}
