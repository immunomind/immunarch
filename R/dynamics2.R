#' Track clonotypes across time and data points
#'
#' @concept dynamics
#'
#' @importFrom dplyr arrange desc distinct rowwise
#' @importFrom tibble tibble
#' @importFrom factoextra hcut
#' @importFrom tibble rownames
#'
#' @aliases trackClonotypes2
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
#' @param .num An integer that specifies how many of the top abundant clonotypes from each sample to include. Default is 1000. 
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
#' @export trackClonotypes2
#' 
trackClonotypes2 <- function(.data, .which = list(), .col = "aa", .norm = TRUE, .k = tibble) {

  if (!has_class(.data, "list")) {
    stop("Error: please pass a list with immune repertoires to track clonotypes.")
  }
  if (length(.data) < 2) {
    stop("Error: please pass a list with 2 or more immune repertoires to track clonotypes.")
  }

  col_param <- .col

  .col <- unlist(strsplit(.col, split = "\\+"))
  .col <- sapply(.col, switch_type, USE.NAMES = FALSE)

  track <- FALSE

  if (has_class(.which, "list")) {

    if (length(.which) != 2) {
      stop("Error: please pass a list with two elements for the .which argument. Run ?trackClonotypes for more details in the documentation.")
    }

    ########## BEGINNING OF TRACK CLONOTYPES ###########

    if (has_class(.which[[1]], "numeric") & has_class(.which[[2]], "tbl")) {
      print("IN IF LOOP")
      track <- TRUE
      n_clonotypes <- .which[[1]]
      meta <- .which[[2]]

      if (!("Timepoints" %in% colnames(meta))) {
        stop("Error: please pass a meta data with Timepoints column to .which argument. Run ?trackClonotypes for more details in the documentation.")

      }

      timepoints <- meta$Sample

      if (n_clonotypes < 1) {
        stop("Error: please enter a non-zero positive number of clonotypes to the .which argument.")
      }

      #### PRE-PROCESSING
      # subset data by the timepoints
      subset_df <- .data[timepoints]

      # for each sample in the list in .timepoint, get top .num clonotypes
      subset_df <- top(subset_df, .n = .num)

      subset_df <- lapply(subset_df, head, .num)

      top_df <- pubRep(subset_df, col_param, .verbose = FALSE)

      row.names(top_df) <- as.vector(unlist(top_df[, 1]))

      top_df <- top_df[, c(.col, "Samples") := NULL]

      top_df <- top_df / colSums(top_df, na.rm = TRUE)

      #### NORMALISATION
      # adding column with max proportion for each clonotype across time points
      m <- top_df %>% rowwise() %>% mutate(max_prop = max(across(timepoints), na.rm = TRUE))

      # dividing proportions by the max proportion across all time points
      m <- m %>%
      ungroup() %>%
      mutate(across(timepoints, ~ . / max_prop))

      # removing max prop column
      m$max_prop <- NULL

      m[is.na(m)] <- 0

      m.matrix <- data.matrix(m)

      result_df <- NULL

      result_df$pca <- immunr_pca(m.matrix)

      ### GENERATE TRAJECTORIES

      res <- hcut(as.matrix(m), k = 3, stand = FALSE)

      m$clust <- res$cluster

      result_df$cluster <- res

      result_df$data <- m

      m <- m %>% reshape2::melt(id.vars = "clust")

      names(m)[2] <- "Sample"

    #   print("Before merge")

      m <- merge(m, meta, by=c("Sample"))

      print(head(m))

      result_df$melt <- m

      result_df$traj <- m %>% group_by(clust, Timepoints) %>% summarise(mean_traj = mean(value), n_traj = n(), SE = sd(value) / sqrt(n()), SE_scaled = 2.96 * sd(value) / sqrt(n()))

      add_class(result_df, "immunr_trajectories")

    } else {
      # Option 1
      print("HERE OPTION 1")
    #   target_df <- .which[[1]]
    #   n_clonotypes <- .which[[2]]

    #   # ToDo: replace Clones with IMMCOL$count
    #   target_seq <- .data[[target_df]] %>%
    #   arrange(desc(Clones)) %>%
    #   select(.col) %>%
    #   head(n_clonotypes) %>%
    #   distinct() %>%
    #   collect(n = Inf)
    }

  } else if (has_class(.which, "character")) {
            print("HERE OPTION 2")

    # Option 2
    # target_seq <- tibble(.which) %>% distinct()

    # ToDo: extract sequence column names
    # names(target_seq)[1] <- .col[1]
  } else {
    # Option 3
          print("HERE OPTION 3")

    # target_seq <- .which %>% distinct()
    # .col <- colnames(target_seq)

    # for (col_name in .col) {
    #   if (!(col_name %in% names(.data[[1]]))) {
    #     stop("Error: can't find a column with the name \"", col_name, "\" in the input data. Please check column names for the .which argument.")
    #   }
    }

     if (track == FALSE) {
      print("HERE")
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
     

  } else {
  add_class(result_df, "immunr_trajectories")

  }

}

 


