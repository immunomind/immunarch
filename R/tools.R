#' Check and normalise distributions
#'
#' @concept utility_private
#'
#' @description
#' Check if the given .data is a distribution and normalise it if necessary with an optional laplace correction.
#'
#' @param .data Numeric vector of values.
#' @param .do.norm One of the three values - NA, TRUE or FALSE. If NA then check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. if TRUE then do normalisation and laplace
#' correction. If FALSE then don't do normalisaton and laplace correction.
#' @param .laplace Value for the laplace correction.
#' @param .na.val Replace all NAs with this value.
#' @param .warn.zero if TRUE then the function checks if in the resulted vector (after normalisation)
#' are any zeros, and prints a warning message if there are some.
#' @param .warn.sum if TRUE then the function checks if the sum of resulted vector (after normalisation)
#' is equal to one, and prints a warning message if not.
#'
#' @return Numeric vector.
#'
#' @section Developer Examples:
#' immunarch:::check_distribution(c(1, 2, 3))
#' immunarch:::check_distribution(c(1, 2, 3), TRUE)
#' immunarch:::check_distribution(c(1, 2, 3), FALSE)
check_distribution <- function(.data, .do.norm = NA, .laplace = 1, .na.val = 0, .warn.zero = FALSE, .warn.sum = TRUE) {
  if (sum(is.na(.data)) == length(.data)) {
    warning("Error! Input vector is completely filled with NAs. Check your input data to avoid this. Returning a vector with zeros.\n")
    return(rep.int(0, length(.data)))
  }

  if (is.na(.do.norm)) {
    .data[is.na(.data)] <- .na.val
    if (sum(.data) != 1) {
      .data <- .data + .laplace
      .data <- prop.table(.data + .laplace)
    }
  } else if (.do.norm) {
    .data[is.na(.data)] <- .na.val
    .data <- prop.table(.data + .laplace)
    .warn.sum <- FALSE
  }

  if (.warn.zero && (0 %in% .data)) {
    warningText <- paste("Warning! There are", sum(which(.data == 0)), "zeros in the input vector. Function may produce incorrect results.\n")
    if (.laplace != 1) {
      warningText <- paste(warningText, "To fix this try to set .laplace = 1 or any other small number in the function's parameters\n")
    } else {
      warningText <- paste(warningText, "To fix this try to set .laplace to any other small number in the function's parameters\n")
    }
    warning(warningText)
  }

  if (.warn.sum && sum(.data) != 1) {
    if (abs(sum(.data) - 1) < 1e-14) {
      # message("Note: difference between the sum of the input vector and 1 is ", (sum(.data) - 1), ", which may be caused by internal R subroutines and may not affect the result at all.\n")
      # Just skip this case - it means all is OK
    } else {
      warningText <- "Warning! Sum of the input vector is NOT equal to 1. Function may produce incorrect results.\n"
      if (!isTRUE(.do.norm)) warningText <- paste(warningText, "To fix this try to set .do.norm = TRUE in the function's parameters.\n")
      warning(warningText)
    }
  }

  .data
}


#' Add a new class attribute
#'
#' @concept utility_private
#'
#' @param .obj R object.
#' @param .class String with the desired class name.
#'
#' @return
#' Input object with additional class \code{.class}.
#'
#' @section Developer Examples:
#' tmp <- "abc"
#' class(tmp)
#' tmp <- immunarch:::add_class(tmp, "new_class")
#' class(tmp)
add_class <- function(.obj, .class) {
  class(.obj) <- c(.class, class(.obj))
  .obj
}


#' Check for the specific class
#'
#' @concept utility_private
#'
#' @description
#' A function to check if an input object has a specific class name.
#'
#' @param .data Any R object.
#' @param .class Character vector. Specifies a class name to check against.
#'
#' @return
#' Logical value.
#'
#' @section Developer Examples:
#' tmp <- "abc"
#' immunarch:::has_class(tmp, "new_class")
#' tmp <- immunarch:::add_class(tmp, "new_class")
#' immunarch:::has_class(tmp, "new_class")
has_class <- function(.data, .class) {
  .class %in% class(.data)
}


#' Set and update progress bars
#'
#' @concept utility_private
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @aliases set_pb add_pb
#'
#' @usage
#' set_pb(.max)
#'
#' add_pb(.pb, .value = 1)
#'
#' @param .max Integer. Maximal value of the progress bar.
#' @param .pb Progress bar object from \code{set_pb}.
#' @param .value Numeric. Value to add to the progress bar at each step.
#'
#' @return
#' An updated progress bar.
#'
#' @section Developer Examples:
#' pb <- immunarch:::set_pb(100)
#' immunarch:::add_pb(pb, 25)
#' immunarch:::add_pb(pb, 25)
#' immunarch:::add_pb(pb, 25)
#' immunarch:::add_pb(pb, 25)
#' close(pb)
set_pb <- function(.max) {
  txtProgressBar(min = 0, max = .max, style = 3)
}

add_pb <- function(.pb, .value = 1) {
  setTxtProgressBar(.pb, .pb$getVal() + .value)
}


#' Get a column's name using the input alias
#'
#' @concept utility_private
#'
#' @param x Character vector of length 1.
#' @param .verbose If TRUE then print any issues or wanings.
#'
#' @return
#' A string with the column name.
#'
#' @section Developer Examples:
#' immunarch:::.quant_column_choice("count")
#' immunarch:::.quant_column_choice("freq")
.quant_column_choice <- function(x, .verbose = TRUE) {
  x <- switch(x[1],
    count = IMMCOL$count,
    Count = IMMCOL$count,
    prop = IMMCOL$prop,
    Prop = IMMCOL$prop,
    proportion = IMMCOL$prop,
    Proportion = IMMCOL$prop,
    freq = IMMCOL$prop,
    {
      message("You have specified an invalid column identifier. Active column: Clones")
      IMMCOL$count
    }
  )
  x
}


#' Copy the upper matrix triangle to the lower one
#'
#' @concept utility_private
#'
#' @param .mat Matrix.
#'
#' @return
#' Matrix with its upper tri part copied to the lower tri part.
#'
#' @section Developer Examples:
#' mat <- matrix(0, 3, 3)
#' mat
#' mat[1, 3] <- 1
#' mat <- immunarch:::matrixdiagcopy(mat)
#' mat
matrixdiagcopy <- function(.mat) {
  .mat[lower.tri(.mat)] <- t(.mat)[lower.tri(.mat)]
  .mat
}


#' Nucleotide to amino acid sequence translation
#'
#' @concept preprocessing
#'
#' @aliases bunch_translate translate_bunch
#'
#' @param .seq Vector or list of strings.
#' @param .two.way Logical. If TRUE (default) then translate from the both ends (like MIXCR).
#'
#' @return
#' Character vector of translated input sequences.
#'
#' @examples
#' data(immdata)
#' head(bunch_translate(immdata$data[[1]]$CDR3.nt))
#' @export
bunch_translate <- function(.seq, .two.way = TRUE) {
  .seq <- toupper(.seq)
  .seq[grepl("N", .seq)] <- NA

  sapply(.seq, function(y) {
    if (!is.na(y)) {
      ny <- nchar(y)
      ny3 <- ny %/% 3
      tmp <- ""
      if (.two.way) {
        if (ny %% 3 != 0) {
          tmp <- paste0(rep("N", times = 3), collapse = "")
        }
        y <- paste0(substr(y, 1, 3 * ((ny3 %/% 2) + (ny %% 2))),
          tmp,
          substr(y, 3 * ((ny3 %/% 2) + (ny3 %% 2)) + (ny %% 3) + 1, ny),
          collapse = ""
        )
      } else {
        y <- substring(y, seq(1, nchar(y) - 2, 3), seq(3, nchar(y), 3))
      }
      paste0(AA_TABLE[unlist(strsplit(gsub("(...)", "\\1_", y), "_"))], collapse = "")
    } else {
      NA
    }
  }, USE.NAMES = FALSE)
}

translate_bunch <- bunch_translate


check_group_names <- function(.meta, .by) {
  names_to_check <- c()
  if (is.null(.by)) {
    names_to_check <- .by
  } else {
    names_to_check <- names(.by)
  }

  for (i in 1:length(names_to_check)) {
    if (!(names_to_check[i] %in% colnames(.meta))) {
      message("Check failed: '", names_to_check[i], "' not in the metadata table!")
      return(FALSE)
    }
  }

  return(TRUE)
}


#' Get a character vector of samples' groups from the input metadata file
#'
#' @concept utility_private
#'
#' @importFrom stringr str_sort str_order
#'
#' @param .by Character vector. Specify a column or columns in the input metadata to group by.
#' @param .metadata Metadata object.
#' @param .sep Character vector. Defines a separator between groups if more than one group passed in \code{.by}.
#'
#' @return
#' Character vector with group names.
#'
#' @section Developer Examples:
#' immunarch:::group_from_metadata("Status", data.frame(Status = c("A", "A", "B", "B", "C")))
group_from_metadata <- function(.by, .metadata, .sep = "; ") {
  if (length(.by) == 1) {
    collect(select(.metadata, .by))[[1]]
  } else {
    do.call(paste, c(list(sep = .sep), lapply(1:length(.by), function(i) {
      collect(select(.metadata, .by[i]))[[1]]
    })))
  }
}


# .meta == NA => .by is a vector of values to group by
# .meta != NA => .by is a name of the column in the metadata file
# .meta == NA & .by == NA => just choose the default column .data.sample.col for grouping
process_metadata_arguments <- function(.data, .by, .meta = NA, .data.sample.col = "Sample",
                                       .meta.sample.col = "Sample") {
  .data[[.data.sample.col]] <- as.character(.data[[.data.sample.col]])
  if (!is.na(.by)[1]) {
    if (!is.na(.meta)[1]) {
      data_groups <- group_from_metadata(.by, .meta)
      # group_name = .by
      group_name <- paste0(.by, collapse = "; ")
      is_grouped <- TRUE
      data_group_names <- .meta[[.meta.sample.col]]
    } else {
      if (length(.by) == length(.data[[.data.sample.col]])) {
        data_groups <- as.character(.by)
        data_group_names <- .data[[.data.sample.col]]
        group_name <- "Group"
        is_grouped <- TRUE
      } else {
        stop("Error: length of the input vector '.by' isn't the same as the length of the input data. Please provide vector of the same length.")
      }
    }
  } else {
    data_groups <- unique(.data[[.data.sample.col]])
    group_name <- .data.sample.col
    is_grouped <- FALSE
    data_group_names <- unique(.data[[.data.sample.col]])

    if (length(data_groups) != length(data_group_names)) {
      stop("Error: number of samples doesn't equal to the number of samples in the metadata")
    }
  }
  names(data_groups) <- data_group_names
  group_vec <- data_groups[.data[[.data.sample.col]]]
  # group_column = stringr::str_sort(data_groups[.data[[.data.sample.col]]], numeric = TRUE)
  group_vec_sorted <- stringr::str_sort(group_vec, numeric = TRUE)
  group_column <- factor(group_vec, levels = unique(group_vec_sorted))

  list(groups = data_groups, group_column = group_column, group_names = data_group_names, name = group_name, is_grouped = is_grouped)
}


rename_column <- function(.data, .old, .new) {
  colnames(.data)[match(.old, colnames(.data))] <- .new
  .data
}


#' Apply function to each pair of data frames from a list.
#'
#' @concept utility_public
#'
#' @aliases apply_symm apply_asymm
#'
#' @description
#' Apply the given function to every pair in the given datalist. Function either
#' symmetrical (i.e. fun(x,y) == fun(y,x)) or assymmetrical (i.e. fun(x,y) != fun(y,x)).
#'
#' @usage
#' apply_symm(.datalist, .fun, ..., .diag = NA, .verbose = TRUE)
#'
#' apply_asymm(.datalist, .fun, ..., .diag = NA, .verbose = TRUE)
#'
#' @param .datalist List with some data.frames.
#' @param .fun Function to apply, which return basic class value.
#' @param ... Arguments passsed to .fun.
#' @param .diag Either NA for NA or something else != NULL for .fun(x,x).
#' @param .verbose if TRUE then output a progress bar.
#'
#' @return Matrix with values M[i,j] = fun(datalist[i], datalist[j])
#'
#' @examples
#' data(immdata)
#' apply_symm(immdata$data, function(x, y) {
#'   nrow(x) + nrow(y)
#' })
#' @export apply_symm apply_asymm
apply_symm <- function(.datalist, .fun, ..., .diag = NA, .verbose = TRUE) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set_pb(length(.datalist)^2 / 2 + length(.datalist) / 2)
  for (i in 1:length(.datalist)) {
    for (j in i:length(.datalist)) {
      if (i == j && is.na(.diag)) {
        res[i, j] <- NA
      } else {
        res[i, j] <- .fun(.datalist[[i]], .datalist[[j]], ...)
      }
      if (.verbose) add_pb(pb)
    }
  }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  matrixdiagcopy(res)
}

apply_asymm <- function(.datalist, .fun, ..., .diag = NA, .verbose = TRUE) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set_pb(length(.datalist)^2)
  for (i in 1:length(.datalist)) {
    for (j in 1:length(.datalist)) {
      if (i == j && is.na(.diag)) {
        res[i, j] <- NA
      } else {
        res[i, j] <- .fun(.datalist[[i]], .datalist[[j]], ...)
      }
      if (.verbose) add_pb(pb)
    }
  }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  res
}


#' Return a column's name
#'
#' @concept utility_private
#'
#' @aliases switch_type process_col_argument
#'
#' @usage
#' switch_type(type)
#'
#' process_col_argument(.col)
#'
#' @param .col A string that specifies the column(s) to be processed. Pass one of the
#' following strings, separated by the plus sign: "nt" for nucleotide sequences,
#' "aa" for amino acid sequences, "v" for V gene segments, "j" for J gene segments.
#' @param type Character. Specifies the column to choose:
#' "nt" chooses the CDR3 nucleotide column,
#' "aa" chooses the CDR3 amino acid column,
#' "v" chooses the V gene segment column,
#' "j" chooses the J gene segment column.
#'
#' @return
#' A column's name.
#'
#' @section Developer Examples:
#' immunarch:::switch_type("nuc")
#' immunarch:::switch_type("v")
switch_type <- function(type) {
  switch(tolower(type),
    nuc = IMMCOL$cdr3nt, # "nuc" left for the compatability with older versions
    nt = IMMCOL$cdr3nt,
    v = IMMCOL$v,
    j = IMMCOL$j,
    aa = IMMCOL$cdr3aa,
    stop("Error: unknown column identifier ", type, '. Please pass one of the following: "nt", "aa", "v" or "j".')
  )
}

process_col_argument <- function(.col) {
  .col <- unlist(strsplit(gsub("nuc", "nt", .col), split = "\\+"))
  sapply(names(sort(c(nt = 1, aa = 2, v = 3, j = 4)[.col])), switch_type, USE.NAMES = FALSE)
}


return_segments <- function(.gene) {
  stringr::str_replace_all(.gene, "\\*[[:digit:]]+", "")
}

return_families <- function(.gene) {
  stringr::str_replace_all(return_segments(.gene), "\\-[[:digit:]]+", "")
}

load_segments <- function(.path, .alias, .gene_df = GENE_SEGMENTS, .filter = NA) {
  segm <- readr::read_tsv(.path)
  setnames(segm, "#species", "species", skip_absent = TRUE)
  if (!is.na(.filter)) {
    segm <- segm %>% filter(species == .filter)
  }
  # print(table(segm$species))
  # print(segm %>% filter(gene == "trdv", species == "BosTaurus"))
  setnames(segm, "id", "allele_id", skip_absent = TRUE)
  segm$gene <- tolower(paste0(segm$gene, substr(segm$segment, 1, 1)))
  segm$family_id <- return_families(segm$allele_id)
  segm$segment_id <- return_segments(segm$allele_id)
  segm$alias <- .alias
  segm <- segm[c("alias", "species", "gene", "family_id", "segment_id", "allele_id", "reference_point", "sequence")]
  segm <- segm[order(segm$gene, segm$allele_id), ]
  new_gs <- rbind(segm, .gene_df[.gene_df$alias != .alias, ])
}

as_numeric_or_fail <- function(.string) {
  result <- as.numeric(.string)
  if (is.na(result)) {
    stop(paste0("\"", .string, "\" is not a valid number."))
  }
  return(result)
}

all_bcr_columns_present <- function(.data) {
  expected_columns <- c(
    IMMCOL$order,
    IMMCOL_EXT$bestv, IMMCOL_EXT$bestj, IMMCOL_EXT$cdr3s, IMMCOL_EXT$cdr3e, IMMCOL_EXT$c
  )
  return(all(expected_columns %in% colnames(.data)))
}

# get genes from original column, remove alleles and write to target column
add_column_without_alleles <- function(.data, .original_colname, .target_colname) {
  if (.original_colname %in% colnames(.data)) {
    .data[[.target_colname]] <- gsub(
      ", ", ",",
      .data[[.original_colname]]
    )
    .data[[.target_colname]] <- gsub(
      "([*][[:digit:]]*)([(][[:digit:]]*[.,]*[[:digit:]]*[)])", "",
      .data[[.target_colname]]
    )
    .data[[.target_colname]] <- gsub(
      ",", ", ",
      .data[[.target_colname]]
    )
    .data[[.target_colname]] <- sapply(
      strsplit(.data[[.target_colname]], ", ", useBytes = TRUE),
      # No sorting because MiXCR outputs segments in a specific order
      function(x) paste0(unique(x), collapse = ", ")
    )
    .data[[.target_colname]] <- gsub(
      "[*][[:digit:]]*", "",
      .data[[.target_colname]]
    )
  } else {
    stop(paste0(
      "Trying to remove alleles from missing column ", .original_colname,
      ", available columns: ", colnames(.data)
    ))
  }
  .data
}
