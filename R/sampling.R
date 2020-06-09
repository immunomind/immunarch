#######
# WIP #
#######

# error_correction <- function () {
#
# }

# decontamination <- function () {
#
# }


#' Downsampling and resampling of immune repertoires
#'
#' @concept preprocessing
#'
#' @importFrom Rcpp cppFunction
#' @importFrom stats rmultinom
#' @importFrom dplyr tally
#'
#' @aliases repSample
#'
#' @description
#' Sample (downsample) repertoires using different approches.
#'
#' @param .data The data to be processed. Can be \link{data.frame},
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
#' @param .method Character. Name of a sampling method. See "Description" for more details. Default value is "downsample"
#' that downsamples repertoires to the number of clones (i.e., reads / UMIs) that the smallest repertoire has, if user
#' doesn't pass any value to the ".n" argument.
#'
#' @param .n Integer. Number of clones / clonotypes / reads / UMIs to choose, depending on the method. Pass NA to sample
#' repertoires to the size of the smallest repertoire in the ".data".
#'
#' @param .prob Logical. If TRUE then sample clonotypes with probability weights equal to their number of clones. Used only if
#' ".method" is "sample".
#'
#' @return Subsampled immune repertoires or a list of subsampled immune repertoires.
#'
#' @details
#' If \code{.method} is "downsample" then \code{repSample} chooses \code{.n} clones (not clonotypes!) from the input repertoires without any probabilistic simulation,
#' but exactly computing each choosed clones. Such approach is is more consistent and biologically pleasant than
#' an output from the function if \code{.method} is "resample".
#'
#' If \code{.method} is "resample" then \code{repSample} uses multinomial distribution to compute the number of occurences for each cloneset.
#' then it removes zero-number clonotypes and return the resulting data frame. Probabilities for \code{rmultinom} for each cloneset
#' is a percentage of this cloneset in the "Proportion" column. It's a some sort of simulation of how clonotypes are chosen from the organisms.
#'
#' if \code{.method} is "sample" then \code{repSample} chooses \code{.n} clonotypes (not clones!) randomly. Depending on the
#' \code{.prob} argument, the function chooses clonotypes either according to their size (if \code{.prob} is TRUE, by default),
#' or each clonotype has an equal chance to be choosed (if \code{.prob} is FALSE). Note that sampling is done without replacing.
#'
#' @seealso \link{rmultinom}, \link{clonal_proportion}
#'
#' @examples
#' data(immdata)
#' # Downsampling to 1000 clones (not clonotypes!)
#' tmp <- repSample(immdata$data[[1]], .n = 1000)
#' sum(tmp$Clones)
#'
#' # Downsampling to 1000 clonotypes
#' tmp <- repSample(immdata$data[[1]], "sample", .n = 1000)
#' nrow(tmp)
#'
#' # Downsampling to the smallest repertoire by clones (not clonotypes!)
#' tmp <- repSample(immdata$data[c(1, 2)])
#' sum(tmp[[1]]$Clones)
#' sum(tmp[[2]]$Clones)
#'
#' # Downsampling to the smallest repertoire by clonotypes
#' tmp <- repSample(immdata$data[c(1, 2)], "sample")
#' nrow(tmp[[1]]$Clones)
#' nrow(tmp[[2]]$Clones)
#' @export repSample
repSample <- function(.data, .method = c("downsample", "resample", "sample"), .n = NA, .prob = TRUE) {
  .method <- .method[1]

  if (is.na(.n)) {
    if ((has_class(.data, "list") && (length(.data) == 1)) || (!has_class(.data, "list"))) {
      warning("Warning: no value passed for argument .n when input data is a single data frame. Setting .n to 1000.")
      .n <- 1000
    }
  }

  if (!(.method %in% c("downsample", "resample", "sample"))) {
    stop('Error: wrong resampling method. Please provide one of the following: "downsample", "resample" or "sample".')
  }

  if (.method == "downsample") {
    # choose the repertoire with the smallest amount of clones
    if (is.na(.n)) {
      clones_vec <- sapply(.data, function(df) sum(df[[IMMCOL$count]]))
      i <- which.min(clones_vec)
      .n <- sum(.data[[i]][[IMMCOL$count]])
    }

    downsample(.data, .n = .n)
  } else if (.method == "resample") {
    # choose the repertoire with the smallest amount of clones
    if (is.na(.n)) {
      i <- which.min(sapply(.data, function(df) sum(df[[IMMCOL$count]])))
      .n <- sum(.data[[i]][[IMMCOL$count]])
    }

    resample(.data, .n = .n)
  } else if (.method == "sample") {
    # choose the repertoire with the smallest amount of clonotypes
    if (is.na(.n)) {
      rows_vec <- unlist(sapply(.data, tally))
      i <- which.min(rows_vec)
      .n <- rows_vec[i]
    }

    sample_clonotypes(.data, .n = .n, .prob = .prob)
  }
}


resample_col <- function(.data, .n) {
  # col_vec <- collect(select(.data, IMMCOL$prop), n = Inf)[[1]]
  col_vec <- collect(select(.data, IMMCOL$count), n = Inf)[[1]]
  col_vec <- col_vec / sum(col_vec)

  rmultinom(1, .n, col_vec)
}

resample <- function(.data, .n) {
  if (has_class(.data, "list")) {
    if (length(.n) != length(.data)) {
      .n <- c(.n, rep.int(-1, length(.data) - length(.n)))
    }
    return(lapply(.data, resample, .n = .n))
  }
  new_col <- resample_col(.data, .n)

  .data <- collect(.data, n = Inf)[new_col != 0, ]
  new_col <- new_col[new_col != 0]
  .data[[IMMCOL$count]] <- new_col
  .data[[IMMCOL$prop]] <- .data[[IMMCOL$count]] / sum(.data[[IMMCOL$count]])
  .data[order(.data[[IMMCOL$prop]], decreasing = TRUE), ]
}

downsample_col <- function(.data, .n) {
  read_vec <- collect(select(.data, IMMCOL$count), n = Inf)[[1]]
  read_indices <- rep(0, sum(read_vec))
  # Rcpp::cppFunction(
  #   '
  #   NumericVector fill_vec(NumericVector read_vec, NumericVector read_indices) {
  #     int dummy = 0;
  #     for (int i = 0; i < read_vec.size(); i++) {
  #     for (int j = dummy; j < (read_vec[i] + dummy); j++) {
  #     read_indices[j] = i;
  #     }
  #     dummy = dummy + read_vec[i];
  #     }
  #     return read_indices;
  #   }
  #   '
  # )
  read_indices <- fill_vec(read_vec, read_indices)
  new_counts <- sample(read_indices, .n)
  new_reads <- rep(0, length(read_vec))
  # Rcpp::cppFunction(
  #   '
  #   NumericVector fill_reads(NumericVector new_reads, NumericVector new_counts) {
  #     for (int i = 0; i < new_counts.size(); i++) {
  #     new_reads[new_counts[i]] = new_reads[new_counts[i]] + 1;
  #     }
  #     return new_reads;
  #   }
  #   '
  # )
  fill_reads(new_reads, new_counts)
}

downsample <- function(.data, .n) {
  if (has_class(.data, "list")) {
    return(lapply(.data, downsample, .n = .n))
  }

  new_col <- downsample_col(.data, .n)

  .data <- collect(.data, n = Inf)[new_col > 0, ]
  new_col <- new_col[new_col > 0]
  .data[[IMMCOL$count]] <- new_col
  .data[[IMMCOL$prop]] <- new_col / sum(new_col)
  .data[order(.data[[IMMCOL$prop]], decreasing = TRUE), ]
}

sample_clonotypes <- function(.data, .n, .prob) {
  if (has_class(.data, "list")) {
    return(lapply(.data, sample_clonotypes, .n = .n, .prob = .prob))
  }

  # ToDo: figure out how to pass a string variable to sample_n weights
  .data <- .data %>% collect(n = Inf)

  weights <- NULL
  if (.prob) {
    weights <- .data[[IMMCOL$count]]
  }

  indices <- sample.int(nrow(.data), .n, replace = FALSE, prob = weights)
  .data <- .data[indices, ]

  .data[[IMMCOL$prop]] <- .data[[IMMCOL$count]] / sum(.data[[IMMCOL$count]])
  .data[order(.data[[IMMCOL$prop]], decreasing = TRUE), ]
}
