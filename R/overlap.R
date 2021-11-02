#' Main function for public clonotype statistics calculations
#'
#' @concept overlap
#'
#' @importFrom data.table data.table setDT set
#' @importFrom magrittr "%>%"
#' @importFrom dtplyr lazy_dt
#'
#' @description The \code{repOverlap} function is designed to analyse the overlap between
#' two or more repertoires. It contains a number of methods to compare immune receptor
#' sequences that are shared between individuals.
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
#' @param .method A string that specifies the method of analysis or a combination of
#' methods. The \code{repOverlap} function supports following basic methods:
#' "public", "overlap", "jaccard", "tversky", "cosine", "morisita".
#'
#' @param .col A string that specifies the column(s) to be processed. Pass one of the
#' following strings, separated by the plus sign: "nt" for nucleotide sequences,
#' "aa" for amino acid sequences, "v" for V gene segments, "j" for J gene segments. E.g.,
#' pass "aa+v" to compute overlaps on CDR3 amino acid sequences paired with V gene segments, i.e.,
#' in this case a unique clonotype is a pair of CDR3 amino acid and V gene segment.
#' Clonal counts of equal clonotypes will be summed up.
#'
#' @param .a,.b Alpha and beta parameters for Tversky Index. Default values give
#' the Jaccard index measure.
#'
#' @param .verbose if TRUE then output the progress.
#'
#' @param .step Either an integer or a numeric vector.
#'
#' In the first case, the integer defines the step of incremental overlap.
#'
#' In the second case, the vector encodes all repertoire sampling depths.
#'
#' @param .n.steps Something. Skipped if ".step" is a numeric vector.
#'
#' @param .downsample If TRUE then perform downsampling to N clonotypes at each step instead of choosing the
#' top N clonotypes in incremental overlaps. Change nothing for conventional methods.
#'
#' @param .bootstrap Pass NA to turn off any bootstrapping, pass a number to perform bootstrapping with this number of tries.
#'
#' @param .verbose.inc Logical. If TRUE then show output from the computation process.
#'
#' @param .force.matrix Logical. If TRUE than always force the matrix output even in case of two input repertoires.
#'
#' @details "public" and "shared" are synonyms that exist
#' for the convenience of researchers.
#'
#' The "overlap" coefficient is a similarity measure that measures the overlap between two finite sets.
#'
#' The "jaccard" index is conceptually a percentage of how many objects two sets
#' have in common out of how many objects they have total.
#'
#' The "tversky" index is an asymmetric similarity measure on sets that
#' compares a variant to a prototype.
#'
#' The "cosine" index is a measure of similarity between two non-zero vectors
#' of an inner product space that measures the cosine of the angle between them.
#'
#' The "morisita" index measures how many times it is more likely to randomly
#' select two sampled points from the same quadrat (the dataset is covered by a
#' regular grid of changing size) then it would be in the case of a random
#' distribution generated from a Poisson process. Duplicate objects are merged
#' with their counts are summed up.
#'
#' @return
#' In most cases the return value is a matrix with overlap values for each pair of repertoires.
#'
#' If only two repertoires were provided, return value is single numeric value.
#'
#' If one of the incremental method is chosen, return list of overlap matrix.
#'
#' @seealso \link{inc_overlap}, \link{vis}
#'
#' @examples
#' data(immdata)
#'
#' # Make data smaller for testing purposes
#' immdata$data <- top(immdata$data, 4000)
#'
#' ov <- repOverlap(immdata$data, .verbose = FALSE)
#' vis(ov)
#'
#' ov <- repOverlap(immdata$data, "jaccard", .verbose = FALSE)
#' vis(ov, "heatmap2")
#' @export repOverlap
repOverlap <- function(.data,
                       .method = c("public", "overlap", "jaccard", "tversky", "cosine", "morisita", "inc+public", "inc+morisita"),
                       .col = "aa",
                       .a = .5,
                       .b = .5,
                       .verbose = TRUE,
                       .step = 1000,
                       .n.steps = 10,
                       .downsample = FALSE,
                       .bootstrap = NA,
                       .verbose.inc = NA,
                       .force.matrix = FALSE) {
  .validate_repertoires_data(.data)

  .method <- .method[1]

  if (stringr::str_starts(.method, "inc")) {
    if (is.na(.verbose.inc)) {
      .verbose.inc <- TRUE
      .verbose <- FALSE
    }

    .method <- substr(.method, 5, nchar(.method))

    inc_overlap(.data, repOverlap,
      .step = .step, .n.steps = .n.steps,
      .verbose.inc = .verbose.inc, .downsample = .downsample,
      .bootstrap = .bootstrap,
      .method = .method, .col = .col,
      .a = .a, .b = .b, .verbose = .verbose
    )
  } else {
    .col <- process_col_argument(.col)

    if (.method == "cosine" || .method == "morisita") {
      .col <- c(.col, IMMCOL$count)
    }

    for (i in 1:length(.data)) {
      if (has_class(.data[[i]], "data.table")) {
        .data[[i]] <- .data[[i]] %>%
          lazy_dt() %>%
          select(.col) %>%
          collect(n = Inf)
      } else {
        .data[[i]] <- .data[[i]] %>%
          select(.col) %>%
          collect(n = Inf)
      }
    }

    res <- switch(.method,
      shared = num_shared_clonotypes(.data),
      public = num_shared_clonotypes(.data),
      overlap = apply_symm(.data, overlap_coef, .diag = NA, .verbose = .verbose),
      jaccard = apply_symm(.data, jaccard_index, .diag = NA, .verbose = .verbose),
      tversky = apply_symm(.data, tversky_index, .a = .a, .b = .b, .diag = NA, .verbose = .verbose),
      cosine = apply_symm(.data, cosine_sim, .quant = IMMCOL$count, .diag = NA, .verbose = .verbose),
      morisita = apply_symm(.data, morisita_index, .quant = IMMCOL$count, .diag = NA, .verbose = .verbose),
      stop("You entered the wrong method! Please, try again.")
    )

    if (length(.data) == 2) {
      if (.force.matrix) {
        res
      } else {
        res <- res[1, 2]
      }
    }

    add_class(res, "immunr_ov_matrix")
  }
}


num_shared_clonotypes <- function(.data) {
  res <- matrix(0, length(.data), length(.data))

  for (i in 1:(length(.data) - 1)) {
    data_i <- collect(.data[[i]], n = Inf)
    for (j in (i + 1):length(.data)) {
      data_j <- collect(.data[[j]], n = Inf)
      res[i, j] <- nrow(dplyr::intersect(data_i, data_j))
      res[j, i] <- res[i, j]
    }
  }

  diag(res) <- NA

  row.names(res) <- names(.data)
  colnames(res) <- names(.data)

  res
}

overlap_coef <- function(.x, .y) {
  UseMethod("overlap_coef")
}

overlap_coef.default <- function(.x, .y) {
  .x <- collect(.x, n = Inf)
  .y <- collect(.y, n = Inf)
  nrow(dplyr::intersect(.x, .y)) / min(nrow(.x), nrow(.y))
}

overlap_coef.character <- function(.x, .y) {
  length(dplyr::intersect(.x, .y)) / min(length(.x), length(.y))
}


jaccard_index <- function(.x, .y) {
  UseMethod("jaccard_index")
}

jaccard_index.default <- function(.x, .y) {
  .x <- collect(.x, n = Inf)
  .y <- collect(.y, n = Inf)
  intersection <- nrow(dplyr::intersect(.x, .y))
  intersection / (nrow(.x) + nrow(.y) - intersection)
}

jaccard_index.character <- function(.x, .y) {
  intersection <- length(dplyr::intersect(.x, .y))
  intersection / (length(.x) + length(.y) - intersection)
}

tversky_index <- function(.x, .y, .a = .5, .b = .5) {
  UseMethod("tversky_index")
}

tversky_index.default <- function(.x, .y, .a = .5, .b = .5) {
  .x <- collect(.x, n = Inf)
  .y <- collect(.y, n = Inf)
  intersection <- nrow(dplyr::intersect(.x, .y))
  intersection / (.a * nrow(dplyr::setdiff(.x, .y)) + .b * nrow(dplyr::setdiff(.y, .x)) + intersection)
}

tversky_index.character <- function(.x, .y, .a = .5, .b = .5) {
  intersection <- length(dplyr::intersect(.x, .y))
  intersection / (.a * length(dplyr::setdiff(.x, .y)) + .b * length(dplyr::setdiff(.y, .x)) + intersection)
}

cosine_sim <- function(.x, .y, .quant) {
  UseMethod("cosine_sim")
}

cosine_sim.default <- function(.x, .y, .quant) {
  .x <- collect(.x, n = Inf)
  .y <- collect(.y, n = Inf)
  col_name <- colnames(.x)[colnames(.x) != .quant]
  joined_set <- full_join(.x, .y, by = col_name, suffix = c("_a", "_b"))
  alpha_col <- paste0(.quant, "_a")
  beta_col <- paste0(.quant, "_b")
  new_set <- collect(select(joined_set, c(alpha_col, beta_col)))
  first_col <- new_set[, 1]
  second_col <- new_set[, 2]
  first_col[is.na(first_col)] <- 0
  second_col[is.na(second_col)] <- 0
  sum(first_col * second_col) / (sqrt(sum(first_col * first_col)) * sqrt(sum(second_col * second_col)))
}

cosine_sim.numeric <- function(.x, .y, .quant) {
  df <- rbind(.x, .y)
  sum(.x * .y) / (sqrt(rowSums(df^2))[1] * sqrt(rowSums(df^2))[2])[[1]]
}

morisita_index <- function(.x, .y, .quant) {
  .quant <- .quant[1]
  # assumption #1: both objects have same column names
  # ToDo: make an assert to check this
  alpha <- collect(.x, n = Inf)
  beta <- collect(.y, n = Inf)

  alpha <- collect(.x, n = Inf)
  beta <- collect(.y, n = Inf)
  setDT(alpha)
  setDT(beta)

  obj_columns <- setdiff(names(alpha), .quant)

  alpha <- alpha[, sum(get(.quant)), by = obj_columns]
  beta <- beta[, sum(get(.quant)), by = obj_columns]
  joined <- merge(alpha, beta, by = obj_columns, all = TRUE)

  col_a <- ncol(joined) - 1
  col_b <- ncol(joined)
  set(joined, which(is.na(joined[[col_a]])), as.integer(col_a), 0)
  set(joined, which(is.na(joined[[col_b]])), as.integer(col_b), 0)

  X <- sum(joined[[col_a]])
  Y <- sum(joined[[col_b]])
  sum_alp_bet <- sum(joined[[col_a]] * joined[[col_b]])
  sum_alp_sq <- sum(joined[[col_a]]^2)
  sum_bet_sq <- sum(joined[[col_b]]^2)

  2 * sum_alp_bet / ((sum_alp_sq / (X^2) + sum_bet_sq / (Y^2)) * X * Y)
}

horn_index <- function(.x, .y) {
  #
  # ToDo: rewrite
  #
  .do.unique <- TRUE
  .x[, 1] <- as.character(.x[, 1])
  .y[, 1] <- as.character(.y[, 1])
  colnames(.x) <- c("Species", "Count")
  colnames(.y) <- c("Species", "Count")
  if (.do.unique) {
    .x <- .x[!duplicated(.x[, 1]), ]
    .y <- .y[!duplicated(.y[, 1]), ]
  }
  .x[, 2] <- as.numeric(.x[, 2]) / sum(.x[, 2])
  .y[, 2] <- as.numeric(.y[, 2]) / sum(.y[, 2])
  merged <- merge(.x, .y, by = "Species", all = TRUE)
  merged[is.na(merged)] <- 0
  rel.12 <- merged[, 2] / merged[, 3]
  rel.12[merged[, 3] == 0] <- 0
  rel.21 <- merged[, 3] / merged[, 2]
  rel.21[merged[, 2] == 0] <- 0
  1 / log(2) * sum(merged[, 2] / 2 * log(1 + rel.21) + merged[, 3] / 2 * log(1 + rel.12))
}


#' Incremental counting of repertoire similarity
#'
#' @concept overlap
#'
#' @description Like in paper https://www.pnas.org/content/111/16/5980 (Fig. 4).
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
#' @param .fun Function to compute overlaps. e.g., \code{morisita_index}.
#'
#' @param .step Either an integer or a numeric vector.
#'
#' In the first case, the integer defines the step of incremental overlap.
#'
#' In the second case, the vector encodes all repertoire sampling depths.
#'
#' @param .n.steps Integer. Number of steps if \code{.step} is a single integer.
#' Skipped if ".step" is a numeric vector.
#'
#' @param .downsample If TRUE then perform downsampling to N clonotypes at each step instead of choosing the
#' top N clonotypes.
#'
#' @param .bootstrap Pass NA to turn off any bootstrapping, pass a number to perform bootstrapping with this number of tries.
#'
#' @param .verbose.inc Logical. If TRUE then show output from the computation process.
#'
#' @param ... Other arguments passed to \code{.fun}.
#'
#' @return
#' List with overlap matrices.
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data, "inc+overlap", .step = 100, .verbose.inc = FALSE, .verbose = FALSE)
#' vis(ov)
#' @export inc_overlap
inc_overlap <- function(.data, .fun, .step = 1000, .n.steps = 10, .downsample = FALSE, .bootstrap = NA, .verbose.inc = TRUE, ...) {
  .n.steps <- as.integer(.n.steps)
  if (.n.steps < 1) {
    stop("Error: please provide the number of steps greater than 1.")
  }

  if (length(.step) == 1) {
    step_vec <- seq(.step, min(sapply(.data, function(d) {
      if (has_class(d, "data.table")) {
        d %>%
          lazy_dt() %>%
          select(IMMCOL$count) %>%
          collect(n = Inf) %>%
          nrow()
      } else {
        d %>%
          select(IMMCOL$count) %>%
          collect(n = Inf) %>%
          nrow()
      }
    })), .step)
    step_vec <- step_vec[1:min(length(step_vec), .n.steps)]
  } else {
    step_vec <- .step
  }

  args <- list(...)
  if (!(".verbose" %in% names(args))) {
    args[[".verbose"]] <- FALSE
  }

  if (.verbose.inc) {
    if (is.na(.bootstrap)) {
      pb <- set_pb(length(step_vec))
    } else {
      pb <- set_pb(length(step_vec) * .bootstrap)
    }
  }

  res <- lapply(step_vec, function(i) {
    if (is.na(.bootstrap)) {
      if (.downsample) {
        fun_res <- .fun(repSample(.data, .method = "sample", .n = i, .prob = TRUE), ...)
      } else {
        fun_res <- .fun(lapply(.data, head, i), ...)
      }

      if (.verbose.inc) {
        add_pb(pb)
      }
      fun_res
    } else {
      fun_res <- list()
      for (i_boot in 1:.bootstrap) {
        fun_res[[i_boot]] <- .fun(repSample(.data, .method = "sample", .n = i, .prob = TRUE), ...)

        if (.verbose.inc) {
          add_pb(pb)
        }
      }

      fun_res
    }
  })

  if (.verbose.inc) {
    close(pb)
  }

  names(res) <- step_vec

  res <- add_class(res, "immunr_inc_overlap")
  if (!is.na(.bootstrap)) {
    attr(res, "bootstrap") <- .bootstrap
  }

  res
}
