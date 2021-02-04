#' Calculate the kmer tatistics of immune repertoires
#'
#' @importFrom data.table data.table
#' @importFrom dplyr as_tibble
#'
#' @concept kmers
#'
#' @aliases getKmers get.kmers makeKmerTable
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
#' @param .k Integer. Length of kmers.
#' @param .col Character. Which column to use, pass "aa" (by default) for CDR3 amino acid sequence,
#' pass "nt" for CDR3 nucleotide sequences.
#' @param .coding Logical. If TRUE (by default) then remove all non-coding sequences from input data first.
#'
#' @return
#' Data frame with two columns (kmers and their counts).
#'
#' @examples
#' data(immdata)
#' kmers <- getKmers(immdata$data[[1]], 5)
#' kmers %>% vis()
#' @export getKmers
getKmers <- function(.data, .k, .col = c("aa", "nt"), .coding = TRUE) {
  seq_col <- switch_type(.col[1])

  if (.coding) {
    .data <- coding(.data)
  }

  if (has_class(.data, "list")) {
    res <- data.table(split_to_kmers(collect(select(.data[[1]], seq_col), n = Inf)[[1]], .k = .k))
    colnames(res)[2] <- names(.data)[1]
    for (i in 2:length(.data)) {
      tmp <- data.table(split_to_kmers(collect(select(.data[[i]], seq_col), n = Inf)[[1]], .k = .k))
      res <- merge(res, tmp, by = "Kmer", all = TRUE, suffixes = c("", ""))
      colnames(res)[ncol(res)] <- names(.data)[i]
    }
  } else {
    res <- split_to_kmers(collect(select(.data, seq_col), n = Inf)[[1]], .k = .k)
  }

  add_class(as_tibble(as.data.frame(res)), "immunr_kmer_table")
}


# WIP
# kmerAnalysis <- function (.data) {
#   UseMethod("kmerAnalysis")
# }
#
# kmerAnalysis.default <- function (.data, .k, .method = c("profile", "gibbs"), .col = c("aa", "nuc"), .ic = FALSE, .remove.stop = TRUE, .use.counts = FALSE) {
#   kmers = getKmers(.data, .k = .k, .col = .col)
#
#   kmerAnalysis(kmers, .method = .method, .use.counts = .use.counts, .ic = .ic, .remove.stop = .remove.stop)
# }
#
# kmerAnalysis.immunr_kmer_table <- function (.data, .method = c("profile", "gibbs"), ...) {
#   .method = .method[1]
#
#   if (.method == "profile") {
#     stop("Unimplemented")
#   } else if (.method == "gibbs") {
#     stop("Unimplemented")
#   } else {
#     stop("Unknown method")
#   }
#
#   res
# }


#' Analysis immune repertoire kmer statistics: sequence profiles, etc.
#'
#' @importFrom dplyr as_tibble
#'
#' @concept kmers
#'
#' @aliases split_to_kmers kmer_profile
#'
#' @usage
#' split_to_kmers(.data, .k)
#'
#' kmer_profile(.data, .method = c("freq", "prob", "wei", "self"), .remove.stop = TRUE)
#'
#' @param .data Character vector or the output from \code{getKmers}.
#' @param .k Integer. Size of kmers.
#' @param .method Character vector of length one. If "freq" then return a position frequency matrix (PFM) -
#' a matrix with occurences of each amino acid in each position.
#'
#' If "prob" then return a position probability matrix (PPM) - a matrix with probabilities of occurences of
#' each amino acid in each position. This is a traditional representation of sequence motifs.
#'
#' If "wei" then return a position weight matrix (PWM) - a matrix with log likelihoods of PPM elements.
#'
#' If "self" then return a matrix with self-information of elements in PWM.
#'
#' For more information see https://en.wikipedia.org/wiki/Position_weight_matrix.
#' @param .remove.stop Logical. If TRUE (by default) remove stop codons.
#'
#' @return
#' \code{split_to_kmers} - Data frame with two columns (kmers and their counts).
#'
#' \code{kmer_profile} - a matrix with per-position amino acid statistics.
#'
#' @examples
#' data(immdata)
#' kmers <- getKmers(immdata$data[[1]], 5)
#' kmer_profile(kmers) %>% vis()
#' @export split_to_kmers kmer_profile
split_to_kmers <- function(.data, .k) {
  max_len <- max(nchar(.data))
  tab <- table(unlist(lapply(1:(max_len - .k + 1), function(i) substr(.data, i, i + .k - 1))))
  tab <- tab[nchar(names(tab)) == .k]
  tab <- as_tibble(as.data.frame(tab, stringsAsFactors = FALSE))
  names(tab) <- c("Kmer", "Count")
  add_class(tab, "immunr_kmers")
}

kmer_profile <- function(.data, .method = c("freq", "prob", "wei", "self"), .remove.stop = TRUE) {
  # ToDo: add a support for nucletide profiles

  .method <- .method[1]
  if (!(.method %in% c("freq", "prob", "wei", "self"))) {
    stop('Error: wrong name for the .method argument. Please provide one of the following: "freq", "prob", "wei" or "self".')
  }

  if (has_class(.data, "immunr_kmer_table")) {
    # ToDo: make it work for any number of samples
    if (ncol(.data) > 2) {
      stop("Error: currently we are not supporting sequence profiles from more than one sample. Will support it soon!")
    }
    seq_vec <- .data$Kmer
    cnt_vec <- .data$Count
  } else if (length(table(nchar(.data))) > 1) {
    stop("Not all kmers in the input data have the same length.")
  } else {
    seq_vec <- .data
    cnt_vec <- rep.int(1, length(.data))
  }

  k <- nchar(seq_vec[1])
  aas <- sort(unique(AA_TABLE))
  if (.remove.stop) {
    aas <- aas[3:length(aas)]
    cnt_vec <- cnt_vec[grep("[*~]", seq_vec, invert = TRUE)]
    seq_vec <- seq_vec[grep("[*~]", seq_vec, invert = TRUE)]
  }

  res <- matrix(0, length(aas), k)
  row.names(res) <- aas
  bad_symbols <- c()
  for (i in 1:k) {
    tab <- tapply(cnt_vec, substr(seq_vec, i, i), sum, simplify = FALSE)
    for (aa in names(tab)) {
      if (aa %in% row.names(res)) {
        res[aa, i] <- res[aa, i] + tab[[aa]]
      } else {
        bad_symbols <- c(bad_symbols, aa)
      }
    }

    if (.method == "freq") {
      # nothing to do here :(
    } else if (.method == "prob") {
      res[, i] <- res[, i] / sum(res[, i])
    } else if (.method == "wei") {
      res[, i] <- res[, i] / nrow(res)
      res[, i] <- log2(res[, i])
    } else {
      res[, i] <- res[, i] / sum(res[, i])
      res[, i] <- -res[, i] * log2(res[, i])
    }
  }

  if (length(bad_symbols)) {
    warning(
      "Warning: removed ", length(bad_symbols), " non-amino acid symbol(s): ", unique(bad_symbols),
      "\nPlease make sure your data doesn't have them in the future."
    )
  }

  if (.method == "freq") {
    add_class(res, "immunr_kmer_profile_pfm")
  } else if (.method == "prob") {
    add_class(res, "immunr_kmer_profile_ppm")
  } else if (.method == "wei") {
    add_class(res, "immunr_kmer_profile_pwm")
  } else {
    add_class(res, "immunr_kmer_profile_self")
  }
}


#######
# WIP #
#######
# gibbs_sampling <- function (.data, .motif.len = 5, .niter = 500) {
#   .score <- function (.seq, .i, .prof, .background) {
#     kmer_aa = strsplit(substr(.seq, seq_i, seq_i + .motif.len - 1), "")[[1]]
#     prod(sapply(1:.motif.len, function (kmer_pos) {
#       sc = .prof[kmer_aa[kmer_pos], kmer_pos] / .background[kmer_aa[kmer_pos]]
#       if (is.nan(sc)) { sc = 0 }
#       sc
#     }))
#   }
#
#   cat("Removed", sum(nchar(.data) < .motif.len), "sequences with the length less than the length of motifs.\n")
#   seq_vec = .data[nchar(.data) >= .motif.len]
#   background = table(unlist(strsplit(seq_vec, "")))
#   background = background / sum(background)
#
#   # Vector of scores for each position in the each input sequence
#   score_vec = lapply(seq_vec, function (seq_x) rep(1, nchar(seq_x) - .motif.len + 1) )
#   start_pos = sapply(nchar(seq_vec), function (max_pos) sample(1:(max_pos - .motif.len + 1), 1))
#
#   # In the loop:
#   pb = set_pb(.niter)
#   for (iter in 1:.niter) {
#     # Get random kmers
#     prev_start_pos = start_pos
#     start_pos = sapply(nchar(seq_vec), function (max_pos) sample(1:(max_pos - .motif.len + 1), 1))
#
#     for (out_kmer_i in sample(1:length(seq_vec), length(seq_vec))) {
#       max_pos = nchar(seq_vec[out_kmer_i]) - .motif.len + 1
#       kmers <- substr(seq_vec[-out_kmer_i], start_pos[-out_kmer_i], start_pos[-out_kmer_i] + .motif.len - 1)
#       prof = kmer_profile(kmers[-out_kmer_i])
#       for (seq_i in 1:max_pos) {
#         score_vec[[out_kmer_i]][seq_i] = .score(seq_vec[out_kmer_i], seq_i, prof, background)
#       }
#       if (sum(score_vec[[out_kmer_i]]) != 0) {
#         poses = c(1:max_pos)[!is.na(score_vec[[out_kmer_i]])]
#         start_pos[out_kmer_i] = sample(c(1:max_pos), 1, prob = score_vec[[out_kmer_i]][poses] / sum(score_vec[[out_kmer_i]][poses]))
#       }
#     }
#
#     add_pb(pb)
#
#     if (sum(prev_start_pos != start_pos) == 0) {
#       break
#     }
#   }
#   close(pb)
#
#   data.frame(Motif = substr(seq_vec, start_pos, start_pos + .motif.len - 1),
#              Start = start_pos,
#              Score = sapply(1:length(score_vec), function (i) { score_vec[[i]][start_pos[i]] }), stringsAsFactors = FALSE)
# }
