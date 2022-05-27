if (getRversion() >= "2.15.1") {
  utils::globalVariables("Div.count")
}


#' The main function for immune repertoire diversity estimation
#'
#' @concept diversity
#'
#' @aliases repDiversity chao1 hill_numbers diversity_eco gini_simpson inverse_simpson gini_coef rarefaction
#'
#' @importFrom reshape2 melt
#' @importFrom utils tail
#' @importFrom dplyr mutate group_by_at pull
#' @importFrom stats qnorm
#' @importFrom rlang sym
#'
#' @description
#' This is a utility function to estimate the diversity of species or objects in the given distribution.
#'
#' Note: functions will check if .data is a distribution of a random variable (sum == 1) or not.
#' To force normalisation and / or to prevent this, set .do.norm to TRUE (do normalisation)
#' or FALSE (don't do normalisation), respectively.
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
#' @param .method Picks a method used for estimation out of a following list: chao1,
#' hill, div, gini.simp, inv.simp, gini, raref, d50, dxx.
#' @param .col A string that specifies the column(s) to be processed. Pass one of the
#' following strings, separated by the plus sign: "nt" for nucleotide sequences,
#' "aa" for amino acid sequences, "v" for V gene segments, "j" for J gene segments. E.g.,
#' pass "aa+v" to compute diversity estimations on CDR3 amino acid sequences paired with V gene segments, i.e.,
#' in this case a unique clonotype is a pair of CDR3 amino acid and V gene segment.
#' Clonal counts of equal clonotypes will be summed up.
#' @param .min.q Function calculates several hill numbers. Set the min (default: 1).
#' @param .max.q The max hill number to calculate (default: 5).
#' @param .q q-parameter for the Diversity index.
#' @param .step Rarefaction step's size.
#' @param .quantile Numeric vector with quantiles for confidence intervals.
#' @param .extrapolation An integer. An upper limit for the number of clones to extrapolate to.
#' Pass 0 (zero) to turn extrapolation subroutines off.
#' @param .perc Set the percent to dXX index measurement.
#' @param .norm Normalises rarefaction curves.
#' @param .verbose If TRUE then outputs progress.
#' @param .do.norm One of the three values - NA, TRUE or FALSE. If NA then checks for distrubution (sum(.data) == 1)
#' and normalises if needed with the given laplace correction value. if TRUE then does normalisation and laplace
#' correction. If FALSE then doesn't do neither normalisaton nor laplace correction.
#' @param .laplace A numeric value, which is used as a pseudocount for Laplace
#' smoothing.
#'
#' @return div, gini, gini.simp, inv.simp, raref return numeric vector of length 1
#' with value.
#'
#' chao1 returns 4 values: estimated number of species, standart deviation of
#' this number and two 95% confidence intervals for the species number.
#'
#' hill returns a vector of specified length \code{.max.q - .min.q}
#'
#' @details
#' - True diversity, or the effective number of types, refers to the number
#' of equally-abundant types needed for the average proportional abundance
#' of the types to equal that observed in the dataset of interest
#' where all types may not be equally abundant.
#'
#' - Inverse Simpson index is the effective number of types that is obtained when
#' the weighted arithmetic mean is used to quantify average
#' proportional abundance of types in the dataset of interest.
#'
#' - The Gini coefficient measures the inequality among values
#' of a frequency distribution (for example levels of income). A Gini coefficient of zero
#' expresses perfect equality, where all values are the same (for example, where everyone
#' has the same income). A Gini coefficient of one (or 100 percents ) expresses maximal inequality
#' among values (for example where only one person has all the income).
#'
#' - The Gini-Simpson index is the probability of interspecific encounter, i.e., probability that two entities
#' represent different types.
#'
#' - Chao1 estimator is a nonparameteric asymptotic estimator of species richness (number of species in a population).
#'
#' - Rarefaction is a technique to assess species richness from the results of sampling through extrapolation.
#'
#' - Hill numbers are a mathematically  unified family of  diversity indices (differing among themselves
#' only by an exponent q).
#'
#' - d50 is a recently developed immune diversity estimate. It calculates the minimum number of distinct clonotypes
#' amounting to greater than or equal to 50 percent of a total of sequencing reads obtained
#' following amplification and sequencing
#'
#' - dXX is a similar to d50 index where XX corresponds to desirable percent of total sequencing reads.
#'
#' @return
#' For most methods, if input data is a single immune repertoire, then the function returns a numeric vector
#' with diversity statistics.
#'
#' Otherwise, it returns a numeric matrix with diversity statistics for all input repertoires.
#'
#' For Chao1 the function returns a matrix with diversity estimations.
#'
#' For rarefaction the function returns either a matrix with diversity estimatinos
#' on different step of the simulaiton process or a list with such matrices.
#'
#'
#' @seealso \link{repOverlap}, \link{entropy}, \link{repClonality}
#' Rarefaction wiki
#' \url{https://en.wikipedia.org/wiki/Rarefaction_(ecology)}
#' Hill numbers paper
#' \url{https://www.uvm.edu/~ngotelli/manuscriptpdfs/ChaoHill.pdf}
#' Diversity wiki
#' \url{https://en.wikipedia.org/wiki/Measurement_of_biodiversity}
#'
#'
#' @examples
#' data(immdata)
#'
#' # Make data smaller for testing purposes
#' immdata$data <- top(immdata$data, 4000)
#'
#' # chao1
#' repDiversity(.data = immdata$data, .method = "chao1") %>% vis()
#'
#' # Hill numbers
#' repDiversity(
#'   .data = immdata$data, .method = "hill", .max.q = 6,
#'   .min.q = 1, .do.norm = NA, .laplace = 0
#' ) %>% vis()
#'
#' # diversity
#' repDiversity(.data = immdata$data, .method = "div", .q = 5, .do.norm = NA, .laplace = 0) %>%
#'   vis()
#'
#' # Gini-Simpson
#' repDiversity(.data = immdata$data, .method = "gini.simp", .q = 5, .do.norm = NA, .laplace = 0) %>%
#'   vis()
#'
#' # inverse Simpson
#' repDiversity(.data = immdata$data, .method = "inv.simp", .do.norm = NA, .laplace = 0) %>% vis()
#'
#' # Gini coefficient
#' repDiversity(.data = immdata$data, .method = "gini", .do.norm = NA, .laplace = 0)
#'
#' # d50
#' repDiversity(.data = immdata$data, .method = "d50") %>% vis()
#' @export repDiversity
repDiversity <- function(.data, .method = "chao1", .col = "aa", .max.q = 6, .min.q = 1, .q = 5, .step = NA,
                         .quantile = c(.025, .975), .extrapolation = NA, .perc = 50,
                         .norm = TRUE, .verbose = TRUE, .do.norm = NA, .laplace = 0) {
  .method <- .method[1]

  if (.method == "rarefaction") {
    .method <- "raref"
  }

  if (has_class(.data, "list") && .method != "raref") {
    res <- lapply(.data, repDiversity,
      .method = .method, .col = .col,
      .max.q = .max.q, .min.q = .min.q, .q = .q,
      .step = .step, .quantile = .quantile,
      .extrapolation = .extrapolation, .perc = .perc,
      .verbose = .verbose, .do.norm = .do.norm, .laplace = .laplace
    )
    new_class <- head(class(res[[1]]), 1)
    res <- do.call(rbind, res)
    if (.method == "hill") {
      res <- reshape2::melt(res)
      colnames(res) <- c("Sample", "Q", "Value")
      res$Q <- as.numeric(sapply(res$Q, stringr::str_sub, start = 2))
    } else if (.method %in% c("div", "gini.simp", "inv.simp")) {
      res <- reshape2::melt(res)[c(1, 3)]
      colnames(res) <- c("Sample", "Value")
    } else if (.method == "chao1") {
      colnames(res) <- c("Estimator", "SD", "Conf.95.lo", "Conf.95.hi")
    }

    res <- add_class(res, new_class)
    return(res)
  } else if (.method == "raref") {
    if (!has_class(.data, "list")) {
      .data <- list(Sample = .data)
    }


    .col <- process_col_argument(.col) # ToDo: refactor this and the next branches

    vec <- lapply(.data, function(x) {
      if (has_class(x, "data.table")) {
        x <- x %>% lazy_dt()
      }
      x %>%
        select(.col, IMMCOL$count) %>%
        group_by_at(vars(.col)) %>%
        summarise(Div.count = sum(!!sym(IMMCOL$count))) %>%
        pull(Div.count)
    })
  } else {
    .col <- process_col_argument(.col)

    if (has_class(.data, "data.table")) {
      .data <- .data %>% lazy_dt()
    }

    vec <- .data %>%
      select(.col, IMMCOL$count) %>%
      group_by_at(vars(.col)) %>%
      summarise(Div.count = sum(!!sym(IMMCOL$count))) %>%
      pull(Div.count)
  }

  switch(tolower(.method),
    chao1 = chao1(.data = vec),
    hill = hill_numbers(.data = vec, .max.q = .max.q, .min.q = .min.q, .do.norm = TRUE, .laplace = .laplace),
    div = diversity_eco(.data = vec, .q = .q, .do.norm = TRUE, .laplace = .laplace),
    gini.simp = gini_simpson(.data = vec, .do.norm = TRUE, .laplace = .laplace),
    inv.simp = inverse_simpson(.data = vec, .do.norm = TRUE, .laplace = .laplace),
    gini = gini_coef(.data = vec, .do.norm = TRUE, .laplace = .laplace),
    raref = rarefaction(
      .data = vec, .step = .step, .quantile = .quantile,
      .extrapolation = .extrapolation, .norm = .norm, .verbose = .verbose
    ),
    dxx = dXX(.data = vec, .perc = .perc),
    d50 = dXX(.data = vec, .perc = 50),
    stop("You entered the wrong method! Please, try again.")
  )
}

chao1 <- function(.data) {
  counts <- table(.data)
  e <- NA
  v <- NA
  lo <- NA
  hi <- NA
  n <- sum(.data)
  D <- length(.data)
  f1 <- counts["1"]
  f2 <- counts["2"]
  # f1 == 0 && f2 == 0
  if (is.na(f1) && is.na(f2)) {
    e <- D
    i <- unique(.data)
    v <- sum(sapply(i, function(i) sum(.data == i) * (exp(-i) - exp(-2 * i)))) - (sum(sapply(i, function(i) i * exp(-i) * sum(.data == i))))^2 / n
    P <- sum(sapply(i, function(i) sum(.data == i) * exp(-i) / D))
    lo <- max(D, D / (1 - P) - qnorm(1 - .05 / 2) * sqrt(v) / (1 - P))
    hi <- D / (1 - P) + qnorm(1 - .05 / 2) * sqrt(v) / (1 - P)
  }
  # f1 != 0 && f2 == 0
  else if (is.na(f2)) {
    e <- D + f1 * (f1 - 1) / 2 * (n - 1) / n
    v <- (n - 1) / n * f1 * (f1 - 1) / 2 + ((n - 1) / n)^2 * f1 * (2 * f1 - 1)^2 / 4 - ((n - 1) / n)^2 * f1^4 / 4 / e
    t <- e - D
    K <- exp(qnorm(1 - .05 / 2) * sqrt(log(1 + v / t^2)))
    lo <- D + t / K
    hi <- D + t * K
  }
  # f1 != && f2 != 0
  else {
    const <- (n - 1) / n
    e <- D + f1^2 / (2 * f2) * const
    f12 <- f1 / f2
    v <- f2 * (const * f12^2 / 2 + const^2 * f12^3 + const^2 * f12^4 / 4)
    t <- e - D
    K <- exp(qnorm(.975) * sqrt(log(1 + v / t^2)))
    lo <- D + t / K
    hi <- D + t * K
  }
  res <- c(Estimator = e, SD = sqrt(v), Conf.95.lo = lo, Conf.95.hi = hi)
  add_class(res, "immunr_chao1")
}

hill_numbers <- function(.data, .max.q = 6, .min.q = 1,
                         .do.norm = NA, .laplace = 0) {
  .data <- check_distribution(.data, .do.norm = .do.norm, .laplace)
  if (.min.q < 0) {
    .min.q <- 0
  }
  res <- c()

  for (q in .min.q:.max.q) {
    res <- c(res, diversity_eco(.data, q))
  }

  names(res) <- paste0("q", .min.q:.max.q)
  add_class(res, "immunr_hill")
}

diversity_eco <- function(.data, .q = 5, .do.norm = NA, .laplace = 0) {
  .data <- check_distribution(.data, .do.norm = NA, .laplace = 0)
  if (.q == 0) {
    res <- length(.data)
  } else if (.q == 1) {
    res <- exp(-sum(.data * log(.data)))
  } else if (.q > 1) {
    res <- 1 / (sum(.data^.q)^(1 / (.q - 1)))
  } else {
    res <- NA
  }
  add_class(res, "immunr_div")
}

gini_coef <- function(.data, .do.norm = NA, .laplace = 0) {
  .data <- sort(check_distribution(.data, .do.norm, .laplace, .warn.sum = FALSE))
  n <- length(.data)
  res <- 1 / n * (n + 1 - 2 * sum((n + 1 - 1:n) * .data) / sum(.data))
  add_class(res, "immunr_gini")
}

gini_simpson <- function(.data, .do.norm = NA, .laplace = 0) {
  .data <- check_distribution(.data, .do.norm, .laplace)
  res <- 1 - sum(.data^2)
  add_class(res, "immunr_ginisimp")
}

inverse_simpson <- function(.data, .do.norm = NA, .laplace = 0) {
  .data <- check_distribution(.data, .do.norm, .laplace)
  res <- 1 / sum(.data^2)
  add_class(res, "immunr_invsimp")
}

dXX <- function(.data, .perc = 10) {
  if (has_class(.data, "list")) {
    return(t(sapply(.data, dXX, .perc = .perc)))
  }
  prop <- 0
  n <- 0
  col <- .data
  col <- sort(col, decreasing = TRUE)
  col.sum <- sum(col)
  while (prop < col.sum * (.perc / 100)) {
    n <- n + 1
    prop <- prop + col[n]
  }
  res <- c(Clones = n, Percentage = 100 * signif((prop / col.sum), 3), Clonal.count.prop = n / nrow(.data))
  add_class(res, "immunr_dxx")
}

rarefaction <- function(.data, .step = NA, .quantile = c(.025, .975),
                        .extrapolation = NA, .norm = TRUE, .verbose = TRUE) {
  if (has_class(.data, "numeric") || has_class(.data, "integer")) {
    .data <- list(Sample = .data)
  }

  if (is.na(.step)) {
    .step <- min(sapply(.data, function(x) sum(as.numeric(x)))) %/% 50.
  }

  if (is.na(.extrapolation)) {
    .extrapolation <- max(sapply(.data, function(x) sum(as.numeric(x)))) * 20
  }

  .alpha <- function(n, Xi, m) {
    k <- Xi
    return((1 - m / n)^Xi)
  }

  if (.verbose) {
    pb <- set_pb(sum(sapply(1:length(.data), function(i) {
      bc.vec <- .data[[i]]
      bc.sum <- sum(.data[[i]])
      sizes <- seq(.step, bc.sum, .step)
      if (sizes[length(sizes)] != bc.sum) {
        sizes <- c(sizes, bc.sum)
      }
      length(sizes)
    })))
  }

  muc.list <- lapply(1:length(.data), function(i) {
    Sobs <- length(.data[[i]])
    bc.vec <- .data[[i]]
    Sest <- chao1(bc.vec)
    n <- sum(bc.vec)
    sizes <- seq(.step, n, .step)
    # if (sizes[length(sizes)] != n) {
    #   sizes <- c(sizes, n)
    # }
    counts <- table(bc.vec)
    muc.res <- t(sapply(sizes, function(sz) {
      freqs <- as.numeric(names(counts))

      alphas <- sapply(freqs, function(k) .alpha(n, k, sz))

      # poisson
      Sind <- sum(sapply(1:length(freqs), function(k) (1 - alphas[k]) * counts[k]))
      if (Sest[1] == Sobs) {
        SD <- 0
      } else {
        SD <- sqrt(sum(sapply(1:length(freqs), function(k) (1 - alphas[k])^2 * counts[k])) - Sind^2 / Sest[1])
      }
      t <- Sind - Sobs
      if (t != 0) {
        K <- exp(qnorm(.975) * sqrt(log(1 + (SD / t)^2)))
        lo <- Sobs + t * K
        hi <- Sobs + t / K
      } else {
        lo <- Sind
        hi <- Sind
      }
      res <- c(sz, lo, Sind, hi)
      names(res) <- c("Size", paste0("Q", .quantile[1]), "Mean", paste0("Q", .quantile[2]))
      if (.verbose) add_pb(pb)
      res
    }))

    if (.extrapolation > 0) {
      # sizes <- seq(sum(.data[[i]]), .extrapolation + max(sapply(.data, function (x) sum(x))), .step)
      sizes <- seq(
        tail(seq(.step, sum(.data[[i]]), .step), 1) + .step,
        .extrapolation,
        .step
      )
      if (length(sizes) != 1) {
        ex.res <- t(sapply(sizes, function(sz) {
          f0 <- Sest[1] - Sobs
          f1 <- counts["1"]
          if (is.na(f1) || f0 == 0) {
            Sind <- Sobs
          } else {
            Sind <- Sobs + f0 * (1 - exp(-(sz - n) / n * f1 / f0))
          }
          res <- c(sz, Sind, Sind, Sind)
          names(res) <- c("Size", paste0("Q", .quantile[1]), "Mean", paste0("Q", .quantile[2]))
          if (.verbose) add_pb(pb)
          res
        }))
        df1 <- data.frame(muc.res, Sample = names(.data)[i], Type = "interpolation", stringsAsFactors = FALSE)
        df2 <- data.frame(ex.res, Sample = names(.data)[i], Type = "extrapolation", stringsAsFactors = FALSE)
        muc.res <- rbind(df1, df2)
      } else {
        muc.res <- data.frame(muc.res, Sample = names(.data)[i], Type = "interpolation", stringsAsFactors = FALSE)
      }
    } else {
      muc.res <- data.frame(muc.res, Sample = names(.data)[i], stringsAsFactors = FALSE)
    }

    if (.norm) {
      quantile1 <- paste0("Q", .quantile[1])
      quantile2 <- paste0("Q", .quantile[2])
      muc.res <- muc.res %>%
        mutate(
          Size = Size / n,
          !!quantile1 := !!rlang::sym(quantile1) / Sobs,
          Mean = Mean / Sobs,
          !!quantile2 := !!rlang::sym(quantile2) / Sobs
        )
    }

    muc.res
  })
  if (.verbose) close(pb)

  res <- do.call(rbind, muc.list)

  add_class(res, "immunr_rarefaction")
}
