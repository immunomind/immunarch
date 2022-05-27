#' Information measures
#'
#' @concept utility_public
#'
#' @aliases entropy kl_div js_div cross_entropy
#'
#' @description Compute information-based estimates and distances.
#'
#' @usage
#' entropy(.data, .base = 2, .norm = FALSE, .do.norm = NA, .laplace = 1e-12)
#'
#' kl_div(.alpha, .beta, .base = 2, .do.norm = NA, .laplace = 1e-12)
#'
#' js_div(.alpha, .beta, .base = 2, .do.norm = NA, .laplace = 1e-12, .norm.entropy = FALSE)
#'
#' cross_entropy(.alpha, .beta, .base = 2, .do.norm = NA,
#'               .laplace = 1e-12, .norm.entropy = FALSE)
#'
#' @param .data Numeric vector. Any distribution.
#' @param .alpha Numeric vector. A distribution of some random value.
#' @param .beta Numeric vector. A distribution of some random value.
#' @param .base Numeric. A base of logarithm.
#' @param .norm Logical. If TRUE then normalises the entropy by the maximal value of the entropy.
#' @param .do.norm If TRUE then normalises the input distributions to make them sum up to 1.
#' @param .laplace Numeric. A value for the laplace correction.
#' @param .norm.entropy Logical. If TRUE then normalises the resulting value by the average entropy of input distributions.
#'
#' @return
#' A numeric value.
#'
#' @examples
#' P <- abs(rnorm(10))
#' Q <- abs(rnorm(10))
#' entropy(P)
#' kl_div(P, Q)
#' js_div(P, Q)
#' cross_entropy(P, Q)
#' @export entropy kl_div js_div cross_entropy
entropy <- function(.data, .base = 2, .norm = FALSE, .do.norm = NA, .laplace = 1e-12) {
  .data <- check_distribution(.data, .do.norm, .laplace, .warn.zero = TRUE)
  res <- -sum(.data * log(.data, base = .base))
  if (.norm) {
    res / log(length(.data), base = .base)
  } else {
    res
  }
}

kl_div <- function(.alpha, .beta, .base = 2, .do.norm = NA, .laplace = 1e-12) {
  .alpha <- check_distribution(.alpha, .do.norm, .laplace, .warn.zero = TRUE)
  .beta <- check_distribution(.beta, .do.norm, .laplace, .warn.zero = TRUE)
  sum(log(.alpha / .beta, base = .base) * .alpha)
}

js_div <- function(.alpha, .beta, .base = 2, .do.norm = NA, .laplace = 1e-12, .norm.entropy = FALSE) {
  .alpha <- check_distribution(.alpha, .do.norm, .laplace, .warn.zero = TRUE)
  .beta <- check_distribution(.beta, .do.norm, .laplace, .warn.zero = TRUE)
  nrm <- if (.norm.entropy) 0.5 * (entropy(.alpha, .base, FALSE, .do.norm, .laplace) + entropy(.beta, .base, FALSE, .do.norm, .laplace)) else 1
  M <- (.alpha + .beta) / 2
  0.5 * (kl_div(.alpha, M, .base, FALSE) + kl_div(.beta, M, .base, FALSE)) / nrm
}

cross_entropy <- function(.alpha, .beta, .base = 2, .do.norm = NA, .laplace = 1e-12, .norm.entropy = FALSE) {
  .alpha <- check_distribution(.alpha, .do.norm, .laplace, .warn.zero = TRUE)
  .beta <- check_distribution(.beta, .do.norm, .laplace, .warn.zero = TRUE)
  -sum(log(.beta, base = .base) * .alpha)
}
