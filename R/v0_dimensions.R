default_scale_fun <- function(x) {
  # if (sum(x < 1, na.rm = TRUE) == 0) {
  #   -log(x)
  # } else {
  #   x
  # }

  x
}


#' Dimensionality reduction
#'
#' @importFrom stats prcomp
#'
#' @concept post_analysis
#'
#' @aliases immunr_pca immunr_mds immunr_tsne
#'
#' @description
#'
#' `r lifecycle::badge('deprecated')`
#'
#' Collects a set of principal variables, reducing the number of not important variables
#' to analyse. Dimensionality reduction makes data analysis algorithms work faster and
#' sometimes more accurate, since it also reduces noise in the data. Currently available
#' methods are:
#'
#' - `immunr_pca` performs PCA (Principal Component Analysis) using [prcomp];
#'
#' - `immunr_mds` performs MDS (Multi-Dimensional Scaling) using isoMODS from MASS package.
#'
#' - `immunr_tsne` performs tSNE (t-Distributed Stochastic Neighbour Embedding) using Rtsne Rtsne package.
#'
#' @usage
#'
#' immunr_pca(.data, .scale = default_scale_fun, .raw = TRUE, .orig = FALSE, .dist = FALSE)
#'
#' immunr_mds(.data, .scale = default_scale_fun, .raw = TRUE, .orig = FALSE, .dist = TRUE)
#'
#' immunr_tsne(.data, .perp = 1, .dist = TRUE, ...)
#'
#' @param .data A matrix or a data frame with features, distance matrix or output from [repOverlapAnalysis] or [geneUsageAnalysis] functions.
#'
#' @param .scale A function to apply to your data before passing it to any of
#' dimensionality reduction algorithms. There is no scaling by default.
#'
#' @param .perp The perplexity parameter for Rtsne. Specifies the number
#' of neighbors each data point must have in the resulting plot.
#'
#' @param .raw If TRUE then returns the non-processed output from dimensionality reduction
#' algorithms. Pass FALSE if you want to visualise results.
#'
#' @param .orig If TRUE then returns the original result from algorithms. Pass FALSE
#' if you want to visualise results.
#'
#' @param .dist If TRUE then assumes that ".data" is a distance matrix.
#'
#' @param ... Other parameters passed to Rtsne.
#'
#' @return
#' `immunr_pca` - an output from [prcomp].
#'
#' `immunr_mds` - an output from isoMDS.
#'
#' `immunr_tsne` - an output from Rtsne.
#'
#' @seealso [vis.immunr_pca] for visualisations.
#'
#' @examples
#' data(immdata)
#' gu <- geneUsage(immdata$data)
#' gu[is.na(gu)] <- 0
#' gu <- t(as.matrix(gu[, -1]))
#' immunr_pca(gu)
#' immunr_mds(dist(gu))
#' immunr_tsne(dist(gu))
#' @export immunr_pca immunr_mds immunr_tsne
immunr_pca <- function(.data, .scale = default_scale_fun, .raw = TRUE, .orig = FALSE, .dist = FALSE) {
  if (.dist) {
    res <- cmdscale(as.dist(.scale(.data)), list. = TRUE)

    if (!.orig) {
      res$scale <- FALSE
      res$center <- c(0, 0, 0, 0)
      res$x <- res$points
      res$rotation <- res$points
      res$sdev <- c(0, 0, 0, 0)
      colnames(res$x) <- c("DimI", "DimII")
      res <- add_class(res, "prcomp")
    }
  } else {
    res <- prcomp(.scale(.data), scale. = TRUE)

    if (!.raw) {
      res <- data.frame(res$x)[c(1, 2)]
      colnames(res) <- c("V1", "V2")
      res$Sample <- row.names(res)
    }
  }

  add_class(res, "immunr_pca")
}

immunr_mds <- function(.data, .scale = default_scale_fun, .raw = TRUE, .orig = FALSE, .dist = TRUE) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for this function. Please install it first via `pak::pkg_install('MASS')`", call. = FALSE)
  }

  if (.dist) {
    .data <- as.dist(.scale(.data))
  } else {
    .data <- .scale(.data)
  }
  res <- MASS::isoMDS(.data, k = 2, trace = FALSE)
  if (!.raw) {
    res <- data.frame(res$points)
    colnames(res) <- c("V1", "V2")
    res$Sample <- row.names(res)
  } else {
    if (!.orig) {
      # Dirty hack to make factoextra work with MDS objects
      res$scale <- FALSE
      res$center <- c(0, 0, 0, 0)
      res$x <- res$points
      res$rotation <- res$points
      res$sdev <- c(0, 0, 0, 0)
      colnames(res$x) <- c("DimI", "DimII")
      res <- add_class(res, "prcomp")
    }
  }
  add_class(res, "immunr_mds")
}

immunr_tsne <- function(.data, .perp = 1, .dist = TRUE, ...) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required for this function. Please install it first via `pak::pkg_install('Rtsne')`", call. = FALSE)
  }

  if (.dist) {
    data_proc <- as.dist(.data)
  } else {
    data_proc <- .data
  }
  res <- Rtsne::Rtsne(data_proc, perplexity = .perp, is_distance = .dist, ...)$Y
  row.names(res) <- row.names(.data)
  colnames(res) <- c("DimI", "DimII")
  add_class(res, "immunr_tsne")
}
