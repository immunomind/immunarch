#' Clustering of objects or distance matrices
#'
#' @concept post_analysis
#'
#' @aliases immunr_hclust immunr_kmeans immunr_dbscan
#'
#' @importFrom stats kmeans as.dist cmdscale dist
#'
#' @description
#'
#' `r lifecycle::badge('deprecated')`
#'
#' Clusters the data with one of the following methods:
#'
#' - `immunr_hclust` clusters the data using the hierarchical clustering from [hcut][factoextra::hcut];
#'
#' - `immunr_kmeans` clusters the data using the K-means algorithm from [kmeans];
#'
#' - `immunr_dbscan` clusters the data using the DBSCAN algorithm from [dbscan][fpc::dbscan].
#'
#' @usage
#' immunr_hclust(.data, .k = 2, .k.max = nrow(.data) - 1, .method = "complete", .dist = TRUE)
#'
#' immunr_kmeans(.data, .k = 2, .k.max = as.integer(sqrt(nrow(.data))) + 1,
#' .method = c("silhouette", "gap_stat"))
#'
#' immunr_dbscan(.data, .eps, .dist = TRUE)
#'
#' @param .data Matrix or data frame with features, distance matrix or output from [repOverlapAnalysis] or [geneUsageAnalysis] functions.
#'
#' @param .k The number of clusters to create, defined as `k` to [hcut][factoextra::hcut] or as `centers` to [kmeans].
#'
#' @param .k.max Limits the maximum number of clusters. It is passed as `k.max` to [factoextra::fviz_nbclust] for `immunr_hclust` and `immunr_kmeans`.
#'
#' @param .eps Local radius for expanding clusters, minimal distance between points to expand clusters. Passed as `eps` to [dbscan][fpc::dbscan].
#'
#' @param .method Passed to [factoextra::hcut] or as [factoextra::fviz_nbclust].
#'
#' In case of [factoextra::hcut] the agglomeration method is going to be used (argument `hc_method`).
#'
#' In case of [factoextra::fviz_nbclust] it is the method to be used for estimating the optimal number of clusters (argument `method`).
#'
#' @param .dist If TRUE then ".data" is expected to be a distance matrix. If FALSE then the euclidean distance is computed for the input objects.
#'
#' @return
#' `immunr_hclust` - list with two elements. The first element is an output from [factoextra::hcut].
#' The second element is an output from [factoextra::fviz_nbclust]
#'
#' `immunr_kmeans` - list with three elements. The first element is an output from [kmeans].
#' The second element is an output from [factoextra::fviz_nbclust].
#' The third element is the input dataset `.data`.
#'
#' `immunr_dbscan` - list with two elements. The first element is an output from [fpc::dbscan].
#' The second element is the input dataset `.data`.
#'
#' @examples
#' data(immdata)
#' gu <- geneUsage(immdata$data, .norm = TRUE)
#' immunr_hclust(t(as.matrix(gu[, -1])), .dist = FALSE)
#'
#' gu[is.na(gu)] <- 0
#' immunr_kmeans(t(as.matrix(gu[, -1])))
#' @export immunr_hclust immunr_kmeans immunr_dbscan
immunr_hclust <- function(.data, .k = 2, .k.max = nrow(.data) - 1, .method = "complete", .dist = TRUE) {

  if (!requireNamespace("fpc", quietly = TRUE)) {
    stop("Package 'fpc' is required for this function. Please install it first via install.packages() or devtools::install_github().", call. = FALSE)
  }
  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop("Package 'factoextra' is required for this function. Please install it first via install.packages() or devtools::install_github().", call. = FALSE)
  }

  if (.dist) {
    dist_mat <- as.dist(.data)
  } else {
    dist_mat <- dist(.data)
  }
  res <- list(
    hcut = add_class(factoextra::hcut(dist_mat, k = .k, hc_method = .method), "immunr_hcut"),
    nbclust = add_class(factoextra::fviz_nbclust(.data, factoextra::hcut, k.max = .k.max), "immunr_nbclust")
  )
  add_class(res, "immunr_hclust")
}

immunr_kmeans <- function(.data, .k = 2, .k.max = as.integer(sqrt(nrow(.data))) + 1, .method = c("silhouette", "gap_stat")) {
  res <- list(
    kmeans = add_class(kmeans(.data, .k), "immunr_kmeans"),
    nbclust = add_class(factoextra::fviz_nbclust(.data, kmeans, k.max = .k.max, .method[1]), "immunr_nbclust"),
    data = .data
  )
  add_class(res, "immunr_kmeans")
}

immunr_dbscan <- function(.data, .eps, .dist = TRUE) {
  if (.dist) {
    .data <- as.dist(.data)
    method <- "dist"
  } else {
    method <- "hybrid"
  }
  res <- list(dbscan = fpc::dbscan(.data, eps = .eps, method = method), data = .data)
  add_class(res, "immunr_dbscan")
}
