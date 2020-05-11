#' Post-analysis of public clonotype statistics: PCA, clustering, etc.
#'
#' @concept overlap
#'
#' @description The \code{\link{repOverlapAnalysis}} function contains advanced data
#' analysis methods. You can use several clustering and dimensionality reduction
#' techniques in order to investigate further the difference between repertoires
#' provided.
#'
#' To cluster a subset of similar data with \code{\link{repOverlapAnalysis}} you can
#' perform hierarchical clustering, k-means or dbscan ('hclust', 'kmeans', 'dbscan'
#' respectively).
#'
#' To reduce dimensions, for example, to select features for subsequent analysis,
#' you can execute the multidimensional scaling or t-sne algorithms ('mds' and 'tsne'
#' respectively).
#'
#' @param .data Any distance matrix between pairs of repertoires. You can also pass your
#' output from \code{\link{repOverlap}}.
#' @param .method A string that defines the type of analysis to perform.
#' @param .scale A function to scale the data before passing it to the MDS algorithm.
#' @param .raw A logical value. Pass TRUE if you want to receive raw output of clustering
#' or dimensionality reduction function of choice. Pass FALSE if you want to receive
#' processed output that can be subjected to visualisation with \code{\link{vis}} function.
#'
#' @param .perp A numerical value, t-SNE parameter, see \code{\link{immunr_tsne}}.
#' @param .theta A numerical value, t-SNE parameter, see \code{\link{immunr_tsne}}.
#'
#' @param .eps A numerical value, DBscan epsylon parameter, see \code{\link{immunr_dbscan}}.
#'
#' @param .k The number of clusters to create, passed as \code{k} to \link[factoextra]{hcut} or as \code{centers} to \link{kmeans}.
#'
#' @return Depends on the last element in the \code{.method} string. See \link{immunr_tsne} for more info.
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' repOverlapAnalysis(ov, "mds+hclust") %>% vis()
#' @export repOverlapAnalysis
repOverlapAnalysis <- function(.data, .method = ("hclust"), .scale = default_scale_fun, .raw = TRUE,
                               .perp = 1, .theta = .1, .eps = .01, .k = 2) {
  if (length(.method) > 2) {
    stop("Can't process the .method parameter with length more than 2.")
  }

  .method <- strsplit(.method, "+", TRUE)[[1]]

  # Preprocessing
  # return_raw = FALSE
  # if (length(.method) == 2) {
  #   return_raw = TRUE
  # }

  if (.method[1] %in% c("mds", "tsne")) {
    res <- switch(.method[1],
      tsne = immunr_tsne(.data, .dist = TRUE, .perp = .perp, .theta = .theta),
      mds = immunr_mds(.data, .scale = .scale, .raw = .raw)
    )
  } else {
    stop("Unknown command. Please provide 'mds' or 'tsne'.")
  }

  if (length(.method) == 1) {
    res
  } else {
    .method <- .method[length(.method)]
    if (.method %in% c("hclust", "kmeans", "dbscan", "permut")) {
      if (has_class(res, "immunr_mds")) {
        res <- res$x
      }

      res <- switch(.method,
        hclust = immunr_hclust(res, .k = .k, .dist = FALSE),
        kmeans = immunr_kmeans(res, .k = .k),
        dbscan = immunr_dbscan(res, .eps = .eps, .dist = FALSE),
        permut = stop(IMMUNR_ERROR_NOT_IMPL)
      )
    } else {
      stop("Unknown command. Please provide 'hclust', 'kmeans', 'dbscan', or 'permut'.")
    }
    res
  }
}
