#' Post-analysis of V-gene and J-gene statistics: PCA, clustering, etc.
#'
#' @concept gene_usage
#'
#' @importFrom stats cor
#'
#' @aliases geneUsageAnalysis
#'
#' @description The \code{\link{geneUsageAnalysis}} function deploys several
#' data analysis methods, including PCA, multidimensional scaling,
#' Jensen-Shannon divergence, k-means, hierarchical clustering, DBscan, and different
#' correlation coefficients.
#'
#' @param .data The \code{\link{geneUsageAnalysis}} function runs on the output from
#' \code{\link{geneUsage}}.
#'
#' @param .method A string that defines the type of analysis to perform. Can be "pca",
#' "mds", "js", "kmeans", "hclust", "dbscan" or "cor" if you want to calculate
#' correlation coefficient. In the latter case you have to provide \code{.cor} argument.
#'
#' @param .base A numerical value that defines the logarithm base for Jensen-Shannon
#' divergence.
#' @param .norm.entropy A logical value. Set TRUE to normalise your data if you haven't
#' done it already.
#' @param .cor A string that defines the correlation coefficient for analysis. Can be
#' "pearson", "kendall" or "spearman".
#' @param .do.norm A logical value. If TRUE it forces Laplace smoothing, if NA it checks
#' if smoothing is necessary, if FALSE does nothing.
#' @param .laplace The numeric value, which is used as a pseudocount for Laplace
#' smoothing.
#' @param .verbose A logical value.
#' @param .k The number of clusters to create, passed as \code{k} to \link[factoextra]{hcut} or as \code{centers} to \link{kmeans}.
#' @param .eps A numerical value, DBscan epsylon parameter, see
#' \code{\link{immunr_dbscan}}.
#' @param .perp A numerical value, t-SNE perplexity, see \code{\link{immunr_tsne}}.
#' @param .theta A numerical value, t-SNE theta parameter, see \code{\link{immunr_tsne}}.
#'
#' @return Depends on the last element in the \code{.method} string. See \link{immunr_tsne} for more info.
#'
#' @examples
#' data(immdata)
#' gu <- geneUsage(immdata$data, .norm = TRUE)
#' geneUsageAnalysis(gu, "js+hclust", .verbose = FALSE) %>% vis()
#' @export geneUsageAnalysis
geneUsageAnalysis <- function(.data, .method = c("js+hclust", "pca+kmeans", "anova", "js+pca+kmeans"),
                              .base = 2, .norm.entropy = FALSE,
                              .cor = c("pearson", "kendall", "spearman"),
                              .do.norm = TRUE, .laplace = 1e-12, .verbose = TRUE, .k = 2,
                              .eps = .01, .perp = 1, .theta = .1) {
  .command_preproc <- function(command, .inp) {
    # Yes, there are ifs instead of switch.
    if (command == "js") {
      .inp[is.na(.inp)] <- .laplace
      res <- gene_usage_js(.inp, .base = .base, .norm.entropy = .norm.entropy, .do.norm = .do.norm, .laplace = .laplace, .verbose = .verbose)
    } else if (command == "cor") {
      .inp[is.na(.inp)] <- .laplace
      # do norm here?
      res <- gene_usage_cor(.inp, .cor = .cor[1], .verbose = .verbose)
    } else if (command == "cosine") {
      .inp[is.na(.inp)] <- .laplace
      res <- gene_usage_cosine(.inp, .verbose = .verbose)
    } else if (command == "pca") {
      if (has_class(.inp, "immunr_gene_usage")) {
        .inp <- t(.inp[, -1])
        .inp[is.na(.inp)] <- .000001
        is_dist <- FALSE
      } else {
        .inp[is.na(.inp)] <- 1
        is_dist <- TRUE
      }

      res <- immunr_pca(.inp, .dist = is_dist)
    } else if (command == "mds") {
      if (has_class(.inp, "immunr_gene_usage")) {
        .inp <- t(.inp[, -1])
        .inp[is.na(.inp)] <- .000001
      }

      if (has_class(.inp, "immunr_gene_usage")) {
        is_dist <- FALSE
      } else {
        is_dist <- TRUE
      }

      res <- immunr_mds(.inp)
    } else if (command == "tsne") {
      if (has_class(.inp, "immunr_gene_usage")) {
        .inp <- t(.inp[, -1])
        .inp[is.na(.inp)] <- .000001
      } else {
        .inp[is.na(.inp)] <- 1
      }

      if (has_class(.inp, "immunr_gene_usage")) {
        is_dist <- FALSE
      } else {
        is_dist <- TRUE
      }

      res <- immunr_tsne(.inp, .dist = is_dist, .perp = .perp, .theta = .theta)
    } else {
      stop(paste0("Unknown method for preprocessing: ", command))
    }

    res
  }

  .command_analysis <- function(command, .inp) {
    if (has_class(.inp, "immunr_gu_matrix")) {
      is_dist <- TRUE
    } else {
      is_dist <- FALSE

      if (has_class(.inp, "immunr_pca") || has_class(.inp, "immunr_mds")) {
        .inp <- .inp$x
      }
    }

    if (command == "hclust") {
      res <- immunr_hclust(.inp, .k = .k, .dist = is_dist)
    } else if (command == "kmeans") {
      .inp[is.na(.inp)] <- 0
      res <- immunr_kmeans(.inp, .k = .k)
    } else if (command == "dbscan") {
      res <- immunr_dbscan(.inp, .eps = .eps, .dist = is_dist)
    } else {
      stop(paste0("Unknown method for analysis: ", command))
    }

    res
  }

  commands <- strsplit(.method[1], "+", TRUE)[[1]]

  res <- .command_preproc(commands[1], .data)
  if (length(commands) > 1) {
    if (commands[2] %in% c("js", "cor", "cosine", "pca", "mds", "tsne")) {
      res <- .command_preproc(commands[2], res)
    } else {
      res <- .command_analysis(commands[2], res)
    }

    if (length(commands) > 2) {
      res <- .command_analysis(commands[3], res)
    }
  }

  res
}


gene_usage_js <- function(.data, .base = 2, .do.norm = NA, .laplace = 1e-12, .norm.entropy = FALSE, .verbose = TRUE) {
  res <- apply_symm(.data[, -1], js_div, .base = .base, .do.norm = .do.norm, .laplace = .laplace, .norm.entropy = .norm.entropy, .verbose = .verbose)
  add_class(res, "immunr_gu_matrix")
}

gene_usage_cor <- function(.data, .cor = c("pearson", "kendall", "spearman"), .verbose = TRUE) {
  .cor <- .cor[1]
  res <- apply_symm(.data[, -1], cor, method = .cor, .verbose = .verbose)
  add_class(res, "immunr_gu_matrix")
}

gene_usage_cosine <- function(.data, .verbose = TRUE) {
  res <- apply_symm(.data[, -1], cosine_sim, .verbose = .verbose)
  add_class(res, "immunr_gu_matrix")
}
