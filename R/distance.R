#' Function for computing distance for sequences
#'
#' @importFrom stringdist stringdistmatrix
#' @importFrom purrr map pmap map2
#' @importFrom magrittr %>% %<>% set_attr
#' @importFrom tidyr unite
#' @importFrom dplyr select_if group_keys group_map group_by group_by_at

#' @description Computing sequential distances between clonotypes from two repertoires:
#'
#' @usage
#'
#' seqDist(.data, .col = 'CDR3.nt', .method = 'hamming',
#'  .group_by = c("V.first", "J.first"), .group_by_seqLength = TRUE, ...)
#'
#' @param .data The data to be processed. Can be \link{data.frame},
#' \link{data.table}, or a list of these objects.
#'
#' Every object must have columns in the immunarch compatible format \link{immunarch_data_format}
#'
#' @param .col A string that specifies the column name to be processed. Default value is 'CDR3.nt'.
#'
#' @param .method Character value or user-defined function.
#'
#' @param .group_by Character vector of column names to group sequence by. Default value is c("V.first", "J.first"). Columns "V.first" and "J.first" containing first genes without allele suffixes are calculated automatically from "V.name" and "J.name" if absent in the data. Pass NA for no grouping options.
#'
#' @param .group_by_seqLength If TRUE  - add grouping by sequence length of .col argument
#'
#' @param ... Extra arguments for user-defined function.
#'
#' Default value is \code{'hamming'} for Hamming distance which counts the number of character substitutions that turns b into a.
#' If a and b have different number of characters the distance is Inf.
#'
#' Other possible values are:
#'
#' \code{'lv'} for Levenshtein distance which counts the number of deletions, insertions and substitutions necessary to turn b into a.
#'
#' \code{'lcs'} for longest common substring is defined as the longest string can be obtained by pairing characters from a and b while keeping the order of characters intact.
#'
#' In case of user-defined function, it should take x and y parameters as input and return \link{dist} object.
#'
#' @return
#'
#' Named list of list with \link{dist} objects for given repertoires for each combination of .group_by variable(s) and/or sequence length of .col.
#'
#' @examples
#'
#' data(immdata)
#' # Reducing data to save time on examples
#' immdata$data <- purrr::map(immdata$data, ~ .x %>% head(10))
#' # Computing hamming distance for the first two repertoires in \code{'immdata'}
#' seqDist(immdata$data[1:2])
#'
#' # Here we define a custom distance function
#' # that will count the difference in number of characters in sequences.
#'
#' f <- function(x, y) {
#'   res <- matrix(nrow = length(x), ncol = length(y))
#'   for (i in 1:length(x)) {
#'     res[i, ] <- abs(nchar(x[i]) - nchar(y))
#'   }
#'   dimnames(res) <- list(x, y)
#'   return(as.dist(res))
#' }
#'
#' seqDist(immdata$data[1:2], .method = f, .group_by_seqLength = FALSE)
#' @export seqDist

seqDist <- function(.data, .col = "CDR3.nt", .method = "hamming", .group_by = c("V.first", "J.first"), .group_by_seqLength = TRUE, ...) {
  .validate_repertoires_data(.data)
  gr_by_is_na <- all(is.na(.group_by))
  # prepare columns with 1st V and J genes if they are used, but not yet calculated
  if ("V.first" %in% .group_by) {
    .data %<>% apply_to_sample_or_list(
      add_column_with_first_gene,
      .original_colname = "V.name",
      .target_colname = "V.first"
    )
  }
  if ("J.first" %in% .group_by) {
    .data %<>% apply_to_sample_or_list(
      add_column_with_first_gene,
      .original_colname = "J.name",
      .target_colname = "J.first"
    )
  }
  # Since seqDist works with any columns of string type, classic .col values are not suported
  if (.col %in% c("aa", "nt", "v", "j", "aa+v")) stop("Please, provide full column name")
  first_sample <- .data[[1]]
  if (!all(.group_by %in% colnames(first_sample)) && !gr_by_is_na) {
    stop("Expected column(s): ", paste0(.group_by, collapse = ", "), "; some of them are missing in data!")
  }
  if (!.col %in% colnames(first_sample)) {
    stop(paste0("There is no ", .col, " column in data!"))
  } else {
    if (!inherits(first_sample[[.col]], "character")) {
      stop("Computing distance is available only for character columns!")
    } else {
      character_dist <- function(x, col = .col, method = .method) {
        stringdist::stringdistmatrix(unique(x[[col]]),
          method = method, useNames = "strings"
        )
      }
      function_dist <- function(x, col = .col, method = .method, args = list(...)) {
        args[["x"]] <- x[[col]]
        args[["y"]] <- x[[col]]
        return(do.call(method, args))
      }
      if (inherits(.method, "character")) {
        dist_fun <- character_dist
      } else if (inherits(.method, "function")) {
        dist_fun <- function_dist
      } else {
        stop(".method argument is not a string nor a function!")
      }
      res_data <- .data
      if (!gr_by_is_na) {
        res_data %<>% map(., ~ .x %>% group_by_at(.group_by))
      }
      if (.group_by_seqLength) {
        res_data %<>% map(., ~ .x %>% group_by(nchar(.x[[.col]]), .add = TRUE))
      }
      result <- map(res_data, ~ .x %>% group_map(~ dist_fun(.)))
      if (!gr_by_is_na) {
        group_by_values <- map(res_data, ~ .x %>%
          group_keys() %>%
          select_if(is.character) %>%
          unite("values", sep = "/"))
        result <- map2(result, group_by_values, ~ map2(.x, .y$values, function(x, y) set_attr(x, "group_values", y)))
      }
    }
  }
  attributes(result)[["col"]] <- .col
  attributes(result)[["group_by"]] <- .group_by
  attributes(result)[["group_by_length"]] <- .group_by_seqLength
  return(result)
}
