#' Function for computing distance for sequences
#'
#' @importFrom stringdist stringdistmatrix
#' @importFrom purrr map
#' @importFrom magrittr %>%

#' @description Computing distance across repertoirs:
#'
#' \code{seqDist} Computing distance between sequences;
#'
#' @usage
#'
#' seqDist(.data)
#'
#' @param .data The data to be processed. Can be \link{data.frame},
##' \link{data.table}, or a list of these objects.
##'
##' Every object must have columns in the immunarch compatible format.
##' \link{immunarch_data_format}
##'
##' Competent users may provide advanced data representations:
##' DBI database connections, Apache Spark DataFrame from \link{copy_to} or a list
##' of these objects. They are supported with the same limitations as basic objects.
##'
##' Note: each connection must represent a separate repertoire.
#'
#' @param .col A string that specifies the column name to be processed.
#'
#' @param .method Passed to \link[stringdist]{stringdistmatrix} or user-defined function.
#'
#' In case of user-defined function, it should take x and y parameters as input and return \link{dist} object.
#'
#' @return
#'
#' Named list of \link{dist} objects for given repertoirs.
#'
#' @examples
#'
#' data(immdata)
#'  # Reducing data to save time on examples
#' immdata$data<-map(immdata$data,~.x %>% head(10))
#'  # Hamming distance computing for each of two first repertoirs
#' seqDist(immdata$data[1:2])
#'
#' # Let's define custom distance function \n
#'   which will count difference in number of characters in sequences.
#'
#' f <- function(x, y) {
#'  res <- matrix(nrow = length(x), ncol = length(y))
#'  for (i in 1:length(x)) {
#'    res[i, ] <- abs(nchar(x[i]) - nchar(y))
#'  }
#'  dimnames(res) <- list(x, y)
#'  return(as.dist(res))
#' }
#'
#' seqDist(immdata$data[1:2],.method = f) # our custom defined distance result
#'
#' @export seqDist

seqDist <- function(.data, .col = "CDR3.nt", .method = "hamming",...) {
  if (!inherits(.data, "list")) { # TODO: here should be general validation methodwith checks .data is list of repertoirs
    stop(paste0(.data, " is not a list of repertoires!"))
  }
  if (!.col %in% colnames(.data[[1]])) {
    stop(paste0("There is no ", .col, " column in ", .data, " data!"))
  } else {
    if (!inherits(.data[[1]][[.col]], "character")) {
      stop("Distance computing are available only for character columns!")
    } else {
      if (inherits(.method, "character")) {
        result <- purrr::map(.data, ~ stringdist::stringdistmatrix(unique(.x[[.col]]), method = .method, useNames = "strings"))
      } else if (inherits(.method, "function")) {
        args<-list(...)
        dist_fun<-function(x){
          args[['x']]<-x
          args[['y']]<-x
          return(do.call(.method,args))}
        result <- purrr::map(.data, ~ .x[[.col]] %>% dist_fun)
      } else {
        stop(".method argument is not a string or a function!")
      }
    }
  }
  return(result)
}


