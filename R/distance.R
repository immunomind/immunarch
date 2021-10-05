#install.packages('stringdist')
#install.packages('proxy')
library(stringdist)
library(proxy)


repDist <- function(.data, .col = "CDR3.nt", .method = "hamming") {
  if (!inherits(.data, "list")) { # TODO: here should be validation that .data is list of repertoirs
    stop(paste0(.data, " is not a list of repertoires!"))
  }
  if (!.col %in% colnames(.data[[1]])) {
    stop(paste0("There is no ", .col, " column in ", .data, " data!"))
  } else {
    if (!inherits(.data[[1]][[.col]], "character")) {
      stop("Distance computing are available only for character columns!")
    } else {
      if (inherits(.method, "character")) {
        # I dont want to check exact method name here due to it can be enchansed with new by stringdist authors
        # It will throw an error if method not to be recognised
        result <- purrr::map(.data, ~ stringdist::stringdistmatrix(.x[[.col]] %>% unique(), method = .method, useNames = "strings"))
      } else if (inherits(.method, "function")) {
        result <- purrr::map(.data, ~ proxy::dist(.x[[.col]] %>% unique(), method = .method))
      } else {
        stop(".method argument is not a string or a function!")
      }
    }
  }
  return(result) # should results saved as dist object or matrix?
}
