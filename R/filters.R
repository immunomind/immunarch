#' Main function for data filtering
#'
#' @param .data The data to be processed. Must be the list of 2 elements:
#' data table and metadata table.
#' @param .method Method of filtering. Implemented methods: by.meta.
#' @param .query Filtering query. It's a named list of filters that will be applied
#' to data. Elements of the list are specific for selected method:
#' - by.meta: names are metadata column headers; each value can be string
#' (metadata value to select for the column), vector of strings (to select multiple
#' possible values) or numeric expression in format '>#', '<#', '>=#', '<=#' where
#' '#' is a number (integer or float, can be negative).
#'
#' @examples
#' data(immdata)
#'
#' # Select samples with status "MS"
#' repFilter(immdata, "by.meta", list(Status = "MS"))
#'
#' # Select samples from lanes "A" and "B" with age >= 15
#' repFilter(immdata, "by.meta", list(Lane = c("A", "B"), Age = ">=15"))
#' @export repFilter
repFilter <- function(.data, .method, .query, ...) {
  switch(tolower(.method),
    by.meta = filter_by_meta(.data, .query),
    stop("You entered the wrong method! Please, try again.")
  )
}

filter_by_meta <- function(.data, .query) {
  filtered_meta <- .data$meta
  names <- names(.query)
  for (i in seq_along(names)) {
    name <- names[i]
    column_query <- unlist(.query[name])
    filtered_meta <-
      if (length(column_query) > 1) {
        filter(filtered_meta, get(name) %in% column_query)
      } else if (startsWith(column_query, ">=")) {
        filter(filtered_meta, get(name) >= as.numeric(substring(column_query, 3)))
      } else if (startsWith(column_query, "<=")) {
        filter(filtered_meta, get(name) <= as.numeric(substring(column_query, 3)))
      } else if (startsWith(column_query, ">")) {
        filter(filtered_meta, get(name) > as.numeric(substring(column_query, 2)))
      } else if (startsWith(column_query, "<")) {
        filter(filtered_meta, get(name) < as.numeric(substring(column_query, 2)))
      } else {
        filter(filtered_meta, get(name) == column_query)
      }
  }
  filtered_data <- .data$data[names(.data$data) %in% filtered_meta$Sample]

  result <- list(data = filtered_data, meta = filtered_meta)
  class(result) <- class(.data)
  return(result)
}
