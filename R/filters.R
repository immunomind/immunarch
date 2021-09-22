#' Main function for data filtering
#'
#' @importFrom magrittr "%>%" "%<>%"
#'
#' @param .data The data to be processed. Must be the list of 2 elements:
#' data table and metadata table.
#' @param .method Method of filtering. Implemented methods: by.meta.
#' @param .query Filtering query. It's a named list of filters that will be applied
#' to data. Elements of the list are specific for selected method:
#' - by.meta: names are metadata column headers; each value can be string
#' (metadata value to select for the column), metadata value starting from '!'
#' to filter out this value; vector of strings (to select multiple possible values),
#' vector of strings with first value '!' (to filter out all values from vector)
#' or numeric expression in format '>#', '<#', '>=#', '<=#' where '#' is a number
#' (integer or float, can be negative);
#' - by.repertoire: query parameter is integer: minimal repertoire size for a sample
#' that passes the filter;
#' - by.clonotype: .
#'
#' @examples
#' data(immdata)
#'
#' # Select samples with status "MS"
#' repFilter(immdata, "by.meta", list(Status = "MS"))
#'
#' # Select samples without status "MS"
#' repFilter(immdata, "by.meta", list(Status = "!MS"))
#'
#' # Select samples from lanes "A" and "B" with age >= 15
#' repFilter(immdata, "by.meta", list(Lane = c("A", "B"), Age = ">=15"))
#'
#' # Select samples that not from lanes "A" and "B"
#' repFilter(immdata, "by.meta", list(Lane = c("!", "A", "B")))
#'
#' # Select samples with at least 6000 clonotypes
#' repFilter(immdata, "by.repertoire", 6000)
#' @export repFilter
repFilter <- function(.data, .method = "by.clonotype",
                      .query = list(CDR3.aa = c("!partial", "!out_of_frame"))) {
  if (!is.list(.data)) {
    stop(paste0(
      "Input data is not a list; ",
      "please pass Immunarch dataset object as input."
    ))
  }

  switch(tolower(.method),
    by.meta = filter_by_meta(.data, .query),
    by.repertoire = filter_by_repertoire(.data, .query),
    by.clonotype = filter_by_clonotype(.data, .query),
    stop(paste0(
      "You entered wrong method \"", .method, "\"! Supported methods are: ",
      "\"by.meta\", \"by.repertoire\", \"by.clonotype\"."
    ))
  )
}

filter_by_meta <- function(.data, .query) {
  if (!("meta" %in% names(.data))) {
    stop(paste0(
      "Input data doesn't contain meta; ",
      "please pass Immunarch dataset object with data and metadata as input."
    ))
  }

  filtered_meta <- .data$meta
  for (i in seq_along(names(.query))) {
    name <- names(.query)[i]
    if (!(name %in% names(.data$meta))) {
      stop(paste0("Column \"", name, "\" not found in metadata."))
    }
    column_query <- .query[[name]]

    if (nrow(filtered_meta) > 0) {
      filtered_meta <-
        if (length(column_query) > 1) {
          if (column_query[[1]] == "!") {
            filter(filtered_meta, !get(name) %in% tail(column_query, -1))
          } else {
            filter(filtered_meta, get(name) %in% column_query)
          }
        } else if (startsWith(column_query, "!")) {
          filter(filtered_meta, get(name) != substring(column_query, 2))
        } else if (startsWith(column_query, ">=")) {
          filter(filtered_meta, get(name) >= as_numeric_or_fail(substring(column_query, 3)))
        } else if (startsWith(column_query, "<=")) {
          filter(filtered_meta, get(name) <= as_numeric_or_fail(substring(column_query, 3)))
        } else if (startsWith(column_query, ">")) {
          filter(filtered_meta, get(name) > as_numeric_or_fail(substring(column_query, 2)))
        } else if (startsWith(column_query, "<")) {
          filter(filtered_meta, get(name) < as_numeric_or_fail(substring(column_query, 2)))
        } else {
          filter(filtered_meta, get(name) == column_query)
        }
      if (nrow(filtered_meta) == 0) {
        warning(paste0("Filter by column \"", name, "\" removed all remaining samples!"))
      }
    }
  }
  filtered_data <- .data$data[filtered_meta$Sample]

  return(list(data = filtered_data, meta = filtered_meta))
}

filter_by_repertoire <- function(.data, .min_size) {
  good_repertoires <- names(.data$data)[sapply(.data$data, nrow) >= .min_size]
  filtered_data <- .data$data[good_repertoires]
  filtered_meta <-
    if ("meta" %in% names(.data)) {
      filter(.data$meta, Sample %in% names(filtered_data))
    } else {
      warning("No metadata in input dataset!")
      tibble()
    }

  return(list(data = filtered_data, meta = filtered_meta))
}

filter_by_clonotype <- function(.data, .query) {
  return(list(data = .data$data, meta = .data$meta))
}

include <- function(...) {
  args <- unname(list(...))
  if (length(args) == 0) {
    stop("include() expects at least 1 argument!")
  }
  return(args)
}

exclude <- function(...) {
  args <- unname(list(...))
  if (length(args) == 0) {
    stop("exclude() expects at least 1 argument!")
  }
  return(args)
}

lessthan <- function(value) {
  return(value)
}

morethan <- function(value) {
  return(value)
}

interval <- function(from, to) {
  return(list(from, to))
}
