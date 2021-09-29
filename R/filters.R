#' Main function for data filtering
#'
#' @importFrom magrittr "%>%" "%<>%"
#'
#' @param .data The data to be processed. Must be the list of 2 elements:
#' data table and metadata table.
#' @param .method Method of filtering. Implemented methods:
#' by.meta, by.repertoire (by.rep), by.clonotype (by.cl)
#' Default value: 'by.clonotype'.
#' @param .query Filtering query. It's a named list of filters that will be applied
#' to data.
#' Possible values for names in this list are dependent on filter methods:
#' - by.meta: filter by metadata. Names in the named list are metadata column headers.
#' - by.repertoire: filter by number of clonotypes or total number of clones in sample.
#' Possible names in the named list are "n_clonotypes" and "n_clones".
#' - by.clonotype: filter by data in all samples. Names in the named list are
#' data column headers.
#' Elements of the named list for each of the filters are filtering options.
#' Possible values for filtering options:
#' - include("STR1", "STR2", ...): keep only rows with matching values.
#' Available for methods: "by.meta", "by.clonotype".
#' - exclude("STR1", "STR2", ...): remove rows with matching values.
#' Available for methods: "by.meta", "by.clonotype".
#' - lessthan(value): keep rows/samples with numeric values less than specified.
#' Available for methods: "by.meta", "by.repertoire", "by.clonotype".
#' - morethan(value): keep rows/samples with numeric values more than specified.
#' Available for methods: "by.meta", "by.repertoire", "by.clonotype".
#' - interval(from, to): keep rows/samples with numeric values that fits in this interval.
#' from is inclusive, to is exclusive.
#' Available for methods: "by.meta", "by.repertoire", "by.clonotype".
#' Default value: 'list(CDR3.aa = exclude("partial", "out_of_frame"))'.
#' @param .match Matching method for "include" and "exclude" options in query.
#' Possible values:
#' - exact: match only the exact specified string;
#' - startswith: match all strings starting with the specified substring;
#' - substring: match all strings containing the specified substring.
#' Default value: 'exact'.
#'
#' @examples
#' data(immdata)
#'
#' # Select samples with status "MS"
#' repFilter(immdata, "by.meta", list(Status = include("MS")))
#'
#' # Select samples without status "MS"
#' repFilter(immdata, "by.meta", list(Status = exclude("MS")))
#'
#' # Select samples from lanes "A" and "B" with age > 15
#' repFilter(immdata, "by.meta", list(Lane = include("A", "B"), Age = morethan(15)))
#'
#' # Select samples that are not from lanes "A" and "B"
#' repFilter(immdata, "by.meta", list(Lane = exclude("A", "B")))
#'
#' # Select samples with a number of clonotypes from 1000 to 5000
#' repFilter(immdata, "by.repertoire", list(n_clonotypes = interval(1000, 5000)))
#'
#' # Select clonotypes in all samples with alpha chains
#' repFilter(immdata, "by.clonotype",
#'   list(V.name = include("AV"), J.name = include("AJ")),
#'   match = "substring"
#' )
#' @export repFilter
repFilter <- function(.data, .method = "by.clonotype",
                      .query = list(CDR3.aa = exclude("partial", "out_of_frame")),
                      .match = "exact") {
  if (!is.list(.data)) {
    stop(paste0(
      "Input data is not a list; ",
      "please pass Immunarch dataset object as input."
    ))
  }

  if (!(.match %in% c("exact", "startswith", "substring"))) {
    stop(paste0(
      "Unknown matching method \"", .match,
      "\"! Supported matching methods are: ",
      "\"exact\", \"startswith\", \"substring\"."
    ))
  }

  switch(tolower(.method),
    by.meta = filter_by_meta(.data, .query, .match),
    by.repertoire = filter_by_repertoire(.data, .query),
    by.rep = filter_by_repertoire(.data, .query),
    by.clonotype = filter_by_clonotype(.data, .query, .match),
    by.cl = filter_by_clonotype(.data, .query, .match),
    stop(paste0(
      "You entered wrong method \"", .method, "\"! Supported methods are: ",
      "\"by.meta\", \"by.repertoire\", \"by.clonotype\"."
    ))
  )
}

filter_by_meta <- function(.data, .query, .match) {
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
    query_type <- column_query[1]
    query_args <- column_query[-1]
    if (nrow(filtered_meta) > 0) {
      filtered_meta %<>% filter_table(name, query_type, query_args, .match)
      if (nrow(filtered_meta) == 0) {
        warning(paste0("Filter by column \"", name, "\" removed all remaining samples!"))
      }
    }
  }
  filtered_data <- .data$data[filtered_meta$Sample]

  return(list(data = filtered_data, meta = filtered_meta))
}

filter_by_repertoire <- function(.data, .query) {
  filtered_data <- .data$data
  for (i in seq_along(names(.query))) {
    name <- names(.query)[i]
    key_query <- .query[[name]]
    query_type <- key_query[[1]]
    query_args <- key_query[-1]
    if (length(filtered_data) > 0) {
      counts <- switch(name,
        n_clonotypes = lapply(filtered_data, nrow),
        n_clones = lapply(filtered_data, function(sample) {
          sum(sample$Clones)
        }),
        stop(paste0(
          "Wrong filter key for \"by.repertoire\": \"", name,
          "\"! Available keys: \"n_clonotypes\", \"n_clones\"."
        ))
      )

      if (query_type == "lessthan") {
        filtered_data <- filtered_data[
          counts < as_numeric_or_fail(query_args)
        ]
      } else if (query_type == "morethan") {
        filtered_data <- filtered_data[
          counts > as_numeric_or_fail(query_args)
        ]
      } else if (query_type == "interval") {
        from <- as_numeric_or_fail(query_args[[1]])
        to <- as_numeric_or_fail(query_args[[2]])
        filtered_data <- filtered_data[
          counts >= from & counts < to
        ]
      }

      if (nrow(filtered_data) == 0) {
        warning(paste0("Filter by key \"", name, "\" removed all remaining samples!"))
      }
    }
  }
  filtered_meta <-
    if ("meta" %in% names(.data)) {
      filter(.data$meta, Sample %in% names(filtered_data))
    } else {
      warning("No metadata in input dataset!")
      tibble()
    }

  return(list(data = filtered_data, meta = filtered_meta))
}

filter_by_clonotype <- function(.data, .query, .match) {
  filtered_data <- .data$data
  for (i in seq_along(names(.query))) {
    name <- names(.query)[i]
    column_query <- .query[[name]]
    query_type <- column_query[[1]]
    query_args <- column_query[-1]
    for (j in seq_along(names(filtered_data))) {
      sample_name <- names(filtered_data)[j]
      sample <- filtered_data[[sample_name]]
      if (!(name %in% names(sample))) {
        stop(paste0("Column \"", name, "\" not found in sample \"", sample_name, "\"."))
      }
      if (nrow(sample) == 0) {
        warning(paste0("Sample \"", sample_name, "\" was empty, removed!"))
      } else {
        sample %<>% filter_table(name, query_type, query_args, .match)
        if (nrow(sample) == 0) {
          warning(paste0(
            "Filter by column \"", name,
            "\" removed all remaining clonotypes from sample \"",
            sample_name, "\"; sample was removed!"
          ))
        }
      }

      if (nrow(sample) == 0) {
        # removing empty sample
        filtered_data <- filtered_data[names(filtered_data) != sample_name]
      } else {
        # updating sample in filtered_data
        filtered_data[[sample_name]] <- sample
      }
    }
  }
  filtered_meta <-
    if ("meta" %in% names(.data)) {
      filter(.data$meta, Sample %in% names(filtered_data))
    } else {
      warning("No metadata in input dataset!")
      tibble()
    }

  return(list(data = filtered_data, meta = filtered_meta))
}

filter_table <- function(.table, .column_name, .query_type, .query_args, .match) {
  if (.query_type == "include") {
    if (.match == "exact") {
      .table %<>% filter(get(.column_name) %in% .query_args)
    } else if (.match == "startswith") {
      .table <- .table[
        unlist(lapply(lapply(.table[[.column_name]], startsWith, .query_args), any)),
      ]
    } else if (.match == "substring") {
      .table <- .table[
        unique(unlist(lapply(.query_args, grep, .table[[.column_name]], fixed = TRUE))),
      ]
    }
  } else if (.query_type == "exclude") {
    if (.match == "exact") {
      .table %<>% filter(
        !get(.column_name) %in% .query_args
      )
    } else if (.match == "startswith") {
      .table <- .table[
        -unlist(lapply(lapply(.table[[.column_name]], startsWith, .query_args), any)),
      ]
    } else if (.match == "substring") {
      .table <- .table[
        -unique(unlist(lapply(.query_args, grep, .table[[.column_name]], fixed = TRUE))),
      ]
    }
  } else if (.query_type == "lessthan") {
    .table %<>% filter(
      .table, get(.column_name) < as_numeric_or_fail(.query_args)
    )
  } else if (.query_type == "morethan") {
    .table %<>% filter(
      .table, get(.column_name) > as_numeric_or_fail(.query_args)
    )
  } else if (.query_type == "interval") {
    .table %<>% filter(
      .table, get(.column_name) >= as_numeric_or_fail(.query_args[[1]])
    )
    .table %<>% filter(
      .table, get(.column_name) < as_numeric_or_fail(.query_args[[2]])
    )
  }
  return(.table)
}

include <- function(...) {
  args <- unlist(unname(list(...)))
  if (length(args) == 0) {
    stop("include() expects at least 1 argument!")
  }
  return(c("include", args))
}

exclude <- function(...) {
  args <- unlist(unname(list(...)))
  if (length(args) == 0) {
    stop("exclude() expects at least 1 argument!")
  }
  return(c("exclude", args))
}

lessthan <- function(value) {
  return(c("lessthan", unname(value)))
}

morethan <- function(value) {
  return(c("morethan", unname(value)))
}

interval <- function(from, to) {
  return(c("interval", unname(from), unname(to)))
}
