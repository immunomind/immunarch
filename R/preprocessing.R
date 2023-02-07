#' Get the N most abundant clonotypes
#'
#' @concept preprocessing
#'
#' @importFrom dplyr top_n collect
#'
#' @param .data The data to be processed. Can be \link{data.frame},
#' \link{data.table}, or a list of these objects.
#'
#' Every object must have columns in the immunarch compatible format.
#' \link{immunarch_data_format}
#'
#' Competent users may provide advanced data representations:
#' DBI database connections, Apache Spark DataFrame from \link{copy_to} or a list
#' of these objects. They are supported with the same limitations as basic objects.
#'
#' Note: each connection must represent a separate repertoire.
#' @param .n Numeric. Number of the most abundant clonotypes to return.
#'
#' @return
#' Data frame with the \code{.n} most abundant clonotypes only.
#'
#' @examples
#' data(immdata)
#' top(immdata$data)
#' top(immdata$data[[1]])
#' @export top
top <- function(.data, .n = 10) {
  if (has_class(.data, "list")) {
    lapply(.data, top, .n = .n)
  } else {
    # ToDo: fix Clones with IMMDATA$count
    top_n(.data, .n, Clones) %>% collect(n = Inf)
  }
}


#' Filter out coding and non-coding clonotype sequences
#'
#' @concept preprocessing
#'
#' @aliases coding noncoding inframes outofframes
#'
#' @description Filter out clonotypes with non-coding, coding, in-frame or out-of-frame CDR3 sequences:
#'
#' `coding()` - remove all non-coding sequences (i.e., remove all sequences with stop codons and frame shifts);
#'
#' `noncoding()` - remove all coding sequences (i.e., leave sequences with stop codons and frame shifts only);
#'
#' `inframes()` - remove all out-of-frame sequences (i.e., remove all sequences with frame shifts);
#'
#' `outofframes()` - remove all in-frame sequences (i.e., leave sequences with frame shifts only).
#'
#' Note: the function will remove all clonotypes sequences with NAs in the CDR3 amino acid column.
#'
#' @usage
#' coding(.data)
#'
#' noncoding(.data)
#'
#' inframes(.data)
#'
#' outofframes(.data)
#'
#' @param .data The data to be processed. Can be \link{data.frame},
#' \link{data.table}, or a list of these objects.
#'
#' Every object must have columns in the immunarch compatible format.
#' \link{immunarch_data_format}
#'
#' Competent users may provide advanced data representations:
#' DBI database connections, Apache Spark DataFrame from \link{copy_to} or a list
#' of these objects. They are supported with the same limitations as basic objects.
#'
#' Note: each connection must represent a separate repertoire.
#'
#' @return
#' Filtered data frame.
#'
#' @examples
#' data(immdata)
#' immdata_cod <- coding(immdata$data)
#' immdata_cod1 <- coding(immdata$data[[1]])
#' @export coding noncoding inframes outofframes
coding <- function(.data) {
  if (has_class(.data, "list")) {
    lapply(.data, coding)
  } else {
    dt_flag <- FALSE
    if (has_class(.data, "data.table")) {
      dt_flag <- TRUE
      .data %<>% lazy_dt()
    }
    if (has_class(.data, "single_cell")) {
      cdr3_aa_columns <- get_colnames_with(.data, names(IMMCOL_SC["CDR3.aa"]))
    } else {
      cdr3_aa_columns <- IMMCOL$cdr3aa
    }
    d <- collect(.data, n = Inf)
    for (cdr3_aa_column in cdr3_aa_columns) {
      d %<>% filter(!is.na(d[[cdr3_aa_column]]))
      d <- d[grep("[*, ~]", d[[cdr3_aa_column]], invert = TRUE), ]
    }
    if (dt_flag) {
      data.table(d)
    } else {
      d
    }
  }
}

noncoding <- function(.data) {
  if (has_class(.data, "list")) {
    lapply(.data, noncoding)
  } else {
    dt_flag <- FALSE
    if (has_class(.data, "data.table")) {
      dt_flag <- TRUE
      .data <- .data %>% lazy_dt()
    }
    d <- collect(.data, n = Inf)
    d <- filter(d, !is.na(d[[IMMCOL$cdr3aa]]))
    d <- d[grep("[*, ~]", d[[IMMCOL$cdr3aa]], invert = FALSE), ]
    if (dt_flag) {
      data.table(d)
    } else {
      d
    }
  }
}

inframes <- function(.data) {
  if (has_class(.data, "list")) {
    lapply(.data, inframes)
  } else {
    dt_flag <- FALSE
    if (has_class(.data, "data.table")) {
      dt_flag <- TRUE
      .data <- .data %>% lazy_dt()
    }
    d <- collect(.data, n = Inf)
    d <- filter(d, !is.na(d[[IMMCOL$cdr3aa]]))
    # subset(.data, nchar(.data[[IMMCOL$cdr3nt]]) %% 3 == 0)
    d <- d[grep("[~]", d[[IMMCOL$cdr3aa]], invert = TRUE), ]
    if (dt_flag) {
      data.table(d)
    } else {
      d
    }
  }
}

outofframes <- function(.data) {
  if (has_class(.data, "list")) {
    lapply(.data, outofframes)
  } else {
    dt_flag <- FALSE
    if (has_class(.data, "data.table")) {
      dt_flag <- TRUE
      .data <- .data %>% lazy_dt()
    }
    d <- collect(.data, n = Inf)
    d <- filter(d, !is.na(d[[IMMCOL$cdr3aa]]))
    # subset(.data, nchar(.data[[IMMCOL$cdr3nt]]) %% 3 != 0)
    d <- d[grep("[~]", d[[IMMCOL$cdr3aa]], invert = FALSE), ]
    if (dt_flag) {
      data.table(d)
    } else {
      d
    }
  }
}
