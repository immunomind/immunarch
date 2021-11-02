#' Main function for exploratory data analysis: compute the distribution of lengths, clones, etc.
#'
#' @importFrom dplyr select
#'
#' @concept explore
#'
#' @aliases repExplore
#'
#' @description The \code{repExplore} function calculates the basic statistics of
#' repertoire: the number of unique immune receptor clonotypes, their relative abundances,
#' and sequence length distribution across the input dataset.
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
#' @param .method A string that specifies the method of analysis. It can be
#' either "volume", "count", "len" or "clones".
#'
#' When .method is set to "volume" the repExplore calculates the number of unique
#' clonotypes in the input data.
#'
#' When .method is set to "count" the repExplore calculates the distribution of
#' clonotype abundances, i.e., how frequent receptors with different abundances are.
#'
#' When .method is set to "len" the repExplore calculates the distribution of
#' CDR3 sequence lengths.
#'
#' When .method is set to "clones" the repExplore returns the number of clones (i.e., cells)
#' per input repertoire.
#'
#' @param .col A string that specifies the column to be processed. Pass "nt" for
#' nucleotide sequence or "aa" for amino acid sequence.
#'
#' @param .coding If \code{TRUE}, then only coding sequences will be analysed.
#'
#' @return
#' If input data is a single immune repertoire, then the function returns a numeric vector
#' with exploratory analysis statistics.
#'
#' Otherwise, it returns a numeric matrix with exploratory analysis statistics for all input repertoires.
#'
#' @seealso \link{vis.immunr_exp_vol}
#'
#' @examples
#' data(immdata)
#'
#' # Calculate statistics and generate a visual output with vis()
#' repExplore(immdata$data, .method = "volume") %>% vis()
#'
#' repExplore(immdata$data, .method = "count") %>% vis()
#'
#' repExplore(immdata$data, .method = "len") %>% vis()
#' @export repExplore
repExplore <- function(.data, .method = c("volume", "count", "len", "clones"), .col = c("nt", "aa"), .coding = TRUE) {
  if (!has_class(.data, "list")) {
    .data <- list(Sample = .data)
  }

  if (.coding) {
    .data <- coding(.data)
  }

  if (.method[1] == "volume") {
    res <- add_class(
      data.frame(
        Sample = names(.data),
        Volume = sapply(.data, function(df) {
          if (has_class(df, "data.table")) {
            df <- df %>% lazy_dt()
          }
          df %>%
            select(IMMCOL$count) %>%
            collect(n = Inf) %>%
            nrow()
        }), stringsAsFactors = FALSE
      ),
      "immunr_exp_vol"
    )
  } else if (.method[1] == "count") {
    res <- lapply(1:length(.data), function(i) {
      count_table <- .data[[i]] %>%
        select(IMMCOL$count) %>%
        collect(n = Inf) %>%
        table()
      data.frame(
        Sample = names(.data)[i],
        Clone.num = as.numeric(names(count_table)),
        Clonotypes = as.vector(count_table), stringsAsFactors = FALSE
      )
    })

    res <- do.call(rbind, res)
    res <- add_class(res, "immunr_exp_count")
  } else if (.method[1] %in% c("len", "lens", "length")) {
    seq_col <- switch(.col[1],
      nt = IMMCOL$cdr3nt,
      aa = IMMCOL$cdr3aa,
      stop("Unknown sequence column: ", .col, ". Please provide either 'nt' or 'aa'")
    )
    res <- lapply(.data, function(df) {
      if (has_class(df, "data.table")) {
        df <- df %>% lazy_dt()
      }
      df %>%
        select(seq_col) %>%
        collect(n = Inf) %>%
        sapply(nchar) %>%
        table()
    })
    res <- lapply(1:length(.data), function(i) {
      data.frame(
        Sample = names(.data)[i],
        Length = as.numeric(names(res[[i]])),
        Count = as.vector(res[[i]]), stringsAsFactors = FALSE
      )
    })
    res <- do.call(rbind, res)
    res <- add_class(res, "immunr_exp_len")
  } else if (.method[1] %in% c("clones", "clone")) {
    res <- add_class(data.frame(
      Sample = names(.data),
      Clones = sapply(.data, function(df) {
        if (has_class(df, "data.table")) {
          df <- df %>% lazy_dt()
        }
        df %>%
          select(IMMCOL$count) %>%
          collect(n = Inf) %>%
          sum(na.rm = TRUE)
      }), stringsAsFactors = FALSE
    ), "immunr_exp_clones")
  } else {
    stop("Unknown method")
    return(NULL)
  }

  res
}

rep.ex <- repExplore
