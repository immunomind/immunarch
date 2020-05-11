#' Immune repertoire spectratyping
#'
#' @concept vis
#'
#' @importFrom dplyr summarise group_by
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
#' @param .quant Select the column with clonal counts to evaluate. Pass "id" to count every clonotype once.
#' Pass "count" to take into the account number of clones per clonotype.
#' @param .col A string that specifies the column(s) to be processed. Pass one of the
#' following strings, separated by the plus sign: "nt" for nucleotide sequences,
#' "aa" for amino acid sequences, "v" for V gene segments, "j" for J gene segments. E.g.,
#' pass "aa+v" for spectratyping on CDR3 amino acid sequences paired with V gene segments, i.e.,
#' in this case a unique clonotype is a pair of CDR3 amino acid and V gene segment.
#' Clonal counts of equal clonotypes will be summed up.
#'
#' @return
#' Data frame with distributions of clonotypes per CDR3 length.
#'
#' @examples
#' # Load the data
#' data(immdata)
#' sp <- spectratype(immdata$data[[1]], .col="aa+v")
#' vis(sp)
#' @export spectratype
spectratype <- function(.data, .quant = c("id", "count"), .col = "nt") {
  assertthat::assert_that(!has_class(.data, "list"))

  .col <- gsub("nuc", "nt", .col) # ToDo: left it here for the backwards compatability, remove it
  if (grepl("nt", .col) && grepl("aa", .col)) {
    stop("Please provide only one sequence column: either 'nt' or 'aa'")
  }
  if (!grepl("nt", .col) && !grepl("aa", .col)) {
    stop("Please provide one sequence column: either 'nt' or 'aa'")
  }
  if (grepl("v", .col) && grepl("j", .col)) {
    stop("Please provide only one gene segment column: either 'v' or 'j'")
  }

  which_cols <- process_col_argument(.col)

  if (!has_class(.data, "data.table")) {
    .data <- as.data.table(.data %>% collect(n = Inf))
  }

  if (.quant[1] == "id") {
    .data[["Spec.count"]] <- 1
    quant_name <- "Spec.count"
  } else {
    quant_name <- IMMCOL$count
  }

  .data <- .data[, sum(get(quant_name)), which_cols]
  setnames(.data, "V1", "Count")

  if (grepl("nt", .col)) {
    .data$Length <- nchar(.data[[IMMCOL$cdr3nt]])
  } else {
    .data$Length <- nchar(.data[[IMMCOL$cdr3aa]])
  }

  if (grepl("v", .col) || grepl("j", .col)) {
    .data <- .data[, sum(Count), c("Length", which_cols[2])]
    setnames(.data, "V1", "Val")
    setnames(.data, if (grepl("v", .col)) {
      "V.name"
    } else {
      "J.name"
    }, "Gene")
    # df = .data %>%
    #   select(Length, Count, Gene = 2) %>%
    #   group_by(Length, Gene) %>%
    #   summarise(Val = sum(Count))
    .data <- as_tibble(.data)
    .data <- .data[order(.data$Val, decreasing = TRUE), ]
    add_class(.data, "immunr_spectr")
  } else {
    .data <- .data[, sum(Count), "Length"]
    setnames(.data, "V1", "Val")
    # df = .data %>%
    #   select(Length, Count) %>%
    #   group_by(Length) %>%
    #   summarise(Val = sum(Count)) %>%
    #   as_tibble()
    .data <- as_tibble(.data)
    add_class(.data, "immunr_spectr_nogene")
  }
}


#' @export
vis.immunr_spectr <- function(.data, .main = "Spectratype", .legend = "Gene segment", .labs = c("CDR3 length", NA), ...) {
  dup <- which(cumsum(!duplicated(.data$Gene)) == 12)[1]
  if (length(dup)) {
    uniq <- unique(.data$Gene[1:(dup - 1)])
    .data$Gene[!(.data$Gene %in% uniq)] <- "Z" # <- dirty hack to avoid factors
  }

  p <- ggplot() +
    geom_bar(aes(x = Length, y = Val, fill = Gene), data = .data, stat = "identity") +
    scale_fill_manual(
      name = .legend, breaks = c(sort(uniq), "Z"),
      labels = c(sort(uniq), "Other"),
      values = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "grey75")
    ) +
    ggtitle(.main)

  if (!is.na(.labs[2])) {
    p <- p + ylab(.labs[2])
  } else {
    p <- p + ylab("Count")
  }
  p <- p + xlab(.labs[1])
  p + theme_linedraw()
}

#' @export
vis.immunr_spectr_nogene <- function(.data, .main = "Spectratype", .legend = "Gene segment", .labs = c("CDR3 length", NA), ...) {
  p <- ggplot() +
    geom_bar(aes(x = Length, y = Val), data = .data, stat = "identity")

  if (!is.na(.labs[2])) {
    p <- p + ylab(.labs[2])
  } else {
    p <- p + ylab("Count")
  }
  p <- p + xlab(.labs[1])
  p + theme_linedraw()
}
