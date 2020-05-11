.datatable.aware <- TRUE

#' Create a repertoire of public clonotypes
#'
#' @concept pubrep
#'
#' @aliases pubRep publicRepertoire
#'
#' @importFrom data.table setcolorder ":=" as.data.table
#' @importFrom dplyr full_join ungroup
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
#' @param .col A string that specifies the column(s) to be processed. Pass one of the
#' following strings, separated by the plus sign: "nt" for nucleotide sequences,
#' "aa" for amino acid sequences, "v" for V gene segments, "j" for J gene segments. E.g.,
#' pass "aa+v" to compute overlaps on CDR3 amino acid sequences paired with V gene segments, i.e.,
#' in this case a unique clonotype is a pair of CDR3 amino acid and V gene segment.
#' @param .quant A string that specifies the column to be processed. Pass "count" to see
#' public clonotype sharing with the number of clones, pass "prop" to see proportions.
#' @param .coding Logical. If TRUE then preprocess the data to filter out non-coding sequences.
#' @param .min.samples Integer. A minimal number of samples a clonotype must have to be included
#' in the public repertoire table.
#' @param .max.samples Integer. A maxminal number of samples a clonotype must have to be included
#' in the public repertoire table. Pass NA (by default) to have the maximal amount of samples.
#' @param .verbose Logical. If TRUE then output the progress.
#'
#' @return
#' Data table with columns for:
#'
#' - Clonotypes (e.g., CDR3 sequence, or two columns for CDR3 sequence and V gene)
#'
#' - Incidence of clonotypes
#'
#' - Per-sample proportions or counts
#'
#' @examples
#' # Subset the data to make the example faster to run
#' immdata$data <- lapply(immdata$data, head, 2000)
#' pr <- pubRep(immdata$data, .verbose=FALSE)
#' vis(pr, "clonotypes", 1, 2)
#' @export pubRep publicRepertoire
pubRep <- function(.data, .col = "aa+v", .quant = c("count", "prop"), .coding = TRUE, .min.samples = 1, .max.samples = NA, .verbose = TRUE) {
  .preprocess <- function(dt) {
    if (has_class(dt, "data.table")) {
      dt <- dt %>% lazy_dt()
    }
    if (.coding) {
      dt <- coding(dt)
    }
    dt <- as.data.table(dt %>% select(.col, .quant) %>% collect(n = Inf))
    dt[, sum(get(.quant)), by = .col]
  }

  assertthat::assert_that(has_class(.data, "list"))

  .col <- sapply(unlist(strsplit(.col, split = "\\+")), switch_type, USE.NAMES = FALSE)
  .quant <- .quant_column_choice(.quant[1])

  if (.verbose) {
    pb <- set_pb(length(.data))
  }

  res <- .preprocess(.data[[1]])
  setnames(res, colnames(res)[ncol(res)], names(.data)[1])

  if (.verbose) {
    add_pb(pb)
  }

  for (i in 2:length(.data)) {
    res <- merge(res, .preprocess(.data[[i]]),
      all = TRUE,
      by = .col, suffixes = c(as.character(i - 1), as.character(i))
    )

    setnames(res, colnames(res)[ncol(res)], names(.data)[i])

    if (.verbose) {
      add_pb(pb)
    }
  }

  if (.verbose) {
    add_pb(pb)
  }

  # res = res %>% ungroup()

  # Add #samples column after .col
  res[["Samples"]] <- rowSums(!is.na(as.matrix(res[, (length(.col) + 1):ncol(res), with = FALSE])))
  if (is.na(.max.samples)) {
    .max.samples <- max(res[["Samples"]])
  }

  res <- res[(Samples >= .min.samples) & (Samples <= .max.samples)]

  if (.verbose) {
    add_pb(pb)
  }

  col_samples <- colnames(res)[(length(.col) + 1):(ncol(res) - 1)]

  res <- setcolorder(res, c(.col, "Samples", col_samples))

  res <- res[order(res$Samples, decreasing = TRUE), ]

  add_class(res, "immunr_public_repertoire")
}


publicRepertoire <- pubRep


#' Get a matrix with public clonotype frequencies
#'
#' @concept pubrep
#'
#' @param .data Public repertoire, an output from \link{pubRep}.
#'
#' @return
#' Matrix with per-sample clonotype counts / proportions only.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 2000)
#' pr <- pubRep(immdata$data, .verbose=FALSE)
#' pr.mat <- public_matrix(pr)
#' dim(pr.mat)
#' head(pr.mat)
#' @export
public_matrix <- function(.data) {
  sample_i <- match("Samples", colnames(.data)) + 1
  max_col <- dim(.data)[2]
  .data %>% dplyr::select(sample_i:max_col) %>% collect(n = Inf) %>% as.matrix()
}


num2bin <- function(number, n_bits) {
  vec <- as.numeric(intToBits(number))
  if (missing(n_bits)) {
    vec
  } else {
    vec[1:n_bits]
  }
}


get_public_repertoire_names <- function(.pr) {
  sample_i <- match("Samples", names(.pr)) + 1
  names(.pr)[sample_i:ncol(.pr)]
}


#' Filter out clonotypes from public repertoires
#'
#' @concept pubrep
#'
#' @aliases pubRepFilter publicRepertoireFilter
#'
#' @description Filter our clonotypes with low incidence in a specific group.
#'
#' @param .pr Public repertoires, an output from \link{pubRep}.
#' @param .meta Metadata file.
#' @param .by Named character vector. Names of the group to filter by.
#' @param .min.samples Integer. Filter out clonotypes with the number of samples below than this number.
#'
#' @return
#' Data frame with filtered clonotypes.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 2000)
#' pr <- pubRep(immdata$data, .verbose=FALSE)
#' pr1 <- pubRepFilter(pr, immdata$meta, .by = c(Status = "MS"))
#' head(pr1)
#' @export pubRepFilter publicRepertoireFilter
pubRepFilter <- function(.pr, .meta, .by, .min.samples = 1) {
  assertthat::assert_that(has_class(.pr, "immunr_public_repertoire"))
  assertthat::assert_that(.min.samples > 0)

  if (!check_group_names(.meta, .by)) {
    return(NA)
  }

  data_groups <- lapply(1:length(.by), function(i) {
    .meta[["Sample"]][.meta[[names(.by)[i]]] == .by[i]]
  })

  samples_of_interest <- data_groups[[1]]

  if (length(.by) > 1) {
    for (i in 2:length(data_groups)) {
      samples_of_interest <- intersect(data_groups[[i]], samples_of_interest)
    }
  }
  samples_of_interest <- intersect(samples_of_interest, get_public_repertoire_names(.pr))

  if (length(samples_of_interest) == 0) {
    message("Warning: no samples found, check group values in the .by argument!")
    return(NA)
  }

  sample_i <- match("Samples", colnames(.pr))
  indices <- c(1:(match("Samples", colnames(.pr))), match(samples_of_interest, colnames(.pr)))
  new.pr <- .pr %>% lazy_dt() %>% dplyr::select(indices) %>% as.data.table()

  new.pr[["Samples"]] <- rowSums(!is.na(as.matrix(new.pr[, (sample_i + 1):ncol(new.pr), with = FALSE])))
  new.pr <- new.pr %>% lazy_dt() %>% dplyr::filter(Samples >= .min.samples) %>% as.data.table()

  new.pr
}

publicRepertoireFilter <- pubRepFilter


#' Apply transformations to public repertoires
#'
#' @concept pubrep
#'
#' @aliases pubRepApply publicRepertoireApply
#'
#' @description Work In Progress
#'
#' @param .pr1 First public repertoire.
#' @param .pr2 Second public repertoire.
#' @param .fun A function to apply to pairs of frequencies of same clonotypes from "pr1" and "pr2".
#' By default - \code{log(X) / log(Y)} where \code{X,Y} - frequencies of the same clonotype,
#' found in both public repertoires.
#'
#' @return
#' Work in progress.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 2000)
#' pr <- pubRep(immdata$data, .verbose=FALSE)
#' pr1 <- pubRepFilter(pr, immdata$meta, .by = c(Status = "MS"))
#' pr2 <- pubRepFilter(pr, immdata$meta, .by = c(Status = "C"))
#' prapp <- pubRepApply(pr1, pr2)
#' head(prapp)
#' @export pubRepApply publicRepertoireApply
pubRepApply <- function(.pr1, .pr2, .fun = function(x) log10(x[1]) / log10(x[2])) {
  col_before_samples <- names(.pr1)[1:(match("Samples", colnames(.pr1)) - 1)]

  # tmp = apply(public_matrix(.pr1), 1, .inner.fun)
  tmp <- rowMeans(public_matrix(.pr1), na.rm = TRUE)
  .pr1[, (match("Samples", colnames(.pr1)) + 1):ncol(.pr1)] <- NULL
  .pr1[["Quant"]] <- tmp

  # tmp = apply(public_matrix(.pr2), 1, .inner.fun)
  tmp <- rowMeans(public_matrix(.pr2), na.rm = TRUE)
  .pr2[, (match("Samples", colnames(.pr2)) + 1):ncol(.pr2)] <- NULL
  .pr2[["Quant"]] <- tmp

  pr.res <- dplyr::inner_join(lazy_dt(.pr1), lazy_dt(.pr2), by = col_before_samples) %>% as.data.table()
  pr.res[["Samples.x"]] <- pr.res[["Samples.x"]] + pr.res[["Samples.y"]]
  pr.res[, Samples.y := NULL]
  names(pr.res)[match("Samples.x", colnames(pr.res))] <- "Samples"

  pr.res[["Result"]] <- apply(pr.res[, c("Quant.x", "Quant.y")], 1, .fun)

  add_class(pr.res, "immunr_public_repertoire_apply")
}

publicRepertoireApply <- pubRepApply


#' Statistics of number of public clonotypes for each possible combinations of repertoires
#'
#' @concept pubrep
#'
#' @importFrom stringr str_count
#' @importFrom dplyr group_by mutate filter
#' @importFrom stats na.omit
#'
#' @param .data Public repertoire, an output from the \link{pubRep} function.
#' @param .by Work in Progress.
#' @param .meta Work in Progress.
#'
#' @return
#' Data frame with incidence statistics per sample.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 2000)
#' pr <- pubRep(immdata$data, .verbose=FALSE)
#' pubRepStatistics(pr) %>% vis()
#' @export pubRepStatistics
pubRepStatistics <- function(.data, .by = NA, .meta = NA) {
  # ToDo: add a support for grouping
  # ToDo: add number of samples per (!) group somehow...
  # ToDo: add a support for V.segments as well

  if (!has_class(.data, "immunr_public_repertoire")) {
    stop("Error: please apply pubRepStatistics() to public repertoires, i.e., output from the pubRep() function.")
  }

  melted_pr <- reshape2::melt(.data, id.vars = colnames(.data)[1:(match("Samples", colnames(.data)))])
  melted_pr <- na.omit(melted_pr)
  melted_pr <- melted_pr %>%
    group_by(CDR3.aa) %>%
    mutate(Group = paste0(variable, collapse = "&")) %>%
    filter(Samples > 1)
  melted_pr <- as_tibble(table(melted_pr$Group), .name_repair = function(x) c("Group", "Count"))
  # melted_pr = bind_cols(melted_pr, Samples = sapply(melted_pr$Group, function (s) stringr::str_count(s, ";") + 1)) %>%
  # arrange(Count) %>%
  # mutate(Group = factor(Group, levels=rev(Group), ordered=TRUE))

  add_class(melted_pr, "immunr_public_statistics")
}


# pubRepAnalysis <- function (.pr, .method, .by = NA, .meta = NA, .min.samples = 1) {
#
# }
