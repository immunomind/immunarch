#' Select specific clonotypes using barcodes from single-cell metadata
#'
#' @concept single_cell
#'
#' @description Given the input data and a vector of barcodes,
#' the function returns the input with cells which has
#' at least one of the given barcode. Columns with clonotype counts
#' and proportions are changed accordingly to the filtered barcodes.
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
#' @param .barcodes Character vector. Vector of barcodes to use for
#' filtering.
#'
#' @return An immune repertoire or a list of immune repertoires with clonotypes which have
#' one or more barcodes in the input barcode vector and with clonotype abundances corrected.
#'
#' @examples
#' data(immdata)
#' # Create a fake single-cell data
#' df <- immdata$data[[1]]
#' df$Barcode <- "AAAAACCCCC"
#' df$Barcode[51:nrow(df)] <- "GGGGGCCCCC"
#' barcodes <- "AAAAACCCCC"
#' df <- select_barcodes(df, barcodes)
#' nrow(df)
#' @export
select_barcodes <- function(.data, .barcodes) {
  if (!has_class(.data, "list")) {
    if (!("Barcode" %in% colnames(.data))) {
      stop('No "Barcode" column in the input data. Did you apply the function to a paired chain repertoire?')
    }

    # barcode_count <- sapply(
    #   strsplit(.data[["Barcode"]], ";"),
    #   function(x) sum(.barcodes %in% x)
    # )
    barcode_count <- sapply(
      strsplit(.data[["Barcode"]], ";"),
      function(x) sum(.barcodes %in% x)
    )

    nonzero_bc <- barcode_count > 0
    if (has_class(.data, "data.table")) {
      .data <- .data %>%
        lazy_dt() %>%
        filter(nonzero_bc) %>%
        collect()
    } else {
      .data <- .data %>%
        filter(nonzero_bc)
    }

    barcode_count <- barcode_count[nonzero_bc]
    .data[[IMMCOL$count]] <- barcode_count
    .data[[IMMCOL$prop]] <- barcode_count / sum(barcode_count)
    .data
  } else {
    lapply(.data, select_barcodes, .barcodes = .barcodes)
  }
}


#' Split the immune repertoire data to clusters from single-cell preprocessing
#'
#' @concept single_cell
#'
#' @description Given the vector of barcodes from Seurat, split the input repertoires
#' to separate subsets following the barcodes' assigned cluster or sample labels.
#'
#' @param .clusters Factor vector with barcodes as vector names and cluster IDs as vector elements.
#' The output of the Seurat \code{Idents} function works.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' Idents(pbmc_small)
#' new_cluster_ids <- c("A", "B", "C")
#' new_cluster_ids <- levels(pbmc_small)
#' new_cluster_ids
#' pbmc_small <- RenameIdents(pbmc_small, new_cluster_ids)
#' }
#' @export
select_clusters <- function(.data, .clusters) {
  if (!has_class(.data, "list")) {
    stop("Please provide a list with both 'data' and 'meta' elements.")
  } else if (!("data" %in% names(.data) && "meta" %in% names(.data))) {
    stop("Your list is missing one of the elements required.
         Please provide a list with both 'data' and 'meta' elements.")
  } else {
    if (!all(sapply(.data$data, function (x) "Barcode" %in% colnames(x)))) {
      stop('No "Barcode" column in the input data. Did you apply the function to a paired chain repertoire?')
    }
  }

  #
  # TODO: make the efficient version without un-nesting the whole data table
  #

  # Un-nest the data table
  df <- df[, c(strsplit(Barcode, ";", useBytes = TRUE, fixed=TRUE), .SD), by = .(Barcode)]
  df[[IMMCOL$count]] <- 1
  df[, Barcode:=NULL]
  setnames(df, "", "Barcode")

  # Split by barcodes

  # Group clonotypes
  df[, sum(Clones), by = setdiff(names(df), c("Clones", "Proportion", "Barcode"))][, Clones := V1]
}
