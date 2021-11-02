if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Barcode", "V1"
  ))
}


#' Select specific clonotypes using barcodes from single-cell metadata
#'
#' @concept single_cell
#'
#' @description Subset the input immune repertoire by barcodes. Pass a vector of
#' barcodes to subset or a vector cluster IDs and corresponding barcodes to
#' get a list of immune repertoires corresponding to cluster IDs.
#' Columns with clonotype counts
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
#' @param .barcodes Either a character vector with barcodes or a named character/factor vector with
#' barcodes as names and cluster IDs a vector elements. The output of Seurat's \code{Idents} function works.
#'
#' @param .force.list Logical. If TRUE then always return a list, even if the result is one data frame.
#'
#' @return An immune repertoire (if ".barcodes" is a barcode vector) or a list of immune repertoires
#' (if ".barcodes" is named vector or an output from Seurat::Idents()). Each element is an immune repertoire
#' with clonotype barcodes corresponding to the input barcodes. The output list's names are cluster names
#' in the ".barcode" argument (Seurat::Idents() case only).
#'
#' @seealso \link{select_clusters}
#'
#' @examples
#' \dontrun{
#' data(immdata)
#' # Create a fake single-cell data
#' df <- immdata$data[[1]]
#' df$Barcode <- "AAAAACCCCC"
#' df$Barcode[51:nrow(df)] <- "GGGGGCCCCC"
#' barcodes <- "AAAAACCCCC"
#' df <- select_barcodes(df, barcodes)
#' nrow(df)
#' }
#' @export
select_barcodes <- function(.data, .barcodes, .force.list = FALSE) {
  if (has_class(.data, "list")) {
    stop("Please apply the function to a repertoire instead of a list of repertoires.
         In order to update both the data and metadata, consider using the select_clusters() function.")
  }
  if (!("Barcode" %in% colnames(.data))) {
    stop('No "Barcode" column in the input data. Did you apply the function to a single-cell repertoire?')
  }

  if (is.null(names(.barcodes))) {
    tmp <- rep("_", length(.barcodes))
    names(tmp) <- .barcodes
    .barcodes <- tmp
  }

  #
  # TODO: make the efficient version without un-nesting the whole data table
  #

  # TODO: is there a way to optimise this without unnecessary conversion?
  convert_to_df <- FALSE
  if (!has_class(.data, "data.table")) {
    df <- as.data.table(.data)
    convert_to_df <- TRUE
  }

  # Un-nest the data table
  df <- df[, c(strsplit(Barcode, IMMCOL_ADD$scsep, useBytes = TRUE, fixed = TRUE), .SD), by = .(Barcode)]
  df[[IMMCOL$count]] <- 1
  df[, Barcode := NULL]
  setnames(df, "", "Barcode")

  # Split by barcodes
  df_list <- split(df, .barcodes[df$Barcode])

  if (length(df_list)) {
    # TODO: is there a way to optimise this without data copying?
    # Group clonotypes
    for (i in 1:length(df_list)) {
      df <- df_list[[i]]
      df <- df[, .(sum(Clones), paste0(Barcode, collapse = IMMCOL_ADD$scsep)), by = setdiff(names(df), c("Clones", "Proportion", "Barcode"))]
      setnames(df, "V1", "Clones")
      setnames(df, "V2", "Barcode")

      df$Proportion <- df$Clones / sum(df$Clones)
      df <- df[order(-Clones)]

      if (convert_to_df) {
        setDF(df)
      }

      df_list[[i]] <- df
    }
  }

  if (.force.list) {
    df_list
  } else if (length(df_list) == 1) {
    df_list[[1]]
  } else {
    df_list
  }
}


#' Split the immune repertoire data to clusters from single-cell barcodes
#'
#' @concept single_cell
#'
#' @description Given the vector of barcodes from Seurat, split the input repertoires
#' to separate subsets following the barcodes' assigned IDs. Useful when you want to
#' split immune repertoires by patients or clusters.
#'
#' @param .data List of two elements "data" and "meta", with "data" being a list of
#' immune repertoires, and "meta" being a metadata table.
#'
#' @param .clusters Factor vector with barcodes as vector names and cluster IDs as vector elements.
#' The output of the Seurat \code{Idents} function works.
#'
#' @param .field A string specifying the name of the field in the input metadata. New immune
#' repertoire subsets will have cluster IDs in this field.
#'
#' @return A list with two elements "data" and "meta" with updated immune repertoire tables and
#' metadata.
#'
#' @seealso \link{select_barcodes}
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
select_clusters <- function(.data, .clusters, .field = "Cluster") {
  if (!has_class(.data, "list")) {
    stop("Please provide a list with both 'data' and 'meta' elements.")
  } else if (!("data" %in% names(.data) && "meta" %in% names(.data))) {
    stop("Your list is missing one of the elements required.
         Please provide a list with both 'data' and 'meta' elements.")
  } else if (!all(sapply(.data$data, function(x) "Barcode" %in% colnames(x)))) {
    stop('No "Barcode" column in the input data. Did you apply the function to a single-cell repertoire?')
  }

  new_data_list <- list()
  new_metadata <- list()

  for (df_i in 1:length(.data$data)) {
    # Select barcodes
    source_name <- names(.data$data)[df_i]
    df_list <- select_barcodes(.data$data[[df_i]], .clusters, TRUE)

    # Update the data and metadata if everything is OK
    if (length(df_list)) {
      for (new_df_i in 1:length(df_list)) {
        cluster_id <- names(df_list)[new_df_i]
        new_name <- paste0(names(.data$data)[df_i], "_", cluster_id)

        new_data_list[[new_name]] <- df_list[[new_df_i]]

        new_metarow <- .data$meta[.data$meta$Sample == source_name, ]
        new_metarow$Sample <- new_name
        new_metarow[[paste0(.field, ".source")]] <- source_name
        new_metarow[[.field]] <- cluster_id

        new_metadata <- c(new_metadata, list(new_metarow))
      }
    }
  }

  list(data = new_data_list, meta = do.call(rbind, new_metadata))
}
