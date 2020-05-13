if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Species", "Chain", "Pathology"))
}


#' Load clonotype databases such as VDJDB and McPAS into the R workspace
#'
#' @importFrom readxl read_xlsx
#' @importFrom readr read_csv read_tsv
#'
#' @concept annotation
#'
#' @description The function automatically detects the database format and loads it into R.
#' Additionally, the function provides a general query interface to databases that allows
#' filtering by species, chain types (i.e., locus) and pathology (i.e., antigen species).
#'
#' Currently we support three popular databases:
#'
#' VDJDB - \url{https://github.com/antigenomics/vdjdb-db}
#'
#' McPAS-TCR - \url{http://friedmanlab.weizmann.ac.il/McPAS-TCR/}
#'
#' TBAdb from PIRD - \url{https://db.cngb.org/pird/tbadb/}
#'
#' @param .path Character. A path to the database file, e.g., "/Users/researcher/Downloads/McPAS-TCR.csv".
#'
#' @param .db Character. A database type: either "vdjdb", "vdjdb-search", "mcpas" or "tbadb".
#'
#' "vdjdb" for VDJDB; "vdjdb-search" for search table obtained from the web interface of VDJDB;
#' "mcpas" for McPAS-TCR; "tbadb" for PIRD TBAdb.
#'
#' @param .species Character. A string or a vector of strings specifying which species need to be in the database, e.g., "HomoSapiens".
#' Pass NA (by default) to load all available species.
#'
#' @param .chain Character. A string or a vector of strings specifying which chains need to be in the database, e.g., "TRB".
#' Pass NA (by default) to load all available chains.
#'
#' @param .pathology Character. A string or a vector of strings specifying which disease, virus, bacteria or any condition
#' needs to be in the database, e.g., "CMV".
#' Pass NA (by default) to load all available conditions.
#'
#' @return
#' Data frame with the input database records.
#'
#' @examples
#' # Example file path
#' file_path <- paste0(system.file(package = "immunarch"), "/extdata/db/vdjdb.example.txt")
#'
#' # Load the database with human-only TRB-only receptors for all known antigens
#' db <- dbLoad(file_path, "vdjdb", "HomoSapiens", "TRB")
#' db
#' @export
dbLoad <- function(.path, .db, .species = NA, .chain = NA, .pathology = NA) {
  .db <- tolower(.db)
  if (!(.db %in% c("vdjdb", "vdjdb-search", "mcpas", "mcpas-tcr", "pird", "tbadb"))) {
    stop('Unknown .db argument. Please provide one of the following: "vdjdb", "vdjdb-search", "mcpas" or "tbadb"')
  }

  if (.db == "vdjdb") {
    db_file <- read_tsv(.path, col_types = cols())

    db_file$Species <- db_file$species
    db_file$Chain <- db_file$gene
    db_file$Pathology <- db_file$antigen.species
  } else if (.db == "vdjdb-search") {
    db_file <- read_tsv(.path, col_types = cols())

    db_file$Chain <- db_file$Gene
    db_file$Pathology <- db_file$`Epitope species`
  } else if (.db == "mcpas") {
    db_file <- read_csv(.path, col_types = cols())

    db_file$Chain <- !is.na(db_file$CDR3.beta.aa)
    db_file$Chain <- "TRB"
  } else if (.db == "tbadb") {
    # ToDo: check for conflicting chains, such as TRB and BCR
    sheet_index <- 1
    if (is.na(.chain)[1]) {
      stop("TBAdb requires the .chain argument. Please specify it and try again.")
    }

    chain_col <- paste0(.chain, collapse = "")
    if (grepl("TRG", chain_col) || grepl("TRD", chain_col)) {
      sheet_index <- 2
    } else if (grepl("IGH", chain_col) || grepl("IGL", chain_col) || grepl("IGK", chain_col)) {
      sheet_index <- 3
    }
    db_file <- read_xlsx(.path, sheet_index)

    db_file$Chain <- db_file$Locus
    db_file$Pathology <- db_file$Disease.name
  }

  if (!is.na(.species)) {
    tmp_names <- names(table(db_file$Species))

    # ToDo: fix this, so it check for every tmp_names
    if (!all(.species %in% tmp_names)) {
      stop("Species ", .species, " not found in the input database. Check available species in the database.\n")
    }

    db_file <- db_file %>% filter(Species %in% .species)
  }
  if (!is.na(.chain)) {
    tmp_names <- names(table(db_file$Chain))

    if (!all(.chain %in% tmp_names)) {
      stop("Chain ", .chain, " not found in the input database. Check available chains in the database:\n'", paste0(names(table(db_file$Chain)), collapse = "', '"), "'")
    }

    db_file <- db_file %>% filter(Chain %in% .chain)
  }
  if (!is.na(.pathology)) {
    tmp_names <- names(table(db_file$Pathology))

    if (!all(.pathology %in% tmp_names)) {
      stop("Antigen species ", .pathology, " not found in the input database. Check available antigen species in the database.\n")
    }

    db_file <- db_file %>% filter(Pathology %in% .pathology)
  }

  db_file
}


#' Annotate clonotypes in immune repertoires using clonotype databases such as VDJDB and MCPAS
#'
#' @concept annotation
#'
#' @description Annotate clonotypes using immune receptor databases with known condition-associated receptors.
#' Before using this function, you need to download database files first.
#' For more details see the tutorial \url{https://immunarch.com/articles/web_only/v11_db.html}.
#'
#' @param .data The data to process. It can be a \link{data.frame}, a
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
#' @param .db A data frame or a data table with an immune receptor database. See \link{dbLoad} on how to load databases into R.
#'
#' @param .data.col Character vector. Vector of columns in the input repertoires to use for clonotype search. E.g., `"CDR3.aa"` or `c("CDR3.aa", "V.name")`.
#'
#' @param .db.col Character vector. Vector of columns in the database to use for clonotype search. The order must match the order of ".data.col".
#' E.g., if ".data.col" is `c("CDR3.aa", "V.name")`, then ".db.col" must have the exact order of columns. i.e., the first column must correspond
#' to CDR3 amino acid sequences, and the second column must correspond to V gene segment names.
#'
#' @return
#' Data frame with input sequences and counts or proportions for each of the input repertoire.
#'
#' @examples
#' data(immdata)
#'
#' #' # Example file path
#' file_path <- paste0(system.file(package = "immunarch"), "/extdata/db/vdjdb.example.txt")
#'
#' # Load the database with human-only TRB-only receptors for all known antigens
#' db <- dbLoad(file_path, "vdjdb", "HomoSapiens", "TRB")
#'
#' res <- dbAnnotate(immdata$data, db, "CDR3.aa", "cdr3")
#' res
#' @export dbAnnotate
dbAnnotate <- function(.data, .db, .data.col, .db.col) {
  # Check if the number of columns is equal
  if (length(.data.col) != length(.db.col)) {
    stop("Number of columns in .data.col and .db.col doesn't match! Please provide equal number of columns.")
  }

  if (!has_class(.data, "list")) {
    .data <- list(Sample = .data)
  }

  # Check if columns presented in the input data and in the database
  if (!all(.data.col %in% colnames(.data[[1]]))) {
    stop("Can't find column(s) ", .data.col, " in the input data! Please check if you provided correct column names.")
  }
  if (!all(.db.col %in% colnames(.db))) {
    stop("Can't find column(s) ", .db.col, " in the database! Please check if you provided correct column names.")
  }

  # ToDo: make a more optimal way to column naming
  for (i in 1:length(.db.col)) {
    .db[[.data.col[i]]] <- .db[[.db.col[i]]]
  }

  ann_res <- trackClonotypes(.data, .db %>% select(.data.col), .norm = FALSE)

  ann_res$Samples <- rowSums(ann_res[, 2:ncol(ann_res)] > 0)
  ann_res <- ann_res[ann_res$Samples > 0, ]
  ann_res <- ann_res[order(ann_res$Samples, decreasing = TRUE), ]

  setcolorder(ann_res, c(.data.col, "Samples", names(ann_res)[(length(.data.col) + 1):(ncol(ann_res) - 1)]))

  ann_res
}


# ToDo for dbAnnotate:
# .db as a path to the database file
# Levenshtein search
# more than one database as an input
