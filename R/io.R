if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "v_call", "d_call", "j_call", "junction", "junction_aa",
    "CDR3.nt", "D.end", "D.name", "D.start", "DJ.ins", "J.name",
    "J.start", "V.end", "V.name", "VD.ins", "VJ.ins", "barcode",
    "chain", "contig_id", "raw_clonotype_id", "raw_consensus_id"
  ))
}



##### Main IO functions #####


#' Load immune repertoire files into the R workspace
#'
#' @concept io
#'
#' @importFrom readr read_delim read_tsv read_csv col_integer col_character col_double col_logical col_guess cols write_lines
#' @importFrom stringr str_split str_detect str_replace_all
#' @importFrom methods as
#' @importFrom dplyr contains first
#' @importFrom utils read.table
#' @importFrom data.table setDF
#'
#' @description The \code{repLoad} function loads repertoire files
#' into R workspace in the immunarch format where you can immediately use them for
#' the analysis. \code{repLoad} automatically detects the right format for
#' your files, so all you need is simply provide the path to your files.
#'
#' See "Details" for more information on supported formats. See "Examples" for
#' diving right into it.
#'
#' @param .path A character string specifying the path to the input data.
#' Input data can be one of the following:
#'
#' - a single repertoire file.
#' In this case \code{repLoad} returns an R \link{data.frame};
#'
#' - a vector of paths to repertoire files.
#' Same as in the case with no metadata file presented in the next section below;
#'
#' - a path to the folder with repertoire files and, if available, metadata file "metadata.txt".
#' If the metadata file if presented, then the \code{repLoad} returns a list with two elements "data" and "meta".
#' "data" is an another list with repertoire R \link{data.frame}s. "meta" is a data frame with the metadata.
#' If the metadata file "metadata.txt" is not presented, then the \code{repLoad} creates a dummy metadata file with
#' sample names and returns a list with two elements "data" and "meta".
#' If input data has multiple chains or cell types stored in the same file
#' (for example, like in 10xGenomics repertoire files), such repertoire files will be splitted to different
#' R data frames with only one type of chain and cell presented. The metadata file will have additional columns specifying
#' cell and chain types for different samples.
#'
#' @param .format A character string specifying what format to use. Do NOT use it. See "Details" for more information on supported formats.
#'
#' Leave NA (which is default) if you want `immunarch` to detect formats automatically.
#'
#' @param .mode Either "single" for single chain data or "paired" for paired chain data.
#'
#' Currently "single" works for every format, and "paired" works only for 10X Genomics data.
#'
#' By default, 10X Genomics data will be loaded as paired chain data, and other files will be loaded as single chain data.
#'
#' @param .coding A logical value. Pass TRUE to get coding-only clonotypes (by defaul). Pass FALSE to get all clonotypes.
#'
#' @details
#' The metadata has to be a tab delimited file with first column named "Sample".
#' It can have any number of additional columns with arbitrary names.
#' The first column should contain base names of files without extensions in
#' your folder. Example:
#' \tabular{llll}{
#'  Sample \tab Sex \tab Age \tab Status\cr
#'  immunoseq_1 \tab M \tab 1 \tab C\cr
#'  immunoseq_2 \tab M \tab 2 \tab C\cr
#'  immunoseq_3 \tab FALSE \tab 3 \tab A
#' }
#'
#' \code{repLoad} has the ".format" argument that sets the format for input repertoire files.
#' Immunarch detects the file format automatically, and the argument is left only for the compatability
#' purposes. It will be soon removed. Do not pass it or your code will stop working!
#'
#' Currently, Immunarch support the following formats:
#'
#' - "immunoseq" - ImmunoSEQ of any version. http://www.adaptivebiotech.com/immunoseq
#'
#' - "mitcr" - MiTCR. https://github.com/milaboratory/mitcr
#'
#' - "mixcr" - MiXCR (the "all" files) of any version. https://github.com/milaboratory/mixcr
#'
#' - "migec" - MiGEC. http://migec.readthedocs.io/en/latest/
#'
#' - "migmap" - For parsing IgBLAST results postprocessed with MigMap. https://github.com/mikessh/migmap
#'
#' - "tcr" - tcR, our previous package. https://imminfo.github.io/tcr/
#'
#' - "vdjtools" - VDJtools of any version. http://vdjtools-doc.readthedocs.io/en/latest/
#'
#' - "imgt" - IMGT HighV-QUEST. http://www.imgt.org/HighV-QUEST/
#'
#' - "airr" - adaptive immune receptor repertoire (AIRR) data format. http://docs.airr-community.org/en/latest/datarep/overview.html
#'
#' - "10x" - 10XGenomics clonotype annotations tables. https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/annotation
#'
#' - "archer" - ArcherDX clonotype tables. https://archerdx.com/
#'
#' @return
#' A list with two named elements:
#'
#'  - "data" is a list of input samples;
#'
#'  - "meta" is a data frame with sample metadata.
#'
#' @seealso \link{immunr_data_format} for immunarch data format; \link{repSave} for file saving;
#' \link{repOverlap}, \link{geneUsage} and \link{repDiversity} for starting with immune repertoires basic statistics.
#'
#' @examples
#' # To load the data from a single file (note that you don't need to specify the data format):
#' file_path <- paste0(system.file(package = "immunarch"), "/extdata/io/Sample1.tsv.gz")
#' immdata <- repLoad(file_path)
#'
#' # Suppose you have a following structure in your folder:
#' # >_ ls
#' # immunoseq1.txt
#' # immunoseq2.txt
#' # immunoseq3.txt
#' # metadata.txt
#'
#' # To load the whole folder with every file in it type:
#' file_path <- paste0(system.file(package = "immunarch"), "/extdata/io/")
#' immdata <- repLoad(file_path)
#' print(names(immdata))
#'
#' # We recommend creating a metadata file named exactly "metadata.txt" in the folder.
#'
#' # In that case, when you load your data you will see:
#' # > immdata <- repLoad("path/to/your/folder/")
#' # > names(immdata)
#' # [1] "data" "meta"
#'
#' # If you do not have "metadata.txt", you will see the same output,
#' # but your metadata will be almost empty:
#' # > immdata <- repLoad("path/to/your/folder/")
#' # > names(immdata)
#' # [1] "data" "meta"
#' @export repLoad
repLoad <- function(.path, .format = NA, .mode = "paired", .coding = TRUE) {
  if (!is.na(.format)) {
    warning("Please don't provide the .format argument,
            immunarch detects the format automatically.
            The .format argument will soon be removed.")
  }

  exclude_extensions <- c("so", "exe", "bam", "fasta", "fai", "fastq", "bed", "rds", "report", "vdjca")

  # Process a repertoire file: detect format and load the data
  # Return: a named list with a repertoire data frame and it's name
  .read_repertoire <- function(.path, .format, .mode, .coding) {
    parse_res <- list()

    # Detect format
    cur_format <- .format
    if (is.na(.format)) {
      cur_format <- .detect_format(.path)
    }

    # Parse the file
    if (is.na(cur_format)) {
      message("unsupported format, skipping")
    }
    else {
      message(cur_format)

      parse_fun <- switch(cur_format,
        mitcr = parse_mitcr,
        migec = parse_migec,
        migmap = parse_migmap,
        mixcr = parse_mixcr,
        tcr = parse_tcr,
        vdjtools = parse_vdjtools,
        imgt = parse_imgt,
        immunoseq = parse_immunoseq,
        airr = parse_airr,
        `10x (consensus)` = parse_10x_consensus,
        `10x (filt.contigs)` = parse_10x_filt_contigs,
        archer = parse_archer,
        immunarch = parse_immunarch,
        catt = parse_catt,
        rtcr = parse_rtcr,
        imseq = parse_imseq,
        vidjil = parse_vidjil,
        NA
      )

      if (suppressWarnings(is.na(parse_fun))) {
        message("unknown format, skipping")
      }
      else {
        parse_res <- parse_fun(.path, .mode)

        if (is.null(parse_res)) {
          message("  [!] Warning: zero clonotypes found, skipping")
          parse_res <- list()
        } else {
          if (.coding) {
            parse_res <- coding(parse_res)
          }

          if (!has_class(parse_res, "list")) {
            parse_res <- list(parse_res)
            names(parse_res) <- .remove.ext(.path)
          }
        }
      }
    }

    parse_res
  }


  # Process a vector with filepaths to repertoire files, metadata and barcodes:
  # just load all repertoire files.
  # Do NOT (!) create a dummy metadata, return en empty data frame instead
  # Return: list with data, metadata and barcodes (if necessary)
  .process_batch <- function(.files, .format, .mode, .coding) {
    parsed_batch <- list()
    metadata <- tibble()

    for (.filepath in .files) {
      message('  -- [', match(.filepath, .files), '/', length(.files), '] Parsing "', .filepath, '" -- ', appendLF = FALSE)
      if (!file.exists(.filepath)) {
        message('Can\'t find\t"', .filepath, '", skipping')
      }
      else {
        # Check for the type: repertoire, metadata or barcodes
        if (stringr::str_detect(.filepath, "metadata")) {
          message("metadata")

          metadata <- readr::read_tsv(.filepath, trim_ws = TRUE, col_types = cols())

          # The most basic check
          if (!("Sample" %in% colnames(metadata))) {
            stop('No "Sample" column found in the metadata file. The "Sample" column with the names of samples without extensions (e.g., ".txt", ".tsv") is required. Please provide it and run the parsing again.')
          }
        }
        else if (stringr::str_detect(.filepath, "barcode")) {
          # TODO: add the barcode processing subroutine to split by samples
        }
        else {
          repertoire <- .read_repertoire(.filepath, .format, .mode, .coding)
          if (length(repertoire) != 0) {
            parsed_batch <- c(parsed_batch, repertoire)
          }
        }
      }
    }

    list(data = parsed_batch, meta = metadata)
  }

  # Recursively find all subdirectories and create
  # a list of batches with repertoire files to process
  .split_to_batches <- function(.paths, .root_folder) {
    # split .paths to files and folders
    folders <- c()
    files <- c()
    for (path in .paths) {
      if (file.info(path)$isdir) {
        folders <- c(folders, path)
      } else if (!any(endsWith(path, exclude_extensions))) {
        files <- c(files, path)
      }
    }

    batches <- list(files)
    names(batches) <- .root_folder

    # process each subdirectory to batches
    for (folder_path in folders) {
      new_batch <- .split_to_batches(dir(folder_path, full.names = TRUE, recursive = FALSE), folder_path)
      batches <- c(batches, new_batch)
    }

    batches
  }


  .check_metadata <- function(.metadata, .rep_names) {
    error_flag <- FALSE

    # Check for zero .metadata
    if (nrow(.metadata) == 0) {
      message("  -- Metadata file not found; creating a dummy metadata...")
      .metadata <- tibble(Sample = .rep_names)
      error_flag <- TRUE
    }

    # Check for repertoire names
    missed_in_folders <- setdiff(.rep_names, .metadata$Sample)
    missed_in_metadata <- setdiff(.metadata$Sample, .rep_names)
    if (length(missed_in_folders) || length(missed_in_metadata)) {
      if (length(missed_in_metadata)) {
        message("  -- Samples found in the metadata, but not in the folder:\n     ", missed_in_metadata)
        message("  Did you correctly specify all the sample names in the metadata file?")

        error_flag <- TRUE
      }
      if (length(missed_in_folders)) {
        message("  -- Samples found in the folder, but not in the metadata:\n     ", missed_in_folders)
        message("  Did you add all the necessary samples to the metadata file with correct names?")
        message("  Creating dummy sample records in the metadata for now...")

        .metadata <- merge(.metadata, as_tibble(data.frame(Sample = missed_in_folders, stringsAsFactors = FALSE)), all = TRUE)

        error_flag <- TRUE
      }
    }

    # Drop incorrect columns
    to_drop <- c()
    for (col_name in names(.metadata)) {
      if (sum(is.na(.metadata[[col_name]])) == nrow(.metadata)) {
        to_drop <- c(to_drop, col_name)
      }
    }

    if (length(to_drop)) {
      .metadata <- .metadata[-match(to_drop, names(.metadata))]
      message("Dropping ", length(to_drop), "column(s) from .metadata.txt. Do you have spaces or tabs after the name of the last column? Remove them to ensure everything works correctly.")
      error_flag <- TRUE
    }

    if (!error_flag) {
      message("  -- Everything is OK!")
    }

    .metadata
  }

  #
  # Step 1: load repertoire files
  #
  message("\n== Step 1/3: loading repertoire files... ==\n")

  # Split .path to files and folders
  batches <- .split_to_batches(.path, "<initial>")

  # Parse repertoires from each batch
  parsed_batches <- list()
  for (batch_i in 1:length(batches)) {
    if (length(batches[[batch_i]])) {
      message('Processing "', names(batches)[batch_i], '" ...')
      parsed_batches[[names(batches)[batch_i]]] <- .process_batch(batches[[batch_i]], .format, .mode, .coding)
    }
  }

  #
  # Step 2: check metadata files
  #
  message("\n== Step 2/3: checking metadata files and merging files... ==\n")

  for (batch_i in 1:length(parsed_batches)) {
    message('Processing "', names(parsed_batches)[batch_i], '" ...')
    parsed_batches[[batch_i]]$meta <- .check_metadata(parsed_batches[[batch_i]]$meta, names(parsed_batches[[batch_i]]$data))
  }

  names(parsed_batches) <- NA
  parsed_batches <- list(data = do.call(c, lapply(parsed_batches, "[[", "data")), meta = Reduce(function(x, y) merge(x, y, all = TRUE), lapply(parsed_batches, "[[", "meta")))
  names(parsed_batches$data) <- substr(names(parsed_batches$data), 4, nchar(names(parsed_batches$data)))

  if (any(duplicated(names(parsed_batches$data)))) {
    stop("Found non-unique sample names: ", paste0(names(parsed_batches$data)[duplicated(names(parsed_batches$data))], collapse = " "), " \n  Please rename your samples so they all have unique names.")
  }

  #
  # Step 3: split by chains and barcodes
  #
  message("\n== Step 3/3: processing paired chain data... ==\n")

  # Check for mixed chains repertoire files and split them if necessary
  rep_names <- names(parsed_batches$data)
  metadata <- parsed_batches$meta

  for (source_name in rep_names) {
    if ("Chain" %in% colnames(metadata)) {
      if (!is.na(metadata[which(metadata$Sample == source_name), "Chain"])) {
        next
      }
    }

    if ("chain" %in% colnames(parsed_batches$data[[source_name]])) {
      #
      # If mode is "single" then split paired chain data to two chain-specific immune repertoires
      #
      if (.mode == "single") {
        df_split <- split(parsed_batches$data[[source_name]], parsed_batches$data[[source_name]]$chain)
        if (length(df_split) > 1) {
          chain_col <- "Chain"
          source_col <- "Source"

          if (!(chain_col %in% colnames(metadata))) {
            message("Splitting repertoires by their chain type. Columns ", chain_col, " and ", source_col, " added to the metadata.")
            message("  -- Column ", chain_col, " specifies chain types: TRA, TRB, etc.")
            message("  -- Column ", source_col, " specifies the initial sample name for a repertoire after splitting by chain types")
          }

          for (i in 1:length(df_split)) {
            sample_name <- paste0(source_name, "_", names(df_split)[i])
            parsed_batches$data[[sample_name]] <- df_split[[i]]
            new_record <- metadata[which(metadata$Sample == source_name), ]
            new_record$Sample <- sample_name
            new_record[[chain_col]] <- names(df_split)[i]
            new_record[[source_col]] <- source_name
            metadata <- merge(metadata, new_record, all = TRUE)
          }

          parsed_batches$data[[source_name]] <- NULL
          metadata <- metadata[-(which(metadata$Sample == source_name)), ]
        }
      }
      #
      # If mode is "paired" then nothing to do here
      #
      else {
        # pass
      }
    }
  }
  message("Done!\n")

  parsed_batches$meta <- metadata
  parsed_batches
}


#' Save immune repertoires to the disk
#'
#' @concept io
#'
#' @importFrom utils packageVersion
#' @importFrom plyr mapvalues
#'
#' @description The \code{repSave} function is deigned to save your data to the disk
#' in desirable format. Currently supports "immunarch" and "vdjtools" file formats.
#'
#' @param .data An R dataframe, a list of R dataframes or a list with \code{data} and
#' \code{meta} where first element is a list of dataframes and the latter is a dataframe
#' with metadata.
#' @param .path A string with the path to the output directory. It should include file
#' name if a single dataframe is provided to .data argument.
#' @param .format A string with desirable format specification. Current options are
#' "immunarch" and "vdjtools".
#' @param .compress A boolean value. Defines whether the output will be compressed or not.
#'
#' @details It is not necessary to create directories beforehand. If the provided directory
#' does not exist it will be created automatically.
#'
#' @return No return value.
#'
#' @examples
#' data(immdata)
#' dirpath <- tempdir()
#' # Save the list of repertoires
#' repSave(immdata, dirpath)
#' # Load it and check if it is the same
#' new_immdata <- repLoad(dirpath)
#' # sum(immdata$data[[1]] != new_immdata$data[[1]], na.rm = TRUE)
#' # sum(immdata$data[[2]] != new_immdata$data[[2]], na.rm = TRUE)
#' # sum(immdata$meta != new_immdata$meta, na.rm = TRUE)
#' @export repSave
repSave <- function(.data, .path, .format = c("immunarch", "vdjtools"),
                    .compress = TRUE) {
  .format <- .format[1]
  if (!.format %in% c("immunarch", "vdjtools")) {
    stop("Unknown format. Please provide either 'immunarch' or 'vdjtools'")
  }

  if (!has_class(.data, "list")) {
    success <- switch(.format,
      immunarch = save_immunarch(.data, .path, .compress),
      vdjtools = save_vdjtools(.data, .path, .compress),
      NA
    )
  } else {
    if ("meta" %in% names(.data)) {
      ifelse(!dir.exists(.path), dir.create(.path), FALSE)
      readr::write_tsv(.data$meta, path = paste0(.path, "/metadata.txt"))
      for (name in names(.data$data)) {
        success <- switch(.format,
          immunarch = save_immunarch(
            .data$data[[name]],
            paste0(.path, "/", name),
            .compress
          ),
          vdjtools = save_vdjtools(
            .data$data[[name]],
            paste0(.path, "/", name),
            .compress
          ),
          NA
        )
      }
    } else {
      ifelse(!dir.exists(.path), dir.create(.path), FALSE)
      for (name in names(.data)) {
        success <- switch(.format,
          immunarch = save_immunarch(
            .data[[name]],
            paste0(.path, "/", name),
            .compress
          ),
          vdjtools = save_vdjtools(
            .data[[name]],
            paste0(.path, "/", name),
            .compress
          ),
          NA
        )
      }
    }
  }
}
