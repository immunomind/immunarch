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
#' @importFrom dplyr contains
#' @importFrom utils read.table
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
repLoad <- function(.path, .format = NA, .coding = TRUE) {
  if (!is.na(.format)) {
    warning("Please don't provide the .format argument,
            immunarch detects the format automatically.
            The .format argument will soon be removed.")
  }

  exclude_extensions <- c("so", "exe", "bam", "fasta", "fai", "fastq", "bed", "rds", "report", "vdjca")

  # Process a repertoire file: detect format and load the data
  # Return: a named list with a repertoire data frame and it's name
  .read_repertoire <- function(.path, .format, .coding) {
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
        NA
      )

      if (suppressWarnings(is.na(parse_fun))) {
        message("unknown format, skipping")
      }
      else {
        parse_res <- parse_fun(.path)

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
  .process_batch <- function(.files, .format, .coding) {
    parsed_batch <- list()
    metadata <- tibble()

    for (.filepath in .files) {
      message('  -- Parsing "', .filepath, '" -- ', appendLF = FALSE)
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
          repertoire <- .read_repertoire(.filepath, .format, .coding)
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
      parsed_batches[[names(batches)[batch_i]]] <- .process_batch(batches[[batch_i]], .format, .coding)
    }
  }

  #
  # Step 2: check metadata files
  #
  message("\n== Step 2/3: checking metadata files and merging... ==\n")

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
  message("\n== Step 3/3: splitting data by barcodes and chains... ==\n")

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
#' sum(immdata$data[[1]] != new_immdata$data[[1]], na.rm = TRUE)
#' sum(immdata$data[[2]] != new_immdata$data[[2]], na.rm = TRUE)
#' sum(immdata$meta != new_immdata$meta, na.rm = TRUE)
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





##### IO utility functions #####
.remove.ext <- function(.str) {
  # gsub(pattern = '.*/|[.].*$', replacement = '', x = .str)
  gsub(pattern = ".*/|[.](txt|tsv|csv)$|([.](txt|tsv|csv))?[.](gz|bzip|bzip2|bz2)$", replacement = "", x = .str)
}


.detect_format <- function(.filename) {
  res_format <- NA

  f <- file(.filename, "r")
  l <- readLines(f, 1)
  close(f)

  if (any(str_detect(l, c("MiTCRFullExport", "mitcr")))) {
    res_format <- "mitcr"
  } else if (str_detect(l, "CDR3 amino acid sequence") && str_detect(l, "V segment") && !str_detect(l, "Good events")) {
    res_format <- "mitcr"
  } else if (str_detect(l, "CDR3 amino acid sequence") && str_detect(l, "V segment") && str_detect(l, "Good events")) {
    res_format <- "migec"
  } else if (str_detect(l, "v.end.in.cdr3") && str_detect(l, "cdr3aa")) {
    res_format <- "migmap"
  } else if (str_detect(l, "CDR3.amino.acid.sequence") && str_detect(l, "Umi.count")) {
    res_format <- "tcr"
  } else if (str_detect(tolower(l), "cdr3nt") && str_detect(tolower(l), "vend") && str_detect(tolower(l), "v")) {
    res_format <- "vdjtools"
  } else if (str_detect(tolower(l), "count") && str_detect(tolower(l), "sequence") && str_detect(tolower(l), "d segment")) {
    res_format <- "vdjtools"
  } else if (str_detect(tolower(l), "junction start") && str_detect(tolower(l), "v-d-j-region end") && str_detect(tolower(l), "v-region")) {
    res_format <- "imgt"
  } else if (str_detect(tolower(l), "v_resolved") && str_detect(tolower(l), "amino_acid")) {
    res_format <- "immunoseq"
  } else if (str_detect(tolower(l), "maxresolved")) {
    res_format <- "immunoseq"
  } else if (str_detect(tolower(l), "v_gene") && str_detect(tolower(l), "templates") && str_detect(tolower(l), "amino_acid")) {
    res_format <- "immunoseq"
  } else if (str_detect(tolower(l), "allvalignment") && str_detect(tolower(l), "vhit")) {
    res_format <- "mixcr"
  } else if (str_detect(tolower(l), "bestvhit") && str_detect(tolower(l), "bestjhit")) {
    res_format <- "mixcr"
  } else if (str_detect(tolower(l), "clonal sequence")) {
    res_format <- "mixcr"
  } else if (str_detect(tolower(l), "clonalsequence")) {
    res_format <- "mixcr"
  } else if (str_detect(tolower(l), "targetsequences")) {
    res_format <- "mixcr"
  } else if (str_detect(tolower(l), "junction_aa") && str_detect(tolower(l), "cigar")) {
    res_format <- "airr"
  } else if (str_detect(tolower(l), "raw_clonotype_id") && str_detect(tolower(l), "barcode") && str_detect(tolower(l), "v_gene")) {
    res_format <- "10x (filt.contigs)"
  } else if (str_detect(tolower(l), "clonotype_id") && str_detect(tolower(l), "v_gene")) {
    res_format <- "10x (consensus)"
  } else if (str_detect(tolower(l), "clonotype sequence") && str_detect(tolower(l), "v regions")) {
    res_format <- "archer"
  } else if (str_detect(tolower(l), "exported from immunarch")) {
    res_format <- "immunarch"
  } else if (str_detect(tolower(l), "clones") && str_detect(tolower(l), "v.name") && str_detect(tolower(l), "proportion")) {
    res_format <- "immunarch"
  }

  res_format
}


.make_names <- function(.char) {
  if (is.na(.char[1])) {
    NA
  }
  # else { tolower(make.names(.char)) }
  else {
    tolower(.char)
  }
}


.which_recomb_type <- function(.name) {
  recomb_type <- NA

  i <- 1

  while (is.na(recomb_type) && i < 100) {
    if (any(str_detect(.name[i], c("TCRA", "TRAV", "TCRG", "TRGV", "IGKV", "IGLV")))) {
      recomb_type <- "VJ"
    } else if (any(str_detect(.name[i], c("TCRB", "TRBV", "TCRD", "TRDV", "IGHV")))) {
      recomb_type <- "VDJ"
    }
    i <- i + 1
  }
  if (is.na(recomb_type)) {
    warning("Can't determine the type of V(D)J recombination. No insertions will be presented in the resulting data table.")
  }

  recomb_type
}


.get_coltypes <- function(.filename, .nuc.seq, .aa.seq, .count,
                          .vgenes, .jgenes, .dgenes,
                          .vend, .jstart, .dstart, .dend,
                          .vd.insertions, .dj.insertions, .total.insertions,
                          .skip, .sep, .add = NA) {
  table.colnames <- colnames(readr::read_delim(.filename,
    col_types = cols(),
    delim = .sep,
    quote = "",
    escape_double = FALSE,
    comment = "",
    n_max = 1,
    trim_ws = TRUE,
    skip = .skip
  ))

  swlist <- list(
    col_character(), col_character(),
    col_integer(),
    col_character(), col_character(), col_character(),
    col_integer(), col_integer(), col_integer(), col_integer(),
    col_integer(), col_integer(), col_integer()
  )

  names(swlist) <- tolower(c(
    .nuc.seq, ifelse(is.na(.aa.seq), "NA", .aa.seq),
    .count,
    .vgenes, .jgenes, .dgenes,
    .vend, .jstart, .dstart, .dend,
    .vd.insertions, .dj.insertions, .total.insertions
  ))
  if (!is.na(.add[1])) {
    swlist <- c(swlist, rep(col_guess(), length(.add)))
    names(swlist)[tail(1:length(swlist), length(.add))] <- .add
  }

  swlist <- c(swlist, "_")

  if (is.na(.aa.seq)) {
    swlist <- swlist[-2]
  }

  col.classes <- list(sapply(tolower(table.colnames), function(x) {
    do.call(switch, c(x, swlist))
  }, USE.NAMES = FALSE))[[1]]
  names(col.classes) <- table.colnames

  col.classes
}

.remove.alleles <- function(.data) {
  if (has_class(.data, "list")) {
    lapply(.data, .remove.alleles)
  } else {
    .data[[IMMCOL$v]] <- return_segments(.data[[IMMCOL$v]])
    .data[[IMMCOL$j]] <- return_segments(.data[[IMMCOL$j]])
    .data
  }
}

.postprocess <- function(.data) {
  .data[[IMMCOL$cdr3nt]][.data[[IMMCOL$cdr3nt]] == "NONE"] <- NA
  logic <- is.na(.data[[IMMCOL$cdr3aa]]) & !is.na(.data[[IMMCOL$cdr3nt]])
  if (any(logic)) {
    .data[[IMMCOL$cdr3aa]][logic] <- bunch_translate(.data[[IMMCOL$cdr3nt]][logic])
  }

  logic <- is.na(.data[[IMMCOL$cdr3aa]]) & is.na(.data[[IMMCOL$cdr3nt]])
  if (any(logic)) {
    warn_msg <- c("  [!] Removed ", sum(logic))
    warn_msg <- c(warn_msg, " clonotypes with no nucleotide and amino acid CDR3 sequence.")
    # warn_msg = c(warn_msg, "\n      Please check if your files are the final repertoire files, not the intermediate ones before filtering out bad clonotypes.\n")
    message(warn_msg)
  }
  .data <- .data[!logic, ]

  if (nrow(.data)) {
    # Process 10xGenomics filtered contigs files - count barcodes, merge consensues ids, clonotype ids and contig ids
    if (all(c("contig_id", "barcode") %in% names(.data))) {
      .data <- .data %>%
        group_by(CDR3.nt, V.name, J.name) %>%
        summarise(
          Clones = length(unique(barcode)),
          CDR3.aa = head(CDR3.aa, 1),
          D.name = head(D.name, 1),
          Sequence = head(CDR3.nt, 1),
          V.end = head(V.end, 1),
          J.start = head(J.start, 1),
          D.end = head(D.end, 1),
          D.start = head(D.start, 1),
          VD.ins = head(VD.ins, 1),
          DJ.ins = head(DJ.ins, 1),
          VJ.ins = head(VJ.ins, 1),
          chain = head(chain, 1),
          barcode = paste0(unique(barcode), collapse = IMMCOL_ADD$scsep),
          # raw_clonotype_id = gsub("clonotype", "", paste0(raw_clonotype_id, collapse = IMMCOL_ADD$scsep)),
          # raw_consensus_id = gsub("clonotype|consensus", "", paste0(raw_consensus_id, collapse = IMMCOL_ADD$scsep)),
          contig_id = paste0(contig_id, collapse = IMMCOL_ADD$scsep)
        ) %>%
        ungroup()
      .data[[IMMCOL$prop]] <- .data[[IMMCOL$count]] / sum(.data[[IMMCOL$count]])
      setcolorder(.data, IMMCOL$order)
    }

    for (colname in c(IMMCOL$ve, IMMCOL$ds, IMMCOL$de, IMMCOL$js, IMMCOL$vnj, IMMCOL$vnd, IMMCOL$dnj)) {
      if (colname %in% colnames(.data)) {
        logic <- is.na(.data[[colname]])
        .data[[colname]][logic] <- -1

        logic <- .data[[colname]] < 0
        .data[[colname]][logic] <- NA
      }
    }

    for (col_i in 1:length(IMMCOL$order)) {
      colname <- IMMCOL$order[col_i]
      if (colname %in% colnames(.data)) {
        if (!has_class(.data[[colname]], IMMCOL$type[col_i])) {
          .data[[colname]] <- as(.data[[colname]], IMMCOL$type[col_i])
        }
      }
    }

    logic <- is.na(.data[[IMMCOL$count]])
    if (any(logic)) {
      message("  [!] Warning: found NAs in clonal counts. Setting them to 1's.")
      .data[[IMMCOL$count]][logic] <- 1
    }
    .data <- .data[order(.data[[IMMCOL$count]], decreasing = TRUE), ]
  } else {
    .data <- NULL
  }

  .data
}


##### Parsers #####

parse_repertoire <- function(.filename, .nuc.seq, .aa.seq, .count,
                             .vgenes, .jgenes, .dgenes,
                             .vend, .jstart, .dstart, .dend,
                             .vd.insertions, .dj.insertions, .total.insertions,
                             .skip = 0, .sep = "\t", .add = NA) {
  .nuc.seq <- .make_names(.nuc.seq)
  .aa.seq <- .make_names(.aa.seq)
  .count <- .make_names(.count)
  .vgenes <- .make_names(.vgenes)
  .jgenes <- .make_names(.jgenes)
  .dgenes <- .make_names(.dgenes)
  .vend <- .make_names(.vend)
  .jstart <- .make_names(.jstart)
  .vd.insertions <- .make_names(.vd.insertions)
  .dj.insertions <- .make_names(.dj.insertions)
  .total.insertions <- .make_names(.total.insertions)
  .dstart <- .make_names(.dstart)
  .dend <- .make_names(.dend)
  .add <- .make_names(.add)

  col.classes <- .get_coltypes(.filename, .nuc.seq, .aa.seq, .count,
    .vgenes, .jgenes, .dgenes,
    .vend, .jstart, .dstart, .dend,
    .vd.insertions, .dj.insertions, .total.insertions,
    .skip = .skip, .sep = "\t", .add
  )

  # IO_REFACTOR
  suppressMessages(df <- readr::read_delim(.filename,
    col_names = TRUE,
    col_types = col.classes, delim = .sep,
    quote = "", escape_double = FALSE,
    comment = "", trim_ws = TRUE,
    skip = .skip, na = c("", "NA", ".")
  ))
  # suppressMessages(df <- fread(.filename, skip = .skip, data.table = FALSE, na.strings = c("", "NA", ".")))

  names(df) <- tolower(names(df))
  recomb_type <- .which_recomb_type(df[[.vgenes]])

  table.colnames <- names(col.classes)

  df[[.nuc.seq]] <- toupper(df[[.nuc.seq]])

  if (is.na(.aa.seq)) {
    df$CDR3.amino.acid.sequence <- bunch_translate(df[[.nuc.seq]])
    .aa.seq <- "CDR3.amino.acid.sequence"
  }

  if (is.na(.count)) {
    .count <- "Count"
    df$Count <- 1
  }

  df$Proportion <- df[[.count]] / sum(df[[.count]])
  .prop <- "Proportion"

  ins_ok <- FALSE
  if (is.na(.vd.insertions)) {
    .vd.insertions <- "VD.insertions"
    df$VD.insertions <- NA
  }

  if (!(.vd.insertions %in% table.colnames)) {
    .vd.insertions <- "VD.insertions"
    df$VD.insertions <- NA

    if (!is.na(.vend) && !is.na(.dstart)) {
      if (!is.na(recomb_type) && recomb_type == "VDJ") {
        df$VD.insertions <- df[[.dstart]] - df[[.vend]] - 1
        df$VD.insertions[is.na(df[[.dstart]])] <- NA
        df$VD.insertions[is.na(df[[.vend]])] <- NA

        ins_ok <- TRUE
      }
    }
  }

  if (!ins_ok) {
    df$V.end <- NA
    df$D.start <- NA
    df$D.end <- NA
    .vend <- "V.end"
    .dstart <- "D.start"
    .dend <- "D.end"
  }

  ins_ok <- FALSE
  if (is.na(.dj.insertions)) {
    .dj.insertions <- "DJ.insertions"
    df$DJ.insertions <- NA
  }

  if (!(.dj.insertions %in% table.colnames)) {
    .dj.insertions <- "DJ.insertions"
    df$DJ.insertions <- NA

    if (!is.na(.jstart) && !is.na(.dend)) {
      if (!is.na(recomb_type) && recomb_type == "VDJ") {
        df$DJ.insertions <- df[[.jstart]] - df[[.dend]] - 1
        df$DJ.insertions[is.na(df[[.dend]])] <- NA
        df$DJ.insertions[is.na(df[[.jstart]])] <- NA

        ins_ok <- TRUE
      }
    }
  }
  if (!ins_ok) {
    df$J.start <- NA
    df$D.start <- NA
    df$D.end <- NA
    .jstart <- "J.start"
    .dstart <- "D.start"
    .dend <- "D.end"
  }

  ins_ok <- FALSE
  if (is.na(.total.insertions)) {
    .total.insertions <- "Total.insertions"
    df$Total.insertions <- NA
  }

  if (!(.total.insertions %in% table.colnames)) {
    .total.insertions <- "Total.insertions"
    df$Total.insertions <- -1
    if (!is.na(recomb_type)) {
      if (recomb_type == "VJ") {
        df$Total.insertions <- df[[.jstart]] - df[[.vend]] - 1
        df$Total.insertions[df$Total.insertions < 0] <- 0
      } else if (recomb_type == "VDJ") {
        df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
      }
    } else {
      df$Total.insertions <- NA
    }
  }

  vec_names <- c(
    .count, .prop, .nuc.seq, .aa.seq,
    .vgenes, .dgenes, .jgenes,
    .vend, .dstart, .dend, .jstart,
    .total.insertions, .vd.insertions, .dj.insertions
  )
  if (!is.na(.add[1])) {
    vec_names <- c(vec_names, .add)
  }

  df <- df[, vec_names]

  colnames(df)[1] <- IMMCOL$count
  colnames(df)[2] <- IMMCOL$prop
  colnames(df)[3] <- IMMCOL$cdr3nt
  colnames(df)[4] <- IMMCOL$cdr3aa
  colnames(df)[5] <- IMMCOL$v
  colnames(df)[6] <- IMMCOL$d
  colnames(df)[7] <- IMMCOL$j
  colnames(df)[8] <- IMMCOL$ve
  colnames(df)[9] <- IMMCOL$ds
  colnames(df)[10] <- IMMCOL$de
  colnames(df)[11] <- IMMCOL$js
  colnames(df)[12] <- IMMCOL$vnj
  colnames(df)[13] <- IMMCOL$vnd
  colnames(df)[14] <- IMMCOL$dnj

  .postprocess(df)
}

parse_immunoseq <- function(.filename, .wash.alleles = TRUE) {
  .fix.immunoseq.genes <- function(.col) {
    # fix ","
    .col <- gsub(",", ", ", .col, fixed = TRUE, useBytes = TRUE)
    # fix forward zeros
    .col <- gsub("-([0])([0-9])", "-\\2", .col, useBytes = TRUE)
    .col <- gsub("([VDJ])([0])([0-9])", "\\1\\3", .col, useBytes = TRUE)
    # fix gene names
    .col <- gsub("TCR", "TR", .col, fixed = TRUE, useBytes = TRUE)
    .col
  }

  filename <- .filename
  file_cols <- list()
  file_cols[[IMMCOL$count]] <- "templates"
  file_cols[[IMMCOL$cdr3nt]] <- "rearrangement"
  file_cols[[IMMCOL$cdr3aa]] <- "amino_acid"
  file_cols[[IMMCOL$v]] <- "v_resolved"
  file_cols[[IMMCOL$d]] <- "d_resolved"
  file_cols[[IMMCOL$j]] <- "j_resolved"
  file_cols[[IMMCOL$vnj]] <- "n1_insertions"
  file_cols[[IMMCOL$vnd]] <- "n1_insertions"
  file_cols[[IMMCOL$dnj]] <- "n2_insertions"

  v_index_col_name <- "v_index"
  d_index_col_name <- "d_index"
  j_index_col_name <- "j_index"
  n1_index_col_name <- "n1_index"
  n2_index_col_name <- "n2_index"

  #
  # Check for the version of ImmunoSEQ files
  #
  f <- file(.filename, "r")
  l <- readLines(f, 2)
  close(f)
  if (str_detect(l[[1]], "v_gene") && !str_detect(l[[1]], "v_resolved")) {
    file_cols[[IMMCOL$v]] <- "v_gene"
    file_cols[[IMMCOL$d]] <- "d_gene"
    file_cols[[IMMCOL$j]] <- "j_gene"
  } else if (str_detect(l[[1]], "MaxResolved")) {
    file_cols[[IMMCOL$v]] <- "vMaxResolved"
    file_cols[[IMMCOL$d]] <- "dMaxResolved"
    file_cols[[IMMCOL$j]] <- "jMaxResolved"

    file_cols[[IMMCOL$vnj]] <- "n1insertion"
    file_cols[[IMMCOL$vnd]] <- "n1insertion"
    file_cols[[IMMCOL$dnj]] <- "n2insertion"
  }


  l_split <- strsplit(l, "\t")
  if (str_detect(l[[1]], "templates")) {
    if (str_detect(l[[1]], "templates/reads")) {
      file_cols[[IMMCOL$count]] <- "count (templates/reads)"
      file_cols[[IMMCOL$cdr3nt]] <- "nucleotide"
      file_cols[[IMMCOL$cdr3aa]] <- "aminoAcid"
    } else if (l_split[[2]][match("templates", l_split[[1]])] == "null") {
      file_cols[[IMMCOL$count]] <- "reads"
    }
  }

  if (!str_detect(l[[1]], "v_index")) {
    v_index_col_name <- "vindex"
    d_index_col_name <- "dindex"
    j_index_col_name <- "jindex"
    n1_index_col_name <- "n1index"
    n2_index_col_name <- "n2index"
  }

  for (col_name in names(file_cols)) {
    file_cols[[col_name]] <- .make_names(file_cols[[col_name]])
  }

  file_cols[[IMMCOL$prop]] <- IMMCOL$prop
  file_cols[[IMMCOL$ve]] <- IMMCOL$ve
  file_cols[[IMMCOL$ds]] <- IMMCOL$ds
  file_cols[[IMMCOL$de]] <- IMMCOL$de
  file_cols[[IMMCOL$js]] <- IMMCOL$js
  file_cols[[IMMCOL$seq]] <- IMMCOL$seq

  # IO_REFACTOR
  suppressMessages(df <- readr::read_delim(.filename,
    col_names = TRUE, col_types = cols(),
    delim = "\t", quote = "",
    escape_double = FALSE, comment = "",
    trim_ws = TRUE, skip = 0
  ))
  # suppressMessages(df <- fread(.filename, data.table = FALSE))

  names(df) <- tolower(names(df))

  df[[file_cols[[IMMCOL$prop]]]] <- df[[file_cols[[IMMCOL$count]]]] / sum(df[[file_cols[[IMMCOL$count]]]])

  # Save full nuc sequences and cut them down to CDR3
  df[[IMMCOL$seq]] <- df[[file_cols[[IMMCOL$cdr3nt]]]]

  # TODO: what if df[["v_index]] has "-1" or something like that?
  df[[file_cols[[IMMCOL$cdr3nt]]]] <- stringr::str_sub(df[[IMMCOL$seq]], df[[v_index_col_name]] + 1, nchar(df[[IMMCOL$seq]]))

  df[[file_cols[[IMMCOL$v]]]] <- .fix.immunoseq.genes(df[[file_cols[[IMMCOL$v]]]])
  df[[file_cols[[IMMCOL$d]]]] <- .fix.immunoseq.genes(df[[file_cols[[IMMCOL$d]]]])
  df[[file_cols[[IMMCOL$j]]]] <- .fix.immunoseq.genes(df[[file_cols[[IMMCOL$j]]]])

  recomb_type <- "VDJ"
  if (recomb_type == "VDJ") {
    df[[file_cols[[IMMCOL$ve]]]] <- df[[n1_index_col_name]] - df[[v_index_col_name]]
    df[[file_cols[[IMMCOL$ds]]]] <- df[[d_index_col_name]] - df[[v_index_col_name]]
    df[[file_cols[[IMMCOL$de]]]] <- df[[n2_index_col_name]] - df[[v_index_col_name]]
    df[[file_cols[[IMMCOL$js]]]] <- df[[j_index_col_name]] - df[[v_index_col_name]]
    file_cols[[IMMCOL$vnj]] <- IMMCOL$vnj
    df[[IMMCOL$vnj]] <- -1
  }

  sample_name_vec <- NA
  if ("sample_name" %in% colnames(df)) {
    if (length(unique(df[["sample_name"]])) > 1) {
      sample_name_vec <- df[["sample_name"]]
    }
  }

  df <- df[unlist(file_cols[IMMCOL$order])]
  names(df) <- IMMCOL$order

  if (.wash.alleles) {
    df <- .remove.alleles(df)
    df[[IMMCOL$v]] <- gsub("([VDJ][0-9]*)$", "\\1-1", df[[IMMCOL$v]], useBytes = TRUE)
    df[[IMMCOL$j]] <- gsub("([VDJ][0-9]*)$", "\\1-1", df[[IMMCOL$j]], useBytes = TRUE)
  }

  if (nrow(df) > 0) {
    if (has_class(df[[IMMCOL$vnj]], "character")) {
      df[[IMMCOL$vnj]][df[[IMMCOL$vnj]] == "no data"] <- NA
    }
    if (has_class(df[[IMMCOL$vnd]], "character")) {
      df[[IMMCOL$vnd]][df[[IMMCOL$vnd]] == "no data"] <- NA
    }
    if (has_class(df[[IMMCOL$dnj]], "character")) {
      df[[IMMCOL$dnj]][df[[IMMCOL$dnj]] == "no data"] <- NA
    }

    df[[IMMCOL$vnj]] <- as.integer(df[[IMMCOL$vnj]])
    df[[IMMCOL$vnd]] <- as.integer(df[[IMMCOL$vnd]])
    df[[IMMCOL$dnj]] <- as.integer(df[[IMMCOL$dnj]])
  }

  df[[IMMCOL$v]][df[[IMMCOL$v]] == "unresolved"] <- NA
  df[[IMMCOL$d]][df[[IMMCOL$d]] == "unresolved"] <- NA
  df[[IMMCOL$j]][df[[IMMCOL$j]] == "unresolved"] <- NA

  if (!is.na(sample_name_vec[1])) {
    df <- lapply(split(df, sample_name_vec), .postprocess)
  } else {
    .postprocess(df)
  }
}

parse_mitcr <- function(.filename) {
  .skip <- 0
  f <- file(.filename, "r")
  l <- readLines(f, 1)
  # Check for different levels of the MiTCR output
  if (any(stringr::str_detect(l, c("MiTCRFullExport", "mitcr")))) {
    .skip <- 1
  }
  mitcr_format <- 1
  if (stringr::str_detect(l, "MiTCRFullExport") || .skip == 0) {
    mitcr_format <- 2
  }
  close(f)

  if (mitcr_format == 1) {
    filename <- .filename
    .count <- "count"
    nuc.seq <- "cdr3nt"
    aa.seq <- "cdr3aa"
    vgenes <- "v"
    jgenes <- "j"
    dgenes <- "d"
    vend <- "VEnd"
    jstart <- "JStart"
    dstart <- "DStart"
    dend <- "DEnd"
    vd.insertions <- NA
    dj.insertions <- NA
    total.insertions <- NA
    .sep <- "\t"
  } else {
    # Check if there are barcodes
    f <- file(.filename, "r")
    l <- readLines(f, 1 + .skip)[.skip + 1]
    barcodes <- NA
    .count <- "Read count"
    if ("NNNs" %in% strsplit(l, "\t", TRUE)[[1]]) {
      .count <- "NNNs"
    }
    close(f)

    filename <- .filename
    nuc.seq <- "CDR3 nucleotide sequence"
    aa.seq <- "CDR3 amino acid sequence"
    vgenes <- "V segments"
    jgenes <- "J segments"
    dgenes <- "D segments"
    vend <- "Last V nucleotide position"
    jstart <- "First J nucleotide position"
    dstart <- "First D nucleotide position"
    dend <- "Last D nucleotide position"
    vd.insertions <- "VD insertions"
    dj.insertions <- "DJ insertions"
    total.insertions <- "Total insertions"
    .sep <- "\t"
  }

  parse_repertoire(
    .filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
    .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
    .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
    .total.insertions = total.insertions, .skip = .skip, .sep = .sep
  )
}

parse_mixcr <- function(.filename) {
  fix.alleles <- function(.data) {
    .data[[IMMCOL$v]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$v]])
    .data[[IMMCOL$d]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$d]])
    .data[[IMMCOL$j]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$j]])
    .data
  }

  .filename <- .filename
  .count <- "clonecount"
  .sep <- "\t"
  .vend <- "allvalignments"
  .jstart <- "alljalignments"
  .dalignments <- "alldalignments"
  .vd.insertions <- "VD.insertions"
  .dj.insertions <- "DJ.insertions"
  .total.insertions <- "Total.insertions"

  table.colnames <- tolower(make.names(read.table(.filename, sep = .sep, skip = 0, nrows = 1, stringsAsFactors = FALSE, strip.white = TRUE, comment.char = "", quote = "")[1, ]))
  table.colnames <- gsub(".", "", table.colnames, fixed = TRUE)

  # Columns of different MiXCR formats
  # Clone count - Clonal sequence(s) - N. Seq. CDR3
  # cloneCount - clonalSequence - nSeqCDR3
  # cloneCount - targetSequences - nSeqImputedCDR3
  # cloneCount - targetSequences - nSeqCDR3
  if ("targetsequences" %in% table.colnames) {
    if ("nseqimputedcdr3" %in% table.colnames) {
      .nuc.seq <- "nseqimputedcdr3"
    } else {
      .nuc.seq <- "nseqcdr3"
    }

    .big.seq <- "targetsequences"
  } else {
    .nuc.seq <- "nseqcdr3"

    if ("clonalsequences" %in% table.colnames) {
      .big.seq <- "clonalsequences"
    } else if ("clonalsequence" %in% table.colnames) {
      .big.seq <- "clonalsequence"
    } else {
      .big.seq <- NA
    }
  }

  if (!("allvalignments" %in% table.colnames)) {
    if ("allvalignment" %in% table.colnames) {
      .vend <- "allvalignment"
    } else {
      .vend <- NA
    }
  }
  if (!("alldalignments" %in% table.colnames)) {
    if ("alldalignment" %in% table.colnames) {
      .dalignments <- "alldalignment"
    } else {
      .dalignments <- NA
    }
  }
  if (!("alljalignments" %in% table.colnames)) {
    if ("alljalignment" %in% table.colnames) {
      .jstart <- "alljalignment"
    } else {
      .jstart <- NA
    }
  }

  if ("bestvhit" %in% table.colnames) {
    .vgenes <- "bestvhit"
  } else if ("allvhits" %in% table.colnames) {
    .vgenes <- "allvhits"
  } else if ("vhits" %in% table.colnames) {
    .vgenes <- "vhits"
  } else if ("allvhitswithscore" %in% table.colnames) {
    .vgenes <- "allvhitswithscore"
  } else {
    message("Error: can't find a column with V genes")
  }

  if ("bestjhit" %in% table.colnames) {
    .jgenes <- "bestjhit"
  } else if ("alljhits" %in% table.colnames) {
    .jgenes <- "alljhits"
  } else if ("jhits" %in% table.colnames) {
    .jgenes <- "jhits"
  } else if ("alljhitswithscore" %in% table.colnames) {
    .jgenes <- "alljhitswithscore"
  } else {
    message("Error: can't find a column with J genes")
  }

  if ("bestdhit" %in% table.colnames) {
    .dgenes <- "bestdhit"
  } else if ("alldhits" %in% table.colnames) {
    .dgenes <- "alldhits"
  } else if ("dhits" %in% table.colnames) {
    .dgenes <- "dhits"
  } else if ("alldhitswithscore" %in% table.colnames) {
    .dgenes <- "alldhitswithscore"
  } else {
    message("Error: can't find a column with D genes")
  }


  # IO_REFACTOR
  df <- read_delim(
    file = .filename, col_types = cols(),
    delim = .sep, skip = 0, comment = "",
    quote = "", escape_double = FALSE, trim_ws = TRUE
  )
  # df <- fread(.filename, data.table = FALSE)

  #
  # return NULL if there is no clonotypes in the data frame
  #
  if (nrow(df) == 0) {
    return(NULL)
  }

  names(df) <- make.names(names(df))
  names(df) <- tolower(gsub(".", "", names(df), fixed = TRUE))
  names(df) <- str_replace_all(names(df), " ", "")

  # check for VJ or VDJ recombination
  # VJ / VDJ / Undeterm
  recomb_type <- "Undeterm"
  if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRA", "TRAV", "TRGV", "IGKV", "IGLV"))) {
    recomb_type <- "VJ"
  } else if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRB", "TRBV", "TRDV", "IGHV"))) {
    recomb_type <- "VDJ"
  }

  if (!is.na(.vend) && !is.na(.jstart)) {
    .vd.insertions <- "VD.insertions"
    df$VD.insertions <- -1
    if (recomb_type == "VJ") {
      df$VD.insertions <- -1
    } else if (recomb_type == "VDJ") {
      logic <- sapply(strsplit(df[[.dalignments]], "|", TRUE, FALSE, TRUE), length) >= 4 &
        sapply(strsplit(df[[.vend]], "|", TRUE, FALSE, TRUE), length) >= 5
      df$VD.insertions[logic] <-
        as.numeric(sapply(strsplit(df[[.dalignments]][logic], "|", TRUE, FALSE, TRUE), "[[", 4)) -
        as.numeric(sapply(strsplit(df[[.vend]][logic], "|", TRUE, FALSE, TRUE), "[[", 5)) - 1
    }

    .dj.insertions <- "DJ.insertions"
    df$DJ.insertions <- -1
    if (recomb_type == "VJ") {
      df$DJ.insertions <- -1
    } else if (recomb_type == "VDJ") {
      logic <- sapply(strsplit(df[[.jstart]], "|", TRUE, FALSE, TRUE), length) >= 4 &
        sapply(strsplit(df[[.dalignments]], "|", TRUE, FALSE, TRUE), length) >= 5
      df$DJ.insertions[logic] <-
        as.numeric(sapply(strsplit(df[[.jstart]][logic], "|", TRUE, FALSE, TRUE), "[[", 4)) -
        as.numeric(sapply(strsplit(df[[.dalignments]][logic], "|", TRUE, FALSE, TRUE), "[[", 5)) - 1
    }

    # VJ.insertions
    logic <- (sapply(strsplit(df[[.vend]], "|", TRUE, FALSE, TRUE), length) > 4) & (sapply(strsplit(df[[.jstart]], "|", TRUE, FALSE, TRUE), length) >= 4)
    .total.insertions <- "Total.insertions"
    if (recomb_type == "VJ") {
      df$Total.insertions <- NA
      if (length(which(logic)) > 0) {
        df$Total.insertions[logic] <-
          as.numeric(sapply(strsplit(df[[.jstart]][logic], "|", TRUE, FALSE, TRUE), "[[", 4)) - as.numeric(sapply(strsplit(df[[.vend]][logic], "|", TRUE, FALSE, TRUE), "[[", 5)) - 1
      }
    } else if (recomb_type == "VDJ") {
      df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
    } else {
      df$Total.insertions <- NA
    }
    df$Total.insertions[df$Total.insertions < 0] <- -1

    df$V.end <- -1
    df$J.start <- -1
    df[[.vend]] <- gsub(";", "", df[[.vend]], fixed = TRUE)
    logic <- sapply(strsplit(df[[.vend]], "|", TRUE, FALSE, TRUE), length) >= 5
    df$V.end[logic] <- sapply(strsplit(df[[.vend]][logic], "|", TRUE, FALSE, TRUE), "[[", 5)
    logic <- sapply(strsplit(df[[.jstart]], "|", TRUE, FALSE, TRUE), length) >= 4
    df$J.start[logic] <- sapply(strsplit(df[[.jstart]][logic], "|", TRUE, FALSE, TRUE), "[[", 4)
  } else {
    df$V.end <- -1
    df$J.start <- -1
    df$Total.insertions <- -1
    df$VD.insertions <- -1
    df$DJ.insertions <- -1

    .dj.insertions <- "DJ.insertions"
    .vd.insertions <- "VD.insertions"
  }

  .vend <- "V.end"
  .jstart <- "J.start"

  if (!is.na(.dalignments)) {
    logic <- sapply(str_split(df[[.dalignments]], "|"), length) >= 5
    df$D5.end <- -1
    df$D3.end <- -1
    df$D5.end[logic] <- sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 4)
    df$D3.end[logic] <- sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 5)
    .dalignments <- c("D5.end", "D3.end")
  } else {
    df$D5.end <- -1
    df$D3.end <- -1
  }

  .dalignments <- c("D5.end", "D3.end")

  if (!(.count %in% table.colnames)) {
    warn_msg <- c("  [!] Warning: can't find a column with clonal counts. Setting all clonal counts to 1.")
    warn_msg <- c(warn_msg, "\n      Did you apply repLoad to MiXCR file *_alignments.txt?")
    warn_msg <- c(warn_msg, " If so please consider moving all *.clonotypes.*.txt MiXCR files to")
    warn_msg <- c(warn_msg, " a separate folder and apply repLoad to the folder.")
    warn_msg <- c(warn_msg, "\n      Note: The *_alignments.txt file IS NOT a repertoire file suitable for any analysis.")
    message(warn_msg)

    df[[.count]] <- 1
  }
  .freq <- "Proportion"
  df$Proportion <- df[[.count]] / sum(df[[.count]], na.rm = TRUE)

  .aa.seq <- IMMCOL$cdr3aa
  df[[.aa.seq]] <- bunch_translate(df[[.nuc.seq]])

  if (is.na(.big.seq)) {
    .big.seq <- "BigSeq"
    df$BigSeq <- df[[.nuc.seq]]
  }

  df <- df[, make.names(c(
    .count, .freq,
    .nuc.seq, .aa.seq,
    .vgenes, .dgenes, .jgenes,
    .vend, .dalignments, .jstart,
    .total.insertions, .vd.insertions, .dj.insertions, .big.seq
  ))]

  colnames(df) <- IMMCOL$order

  df[[IMMCOL$v]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.,]*[[:digit:]]*[)])", "", df[[IMMCOL$v]])
  df[[IMMCOL$v]] <- gsub(",", ", ", df[[IMMCOL$v]])
  df[[IMMCOL$v]] <- str_replace_all(df[[IMMCOL$v]], '"', "")
  df[[IMMCOL$v]] <- sapply(
    strsplit(df[[IMMCOL$v]], ", ", useBytes = TRUE),
    function(x) paste0(sort(unique(x)), collapse = ", ")
  )

  df[[IMMCOL$d]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.,]*[[:digit:]]*[)])", "", df[[IMMCOL$d]])
  df[[IMMCOL$d]] <- gsub(",", ", ", df[[IMMCOL$d]])
  df[[IMMCOL$d]] <- str_replace_all(df[[IMMCOL$d]], '"', "")
  df[[IMMCOL$d]] <- sapply(
    strsplit(df[[IMMCOL$d]], ", ", useBytes = TRUE),
    function(x) paste0(sort(unique(x)), collapse = ", ")
  )

  df[[IMMCOL$j]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.,]*[[:digit:]]*[)])", "", df[[IMMCOL$j]])
  df[[IMMCOL$j]] <- gsub(",", ", ", df[[IMMCOL$j]])
  df[[IMMCOL$j]] <- str_replace_all(df[[IMMCOL$j]], '"', "")
  df[[IMMCOL$j]] <- sapply(
    strsplit(df[[IMMCOL$j]], ", ", useBytes = TRUE),
    function(x) paste0(sort(unique(x)), collapse = ", ")
  )

  .postprocess(fix.alleles(df))
}

parse_migec <- function(.filename) {
  filename <- .filename
  nuc.seq <- "CDR3 nucleotide sequence"
  aa.seq <- "CDR3 amino acid sequence"
  .count <- "Good events"
  vgenes <- "V segments"
  jgenes <- "J segments"
  dgenes <- "D segments"
  vend <- "Last V nucleotide position"
  jstart <- "First J nucleotide position"
  dstart <- "First D nucleotide position"
  dend <- "Last D nucleotide position"
  vd.insertions <- "VD insertions"
  dj.insertions <- "DJ insertions"
  total.insertions <- "Total insertions"
  .skip <- 0
  .sep <- "\t"

  parse_repertoire(
    .filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count,
    .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
    .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
    .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
    .total.insertions = total.insertions, .skip = .skip, .sep = .sep
  )
}

parse_migmap <- function(.filename) {
  filename <- .filename
  nuc.seq <- "cdr3nt"
  aa.seq <- "cdr3aa"
  .count <- "count"
  vgenes <- "v"
  jgenes <- "j"
  dgenes <- "d"
  vend <- "v.end.in.cdr3"
  jstart <- "j.start.in.cdr3"
  dstart <- "d.start.in.cdr3"
  dend <- "d.end.in.cdr3"
  vd.insertions <- NA
  dj.insertions <- NA
  total.insertions <- NA
  .skip <- 0
  .sep <- "\t"

  parse_repertoire(
    .filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
    .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
    .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
    .total.insertions = total.insertions, .skip = .skip, .sep = .sep
  )
}

parse_tcr <- function(.filename) {
  f <- file(.filename, "r")
  l <- readLines(f, 2)[2]
  close(f)

  nuc.seq <- "CDR3.nucleotide.sequence"
  aa.seq <- "CDR3.amino.acid.sequence"
  .count <- "Read.count"
  vgenes <- "V.gene"
  jgenes <- "J.gene"
  dgenes <- "D.gene"
  vend <- "V.end"
  jstart <- "J.start"
  dstart <- "D5.end"
  dend <- "D3.end"
  vd.insertions <- "VD.insertions"
  dj.insertions <- "DJ.insertions"
  total.insertions <- "Total.insertions"
  .skip <- 0
  .sep <- "\t"

  if (substr(l, 1, 2) != "NA") {
    .count <- "Umi.count"
  }

  parse_repertoire(
    .filename = .filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
    .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
    .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
    .total.insertions = total.insertions, .skip = .skip, .sep = .sep
  )
}

parse_vdjtools <- function(.filename) {
  skip <- 0

  # Check for different VDJtools outputs
  f <- file(.filename, "r")
  l <- readLines(f, 1)
  close(f)

  .skip <- 0
  .count <- "count"
  filename <- .filename
  nuc.seq <- "cdr3nt"
  aa.seq <- "CDR3aa"
  vgenes <- "V"
  jgenes <- "J"
  dgenes <- "D"
  vend <- "Vend"
  jstart <- "Jstart"
  dstart <- "Dstart"
  dend <- "Dend"
  vd.insertions <- NA
  dj.insertions <- NA
  total.insertions <- NA
  .sep <- "\t"

  if (length(strsplit(l, "-", TRUE)) > 0) {
    if (length(strsplit(l, "-", TRUE)[[1]]) == 3) {
      if (strsplit(l, "-", TRUE)[[1]][2] == "header") {
        .count <- "count"
        .skip <- 1
      }
    } else if (tolower(substr(l, 1, 2)) == "#s") {
      .count <- "#Seq. count"
      nuc.seq <- "N Sequence"
      aa.seq <- "AA Sequence"
      vgenes <- "V segments"
      jgenes <- "J segments"
      dgenes <- "D segment"
      vend <- NA
      jstart <- NA
      dstart <- NA
      dend <- NA
    } else if (stringr::str_detect(l, "#")) {
      .count <- "X.count"
    } else {
      .count <- "count"
    }
  }

  parse_repertoire(
    .filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
    .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
    .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
    .total.insertions = total.insertions, .skip = .skip, .sep = .sep
  )
}

parse_imgt <- function(.filename) {
  .fix.imgt.alleles <- function(.col) {
    sapply(strsplit(.col, " "), function(x) {
      if (length(x) > 1) {
        x[[2]]
      } else {
        NA
      }
    })
  }

  f <- file(.filename, "r")
  l <- readLines(f, 2)[2]
  close(f)

  nuc.seq <- "JUNCTION"
  aa.seq <- NA
  .count <- NA
  vgenes <- "V-GENE and allele"
  jgenes <- "J-GENE and allele"
  dgenes <- "D-GENE and allele"
  vend <- "3'V-REGION end"
  jstart <- "5'J-REGION start"
  dstart <- "D-REGION start"
  dend <- "D-REGION end"
  vd.insertions <- NA
  dj.insertions <- NA
  total.insertions <- NA
  .skip <- 0
  .sep <- "\t"
  junc_start <- .make_names("JUNCTION start")

  df <- parse_repertoire(
    .filename = .filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
    .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
    .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
    .total.insertions = total.insertions, .skip = .skip, .sep = .sep, .add = junc_start
  )

  df[[IMMCOL$ve]] <- df[[IMMCOL$ve]] - df[[junc_start]]
  df[[IMMCOL$ds]] <- df[[IMMCOL$ds]] - df[[junc_start]]
  df[[IMMCOL$de]] <- df[[IMMCOL$de]] - df[[junc_start]]
  df[[IMMCOL$js]] <- df[[IMMCOL$js]] - df[[junc_start]]

  df[[IMMCOL$ve]][df[[IMMCOL$ve]] < 0] <- NA
  df[[IMMCOL$ds]][df[[IMMCOL$ds]] < 0] <- NA
  df[[IMMCOL$de]][df[[IMMCOL$de]] < 0] <- NA
  df[[IMMCOL$js]][df[[IMMCOL$js]] < 0] <- NA

  df[[IMMCOL$v]] <- .fix.imgt.alleles(df[[IMMCOL$v]])
  df[[IMMCOL$d]] <- .fix.imgt.alleles(df[[IMMCOL$d]])
  df[[IMMCOL$j]] <- .fix.imgt.alleles(df[[IMMCOL$j]])

  df[[junc_start]] <- NULL

  df
}

# parse_vidjil <- function (.filename) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }
#
# parse_rtcr <- function (.filename) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }
#
# parse_imseq <- function (.filename) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }

parse_airr <- function(.filename) {
  df <- airr::read_rearrangement(.filename)

  df <- df %>%
    select(
      sequence, v_call, d_call, j_call, junction, junction_aa,
      contains("v_germline_end"), contains("d_germline_start"), contains("d_germline_end"),
      contains("j_germline_start"), contains("np1_length"), contains("np2_length"),
      contains("duplicate_count")
    )

  namekey <- c(
    duplicate_count = IMMCOL$count, junction = IMMCOL$cdr3nt, junction_aa = IMMCOL$cdr3aa,
    v_call = IMMCOL$v, d_call = IMMCOL$d, j_call = IMMCOL$j, v_germline_end = IMMCOL$ve,
    d_germline_start = IMMCOL$ds, d_germline_end = IMMCOL$de, j_germline_start = IMMCOL$js,
    np1_length = "unidins", np2_length = IMMCOL$dnj, sequence = IMMCOL$seq
  )

  names(df) <- namekey[names(df)]

  if (!("unidins" %in% colnames(df))) {
    df["unidins"] <- NA
  }

  recomb_type <- .which_recomb_type(df[[IMMCOL$v]])

  if (!is.na(recomb_type)) {
    if (recomb_type == "VJ") {
      df[IMMCOL$vnj] <- df["unidins"]
      df[IMMCOL$vnd] <- NA
      df[IMMCOL$dnj] <- NA
    } else if (recomb_type == "VDJ") {
      df[IMMCOL$vnj] <- NA
      df[IMMCOL$vnd] <- df["unidins"]
    }
  }

  for (column in IMMCOL$order) {
    if (!(column %in% colnames(df))) {
      df[column] <- NA
    }
  }

  df <- df[IMMCOL$order]
  total <- sum(df$Clones)
  df[IMMCOL$prop] <- df[IMMCOL$count] / total
  df[IMMCOL$seq] <- stringr::str_remove_all(df[[IMMCOL$seq]], "N")
  df <- .postprocess(df)
  df
}

parse_immunarch <- function(.filename) {
  df <- readr::read_tsv(.filename, col_types = cols(), comment = "#")
  if (ncol(df) == 1) {
    # "," in the files, parse differently then
    df <- readr::read_csv(.filename, col_types = cols(), comment = "#")
  }
  df <- .postprocess(df)
  df
}

parse_10x_consensus <- function(.filename) {
  df <- parse_repertoire(.filename,
    .nuc.seq = "cdr3_nt", .aa.seq = NA, .count = "umis",
    .vgenes = "v_gene", .jgenes = "j_gene", .dgenes = "d_gene",
    .vend = NA, .jstart = NA, .dstart = NA, .dend = NA,
    .vd.insertions = NA, .dj.insertions = NA, .total.insertions = NA,
    .skip = 0, .sep = ",", .add = c("chain", "clonotype_id", "consensus_id")
  )
  setnames(df, "clonotype_id", "ClonotypeID")
  setnames(df, "consensus_id", "ConsensusID")
  df
}

parse_10x_filt_contigs <- function(.filename) {
  df <- parse_repertoire(.filename,
    .nuc.seq = "cdr3_nt", .aa.seq = NA, .count = "umis",
    .vgenes = "v_gene", .jgenes = "j_gene", .dgenes = "d_gene",
    .vend = NA, .jstart = NA, .dstart = NA, .dend = NA,
    .vd.insertions = NA, .dj.insertions = NA, .total.insertions = NA,
    .skip = 0, .sep = ",", # .add = c("chain", "raw_clonotype_id", "raw_consensus_id", "barcode", "contig_id")
    .add = c("chain", "barcode", "contig_id")
  )
  setnames(df, "contig_id", "ContigID")
  # setnames(df, "raw_clonotype_id", "RawClonotypeID")
  # setnames(df, "raw_consensus_id", "RawConsensusID")
  setnames(df, "barcode", "Barcode")
  df
}

parse_archer <- function(.filename) {
  parse_repertoire(.filename,
    .nuc.seq = "Clonotype Sequence", .aa.seq = NA, .count = "Clone Abundance",
    .vgenes = "Predicted V Region", .jgenes = "Predicted J Region", .dgenes = "Predicted D Region",
    .vend = NA, .jstart = NA, .dstart = NA, .dend = NA,
    .vd.insertions = NA, .dj.insertions = NA, .total.insertions = NA,
    .skip = 0, .sep = "\t"
  )
}

##### Savers #####

save_immunarch <- function(.data, .path, .compress = TRUE) {
  if (.compress) {
    filepath <- gzfile(paste0(.path, ".tsv.gz"), compression = 9)
  } else {
    filepath <- paste0(.path, ".tsv")
  }
  readr::write_lines(paste0("# Exported from immunarch ", packageVersion("immunarch"), " https://immunarch.com"),
    path = filepath
  )
  filepath <- gzfile(paste0(.path, ".tsv.gz"), compression = 9)
  readr::write_tsv(x = .data, path = filepath, append = TRUE, col_names = TRUE)
}

save_vdjtools <- function(.data, .path, .compress = TRUE) {
  if (.compress) {
    filepath <- gzfile(paste0(.path, ".tsv.gz"), compression = 9)
  } else {
    filepath <- paste0(.path, ".tsv")
  }

  old <- c(
    IMMCOL$count,
    IMMCOL$prop,
    IMMCOL$cdr3nt,
    IMMCOL$cdr3aa,
    IMMCOL$v,
    IMMCOL$d,
    IMMCOL$j
  )

  new <- c(
    "#Seq. Count",
    "Percent",
    "N Sequence",
    "AA Sequence",
    "V Segments",
    "D Segment",
    "J Segments"
  )

  names(.data) <- plyr::mapvalues(names(.data),
    from = old,
    to = as.character(new)
  )

  readr::write_tsv(x = .data, path = filepath)
}
