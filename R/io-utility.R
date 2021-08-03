.remove.ext <- function(.str) {
  # gsub(pattern = '.*/|[.].*$', replacement = '', x = .str)
  gsub(pattern = ".*/|[.](txt|tsv|csv)$|([.](txt|tsv|csv))?[.](gz|bzip|bzip2|bz2)$", replacement = "", x = .str)
}


.detect_format <- function(.filename) {
  res_format <- NA

  f <- file(.filename, "r")
  l <- readLines(f, 1)
  # use 2nd line of file for JSON formats
  l2 <- ifelse(str_trim(l) == "{", readLines(f, 1), NA)
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
  } else if (str_detect(l, "AAseq") && str_detect(l, "Vregion") && str_detect(l, "Frequency")) {
    res_format <- "catt"
  } else if (str_detect(l, "Number of reads") && str_detect(l, "Amino acid sequence") && str_detect(l, "V gene")) {
    res_format <- "rtcr"
  } else if (str_detect(l, "seqId") && str_detect(l, "cdrNucSeq") && str_detect(l, "cdrAASeq")) {
    res_format <- "imseq"
  } else if (str_trim(l2) == "\"clones\": [") {
    res_format <- "vidjil"
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

.postprocess <- function(.data, .mode) {
  .data[[IMMCOL$cdr3nt]][.data[[IMMCOL$cdr3nt]] == "NONE"] <- NA
  logic <- is.na(.data[[IMMCOL$cdr3aa]]) & !is.na(.data[[IMMCOL$cdr3nt]])
  if (any(logic)) {
    .data[[IMMCOL$cdr3aa]][logic] <- bunch_translate(.data[[IMMCOL$cdr3nt]][logic])
  }

  logic <- is.na(.data[[IMMCOL$cdr3aa]]) & is.na(.data[[IMMCOL$cdr3nt]])
  if (any(logic)) {
    warn_msg <- c("  [!] Removed ", sum(logic))
    warn_msg <- c(warn_msg, " clonotypes with no nucleotide and amino acid CDR3 sequence.")
    message(warn_msg)
  }
  .data <- .data[!logic, ]

  if (nrow(.data)) {
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
