#' Creates germlines for clonal lineages
#'
#' @concept germline
#'
#' @aliases repGermline
#'
#' @importFrom stringr str_sub str_length str_replace fixed str_extract_all str_extract boundary str_c
#' @importFrom purrr imap map_dfr
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom dplyr filter rowwise
#' @importFrom parallel parApply detectCores makeCluster clusterExport stopCluster
#' @importFrom ape as.DNAbin clustal
#'
#' @description This function creates germlines for clonal lineages. B cell clonal lineage
#' represents a set of B cells that presumably have a common origin (arising from the same VDJ
#' rearrangement event) and a common ancestor. Each clonal lineage has its own germline sequence
#' that represents the ancestral sequence for each BCR in clonal lineage. In other words,
#' germline sequence is a sequence of B-cells immediately after VDJ recombination, before
#' B-cell maturation and hypermutation process. Germline sequence is useful for assessing
#' the degree of mutation and maturity of the repertoire.
#'
#' @usage
#'
#' repGermline(.data,
#' .species, .align_j_gene, .min_nuc_outside_cdr3, .threads)
#'
#' @param .data The data to be processed. Can be \link{data.frame}, \link{data.table}
#' or a list of these objects.
#'
#' It must have columns in the immunarch compatible format \link{immunarch_data_format}.
#'
#' @param .species Species from which the data was acquired. Available options:
#' "HomoSapiens" (default), "MusMusculus", "BosTaurus", "CamelusDromedarius",
#' "CanisLupusFamiliaris", "DanioRerio", "MacacaMulatta", "MusMusculusDomesticus",
#' "MusMusculusCastaneus", "MusMusculusMolossinus", "MusMusculusMusculus", "MusSpretus",
#' "OncorhynchusMykiss", "OrnithorhynchusAnatinus", "OryctolagusCuniculus", "RattusNorvegicus",
#' "SusScrofa".
#'
#' @param .align_j_gene MiXCR provides the number of J indels only for 1 allele of J gene
#' in the output file, and a germline can contain another allele. Therefore, calculation of
#' J gene start in reference based on numbers from input file can be sometimes incorrect.
#' As result, J gene in the germline will be trimmed in the start or will contain some
#' nucleotides from CDR3. Setting this parameter to TRUE enables precise alignment of J genes
#' to detect the correct starting point, but it significantly reduces performance.
#'
#' @param .min_nuc_outside_cdr3 This parameter sets how many nucleotides should have V or J chain
#' outside of CDR3 to be considered good for further alignment.
#'
#' @param .threads Number of threads to use.
#'
#' @return
#'
#' Data with added columns:
#' * V.first.allele, J.first.allele (first alleles of V and J genes),
#' * V.ref.nt, J.ref.nt (V and J reference sequences),
#' * V.germline.nt, J.germline.nt (V and J germline sequences; they are references with
#'   trimmed parts that are from CDR3),
#' * CDR3.germline.length (length of CDR3 in the germline),
#' * Germline.sequence (combined germline sequence)
#'
#' @examples
#'
#' data(bcrdata)
#'
#' bcrdata$data %>%
#'   top(5) %>%
#'   repGermline(.threads = 1)
#' @export repGermline
repGermline <- function(.data,
                        .species = "HomoSapiens",
                        .align_j_gene = FALSE,
                        .min_nuc_outside_cdr3 = 5,
                        .threads = parallel::detectCores()) {
  if (.align_j_gene) {
    require_system_package("clustalw", error_message = paste0(
      "repGermline with .align_j_gene = TRUE requires Clustal W app to be installed!\n",
      "Please download it from here: http://www.clustal.org/download/current/\n",
      "or install it with your system package manager (such as apt or dnf)."
    ))
  }

  # prepare reference sequences for all alleles
  genesegments_env <- new.env()
  data("genesegments", envir = genesegments_env)
  reference <- genesegments_env$GENE_SEGMENTS %>% filter(species == .species)
  rm(genesegments_env)
  reference <- reference[c("sequence", "allele_id")]

  .data %<>%
    apply_to_sample_or_list(
      germline_single_df,
      .with_names = TRUE,
      reference = reference,
      species = .species,
      align_j_gene = .align_j_gene,
      min_nuc_outside_cdr3 = .min_nuc_outside_cdr3,
      threads = .threads
    )
  return(.data)
}

germline_single_df <- function(data,
                               reference,
                               species,
                               align_j_gene,
                               min_nuc_outside_cdr3,
                               threads,
                               sample_name = NA) {
  data %<>%
    validate_genes_edges(sample_name) %>%
    add_column_with_first_gene(
      "V.name",
      "V.first.allele",
      .with_allele = TRUE
    ) %>%
    merge_reference_sequences(reference, "V", species, sample_name) %>%
    add_column_with_first_gene(
      "J.name",
      "J.first.allele",
      .with_allele = TRUE
    ) %>%
    merge_reference_sequences(reference, "J", species, sample_name) %>%
    validate_chains_length(min_nuc_outside_cdr3, sample_name) %>%
    calculate_germlines_parallel(align_j_gene, threads, sample_name) %>%
    filter(!is.na(get("Germline.sequence")))
  return(data)
}

calculate_germlines_parallel <- function(data, align_j_gene, threads, sample_name) {
  if (threads == 1) {
    cluster <- NA
  } else {
    cluster <- makeCluster(threads)
    clusterExport(cluster, c("generate_germline_sequence", "align_and_find_j_start", "sample_name"),
      envir = environment()
    )
  }

  # rowwise parallel calculation of new columns that are added to data
  data <- par_or_normal_apply(cluster, data, 1, function(row) {
    generate_germline_sequence(
      seq = row[["Sequence"]],
      v_ref = row[["V.ref.nt"]],
      j_ref = row[["J.ref.nt"]],
      cdr1_nt = row[["CDR1.nt"]],
      cdr2_nt = row[["CDR2.nt"]],
      fr1_nt = row[["FR1.nt"]],
      fr2_nt = row[["FR2.nt"]],
      fr3_nt = row[["FR3.nt"]],
      fr4_nt = row[["FR4.nt"]],
      cdr3_start = as.integer(row[["CDR3.start"]]),
      cdr3_end = as.integer(row[["CDR3.end"]]),
      j_start = as.integer(row[["J.start"]]),
      j3_del = as.integer(row[["J3.Deletions"]]),
      align_j_gene = align_j_gene,
      sample_name = sample_name
    )
  }) %>%
    map_dfr(~.) %>%
    germline_handle_warnings() %>%
    cbind(data, .)

  if (!has_no_data(cluster)) {
    stopCluster(cluster)
  }
  return(data)
}

generate_germline_sequence <- function(seq, v_ref, j_ref, cdr1_nt, cdr2_nt,
                                       fr1_nt, fr2_nt, fr3_nt, fr4_nt,
                                       cdr3_start, cdr3_end, j_start, j3_del,
                                       align_j_gene, sample_name) {
  if (any(is.na(c(
    seq, v_ref, j_ref, cdr1_nt, cdr2_nt, fr1_nt, fr2_nt, fr3_nt, fr4_nt,
    cdr3_start, cdr3_end, j_start, j3_del
  ))) ||
    (seq == "")) {
    # warnings cannot be displayed from parApply; save them and display after finish
    warn <- paste0(
      "Some of mandatory fields in a row ",
      optional_sample("from sample ", sample_name, " "),
      "contain unexpected NA or empty strings! Found values:\n",
      "Sequence = ",
      seq,
      ",\nV.ref.nt = ",
      v_ref,
      ",\nJ.ref.nt = ",
      j_ref,
      ",\nCDR1.nt = ",
      cdr1_nt,
      ",\nCDR2.nt = ",
      cdr2_nt,
      ",\nFR1.nt = ",
      fr1_nt,
      ",\nFR2.nt = ",
      fr2_nt,
      ",\nFR3.nt = ",
      fr3_nt,
      ",\nFR4.nt = ",
      fr4_nt,
      ",\nCDR3.start = ",
      cdr3_start,
      ", CDR3.end = ",
      cdr3_end,
      ", J.start = ",
      j_start,
      ", J3.Deletions = ",
      j3_del,
      ".\nThe row will be dropped!"
    )
    return(list(
      V.germline.nt = NA,
      J.germline.nt = NA,
      CDR3.germline.length = NA,
      Germline.sequence = NA,
      Warning = warn
    ))
  } else {
    v_end <- str_length(cdr1_nt) + str_length(cdr2_nt) +
      str_length(fr1_nt) + str_length(fr2_nt) + str_length(fr3_nt)
    cdr3_length <- cdr3_end - cdr3_start

    # trim intersection of V and CDR3 from reference V gene
    v_part <- str_sub(v_ref, 1, v_end)

    cdr3_part <- paste(rep("n", cdr3_length), collapse = "")

    # trim intersection of J and CDR3 from reference J gene
    if (align_j_gene) {
      calculated_j_start <- align_and_find_j_start(j_ref, fr4_nt)
    } else {
      calculated_j_start <- max(0, cdr3_end - j_start - j3_del + 1)
    }
    j_part <- str_sub(j_ref, calculated_j_start)

    germline <- paste0(v_part, cdr3_part, j_part) %>%
      toupper()

    # return values for new calculated columns
    return(list(
      V.germline.nt = v_part,
      J.germline.nt = j_part,
      CDR3.germline.length = cdr3_length,
      Germline.sequence = germline,
      Warning = NA
    ))
  }
}

merge_reference_sequences <- function(data, reference, chain_letter, species, sample_name) {
  chain_seq_colname <- paste0(chain_letter, ".ref.nt")
  chain_allele_colname <- paste0(chain_letter, ".first.allele")
  colnames(reference) <- c(chain_seq_colname, chain_allele_colname)

  # check for alleles in data that don't exist in the reference
  alleles_in_data <- unique(data[[chain_allele_colname]])
  alleles_in_ref <- unique(reference[[chain_allele_colname]])
  missing_alleles <- alleles_in_data[!(alleles_in_data %in% alleles_in_ref)]
  if (length(missing_alleles) > 0) {
    warning(
      "Alleles ",
      paste(missing_alleles, collapse = ", "),
      " ",
      optional_sample("from sample ", sample_name, " "),
      "not found in the reference and will be dropped!\n",
      "Probably, species argument is wrong (current value: ",
      species,
      ") or the data contains non-BCR genes."
    )
  }

  data %<>% merge(reference, by = chain_allele_colname)
  if (nrow(data) == 0) {
    stop(
      "After merging with reference, the data ",
      optional_sample("from sample ", sample_name, " "),
      "is empty.\n",
      "There were no valid alleles in the data!"
    )
  }
  return(data)
}

validate_genes_edges <- function(data, sample_name) {
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty!"
    )
  }
  for (column in c("V.end", "J.start")) {
    if (!(column %in% colnames(data))) {
      stop(
        "Missing mandatory ",
        column,
        " column",
        optional_sample(" in sample ", sample_name, ""),
        "!"
      )
    }
    if (all(is.na(data[, column]))) {
      stop(
        "No data in mandatory ",
        column,
        " column",
        optional_sample(" in sample ", sample_name, ""),
        "!"
      )
    }
  }
  old_length <- nrow(data)
  data %<>% filter(!is.na(get("V.end")) & !is.na(get("J.start")))
  dropped_num <- old_length - nrow(data)
  if (dropped_num > 0) {
    warning(
      dropped_num,
      " rows from ",
      old_length,
      optional_sample(" in sample ", sample_name, ""),
      " were dropped because of missing values in mandatory columns V.end and J.start!"
    )
  }
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty after dropping missing values!"
    )
  }
  return(data)
}

validate_chains_length <- function(data, min_nuc_outside_cdr3, sample_name) {
  old_length_v <- nrow(data)
  data %<>% filter(
    v_len_outside_cdr3(
      get("V.end"),
      get("CDR3.start")
    ) >= min_nuc_outside_cdr3
  )
  dropped_v <- old_length_v - nrow(data)
  old_length_j <- nrow(data)
  if (nrow(data) > 0) {
    data %<>% filter(
      j_len_outside_cdr3(
        get("Sequence"),
        get("J.start"),
        get("CDR3.end")
      ) >= min_nuc_outside_cdr3
    )
  }
  dropped_j <- old_length_j - nrow(data)

  warning_prefix <- paste0(
    " clonotype(s) ",
    optional_sample("in sample ", sample_name, " "),
    "were dropped because they have too short part of "
  )
  warning_suffix <- paste0(
    " chain that doesn't intersect with CDR3!\n",
    "Lengths less than ",
    min_nuc_outside_cdr3,
    " are considered too short."
  )
  if (dropped_v > 0) {
    warning(dropped_v, warning_prefix, "V", warning_suffix)
  }
  if (dropped_j > 0) {
    warning(dropped_j, warning_prefix, "J", warning_suffix)
  }
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty after dropping sequences with too short V and J chains!"
    )
  }
  return(data)
}

# align reference J gene and FR4 segment from clonotype to find start of J gene outside of CDR3
align_and_find_j_start <- function(j_ref, fr4_seq, max_len_diff = 10) {
  # max_len_diff is needed to prevent alignment of sequences that are very different in length;
  # we are only interested in the left side of alignment
  j_len <- str_length(j_ref)
  fr4_len <- str_length(fr4_seq)
  min_len <- min(j_len, fr4_len)
  max_len <- max(j_len, fr4_len)
  trim_len <- min(min_len + max_len_diff, max_len)
  j_trimmed <- str_sub(j_ref, 1, trim_len)
  fr4_trimmed <- str_sub(fr4_seq, 1, trim_len)

  # alignment will contain vector of 2 aligned strings where deletions are "-" characters
  alignment <- list(J = j_trimmed, FR4 = fr4_trimmed) %>%
    lapply(
      function(sequence) {
        sequence %>%
          str_extract_all(boundary("character")) %>%
          unlist()
      }
    ) %>%
    ape::as.DNAbin() %>%
    ape::clustal() %>%
    as.character() %>%
    apply(1, paste, collapse = "") %>%
    lapply(toupper)

  num_deletions <- str_length(stringr::str_extract(alignment[["FR4"]], pattern = "^(-+)"))
  if (is.na(num_deletions) | (str_length(alignment[["J"]]) <= num_deletions)) {
    return(1)
  } else {
    return(num_deletions + 1)
  }
}

germline_handle_warnings <- function(df) {
  warnings <- df$Warning
  warnings <- warnings[!is.na(warnings)]
  options(warning.length = 5000L)
  for (warn in warnings) {
    warning(warn)
  }
  df$Warning <- NULL
  return(df)
}
