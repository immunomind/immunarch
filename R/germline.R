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
#' repGermline(.data, .species, .min_nuc_outside_cdr3, .threads)
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
#' @param .min_nuc_outside_cdr3 This parameter sets how many nucleotides should have V or J chain
#' outside of CDR3 to be considered good for further alignment.
#'
#' @param .threads Number of threads to use.
#'
#' @return
#'
#' Data with added columns:
#' * Sequence (FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4 in nucleotides; column will be replaced if exists)
#' * V.allele, J.allele (chosen alleles of V and J genes),
#' * V.aa, J.aa (V and J sequences from original clonotype, outside CDR3, converted to amino acids)
#' * CDR3.length (length of CDR3),
#' * Germline.sequence (combined germline nucleotide sequence)
#'
#' @examples
#'
#' data(bcrdata)
#'
#' bcrdata$data %>%
#'   top(5) %>%
#'   repGermline()
#' @export repGermline
repGermline <- function(.data,
                        .species = "HomoSapiens",
                        .min_nuc_outside_cdr3 = 5,
                        .threads = 1) {
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
      min_nuc_outside_cdr3 = .min_nuc_outside_cdr3,
      threads = .threads
    )
  return(.data)
}

germline_single_df <- function(data,
                               reference,
                               species,
                               min_nuc_outside_cdr3,
                               threads,
                               sample_name = NA) {
  data %<>%
    validate_mandatory_columns(sample_name) %>%
    add_allele_column(reference[["allele_id"]], "V") %>%
    merge_reference_sequences(reference, "V", species, sample_name) %>%
    add_allele_column(reference[["allele_id"]], "J") %>%
    merge_reference_sequences(reference, "J", species, sample_name) %>%
    validate_chains_length(min_nuc_outside_cdr3, sample_name) %>%
    calculate_germlines_parallel(threads, sample_name) %>%
    filter(!is.na(get("Germline.sequence")))
  return(data)
}

calculate_germlines_parallel <- function(data, threads, sample_name) {
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
    calculate_new_columns(row, sample_name)
  }) %>%
    map_dfr(~.) %>%
    germline_handle_warnings() %>%
    merge_germline_results(data)

  if (!has_no_data(cluster)) {
    stopCluster(cluster)
  }
  return(data)
}

calculate_new_columns <- function(row, sample_name) {
  v_ref <- row[["V.ref.nt"]]
  j_ref <- row[["J.ref.nt"]]
  cdr1_nt <- row[["CDR1.nt"]]
  cdr2_nt <- row[["CDR2.nt"]]
  cdr3_nt <- row[["CDR3.nt"]]
  fr1_nt <- row[["FR1.nt"]]
  fr2_nt <- row[["FR2.nt"]]
  fr3_nt <- row[["FR3.nt"]]
  fr4_nt <- row[["FR4.nt"]]

  if (any(is.na(c(
    v_ref, j_ref, cdr1_nt, cdr2_nt, fr1_nt, fr2_nt, fr3_nt, cdr3_nt, fr4_nt
  )))) {
    # warnings cannot be displayed from parApply; save them and display after finish
    warn <- paste0(
      "Some of mandatory fields in a row ",
      optional_sample("from sample ", sample_name, " "),
      "contain unexpected NAs! Found values:\n",
      "V.ref.nt = ",
      v_ref,
      "\nJ.ref.nt = ",
      j_ref,
      "\nCDR1.nt = ",
      cdr1_nt,
      "\nCDR2.nt = ",
      cdr2_nt,
      "\nCDR3.nt = ",
      cdr3_nt,
      "\nFR1.nt = ",
      fr1_nt,
      "\nFR2.nt = ",
      fr2_nt,
      "\nFR3.nt = ",
      fr3_nt,
      "\nFR4.nt = ",
      fr4_nt,
      "\nThe row will be dropped!"
    )
    return(list(
      Sequence = NA,
      V.aa = NA,
      J.aa = NA,
      CDR3.length = NA,
      Germline.sequence = NA,
      Warning = warn
    ))
  } else {
    seq <- paste0(fr1_nt, cdr1_nt, fr2_nt, cdr2_nt, fr3_nt, cdr3_nt, fr4_nt) %>% toupper()
    cdr1_aa <- ifelse(has_no_data(row[["CDR1.aa"]]), bunch_translate(cdr1_nt), row[["CDR1.aa"]])
    cdr2_aa <- ifelse(has_no_data(row[["CDR2.aa"]]), bunch_translate(cdr2_nt), row[["CDR2.aa"]])
    fr1_aa <- ifelse(has_no_data(row[["FR1.aa"]]), bunch_translate(fr1_nt), row[["FR1.aa"]])
    fr2_aa <- ifelse(has_no_data(row[["FR2.aa"]]), bunch_translate(fr2_nt), row[["FR2.aa"]])
    fr3_aa <- ifelse(has_no_data(row[["FR3.aa"]]), bunch_translate(fr3_nt), row[["FR3.aa"]])
    fr4_aa <- ifelse(has_no_data(row[["FR4.aa"]]), bunch_translate(fr4_nt), row[["FR4.aa"]])
    v_aa <- paste0(fr1_aa, cdr1_aa, fr2_aa, cdr2_aa, fr3_aa)
    j_aa <- fr4_aa

    # trim intersection of V and CDR3 from reference V gene
    v_length <- str_length(cdr1_nt) + str_length(cdr2_nt) +
      str_length(fr1_nt) + str_length(fr2_nt) + str_length(fr3_nt)
    v_part <- str_sub(v_ref, 1, v_length)

    cdr3_length <- str_length(cdr3_nt)
    cdr3_part <- paste(rep("n", cdr3_length), collapse = "")

    # trim intersection of J and CDR3 from reference J gene
    j_length <- str_length(fr4_nt)
    j_part <- str_sub(j_ref, str_length(j_ref) - j_length + 1, str_length(j_ref))

    germline <- paste0(v_part, cdr3_part, j_part) %>%
      toupper()

    # return values for new calculated columns
    return(list(
      Sequence = seq,
      V.aa = v_aa,
      J.aa = j_aa,
      CDR3.length = cdr3_length,
      Germline.sequence = germline,
      Warning = NA
    ))
  }
}

add_allele_column <- function(.data, .reference_allele_ids, .gene) {
  raw_genes_colname <- paste0(.gene, ".name")
  target_colname <- paste0(.gene, ".allele")
  if (validate_columns(.data, raw_genes_colname, target_colname)) {
    .data[[target_colname]] <- .data[[raw_genes_colname]] %>% sapply(
      function(genes_string) {
        # first allele is substring until first ',' or '(' in string taken from column with gene names
        first_allele <- strsplit(genes_string, ",|\\(")[[1]][1]
        # MiXCR uses *00 for unknown alleles; drop it to find first matching allele in reference
        name_to_search <- str_replace(first_allele, fixed("*00"), "")
        first_matching_allele <- .reference_allele_ids[which(startsWith(
          .reference_allele_ids,
          name_to_search
        ))][1]
        if (has_no_data(first_matching_allele)) {
          return(name_to_search)
        } else {
          return(first_matching_allele)
        }
      }
    )
  }
  return(.data)
}

merge_reference_sequences <- function(data, reference, chain_letter, species, sample_name) {
  chain_seq_colname <- paste0(chain_letter, ".ref.nt")
  chain_allele_colname <- paste0(chain_letter, ".allele")
  colnames(reference) <- c(chain_seq_colname, chain_allele_colname)

  # check for alleles in data that don't exist in the reference
  alleles_in_data <- unique(data[[chain_allele_colname]])
  alleles_in_ref <- unique(reference[[chain_allele_colname]])
  missing_alleles <- alleles_in_data[!(alleles_in_data %in% alleles_in_ref)]
  if (length(missing_alleles) > 0) {
    warning(
      "Genes or alleles ",
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

validate_mandatory_columns <- function(data, sample_name) {
  if (nrow(data) == 0) {
    stop(
      "Sample ",
      optional_sample("", sample_name, " "),
      "dataframe is empty!"
    )
  }
  old_length <- nrow(data)
  mandatory_columns <- c("FR1.nt", "CDR1.nt", "FR2.nt", "CDR2.nt", "FR3.nt", "CDR3.nt", "FR4.nt")
  for (column in mandatory_columns) {
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
        "Dropped all rows when filtering out NAs from mandatory columns ",
        paste(mandatory_columns, collapse = ", "),
        optional_sample(" in sample ", sample_name, ""),
        "!"
      )
    }
    data %<>% filter(!is.na(get(column)))
  }
  dropped_num <- old_length - nrow(data)
  if (dropped_num > 0) {
    warning(
      dropped_num,
      " rows from ",
      old_length,
      optional_sample(" in sample ", sample_name, ""),
      " were dropped because of missing values in mandatory columns ",
      paste(mandatory_columns, collapse = ", "),
      "!"
    )
  }
  return(data)
}

validate_chains_length <- function(data, min_nuc_outside_cdr3, sample_name) {
  old_length_v <- nrow(data)
  data %<>% filter(
    str_length(get("FR1.nt")) + str_length(get("CDR1.nt")) + str_length(get("FR2.nt"))
      + str_length(get("CDR2.nt")) + str_length(get("FR3.nt"))
    >= min_nuc_outside_cdr3
  )
  dropped_v <- old_length_v - nrow(data)
  old_length_j <- nrow(data)
  if (nrow(data) > 0) {
    data %<>% filter(str_length(get("FR4.nt")) >= min_nuc_outside_cdr3)
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

merge_germline_results <- function(new_columns, data) {
  data %<>%
    subset(select = -c(get("Sequence"), get("V.ref.nt"), get("J.ref.nt"))) %>%
    cbind(new_columns)
  return(data)
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
