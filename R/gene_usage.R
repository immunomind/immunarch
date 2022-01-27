if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("alias", "species", "gene", "Quant", "Names", "Counts"))
}


#' Main function for estimation of V-gene and J-gene statistics
#'
#' @concept gene_usage
#'
#' @aliases geneUsage get_aliases get_genes
#'
#' @description
#' An utility function to analyse the immune receptor gene usage
#' (IGHD, IGHJ, IDHV, IGIJ, IGKJ, IGKV, IGLJ, IGLV, TRAJ, TRAV, TRBD, etc.)
#' and statistics. For gene details run \code{gene_stats()}.
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
#' @param .gene A character vector of length one with the name of the gene you want
#' to analyse of the specific species. If you provide a vector of different length, only first element
#' will be used. The string should also contain the species of interest, for example, valid ".gene" arguments
#' are "hs.trbv", "HomoSapiens.TRBJ" or "macmul.IGHV". For details run \code{gene_stats()}.
#' @param .quant Select the column with data to evaluate.
#' Pass NA if you want to compute gene statistics at the clonotype level without re-weighting.
#' Pass "count" to use the "Clones" column to weight genes by abundance of their corresponding clonotypes.
#' @param .ambig An option to handle ambiguous gene assigments, e.g., "TRAV1,TRAV2".
#'
#'
#' - Pass "inc" to include all possible gene segments, so "TRAV1,TRAV2" is counted as a different gene segment.
#'
#' - Pass "exc" to exclude all ambiguous gene assignments, so "TRAV1,TRAV2" is excluded from the resultant gene table.
#'
#'
#' We recommend to turn in on by passing "inc" (turned on by default).
#' You can exclude data for the cases where
#' there is no clear match for gene, include it for every supplied gene,
#' or pick only first from the set. Set it to "exc", "inc" or "maj", correspondingly.
#' @param .type Set the type of data to evaluate: "segment", "allele", or "family".
#' @param .norm If TRUE then return proportions of genes. If FALSE then return counts of genes.
#'
#' @return
#' A data frame with rows corresponding to gene segments and columns corresponding to the input samples.
#'
#' @examples
#' data(immdata)
#' gu <- geneUsage(immdata$data)
#' vis(gu)
#' @export geneUsage get_genes
geneUsage <- function(.data, # df, list, MonetDB
                      .gene = c("hs.trbv", "HomoSapiens.TRBJ", "macmul.IGHV"),
                      .quant = c(NA, "count"),
                      .ambig = c("inc", "exc", "maj"),
                      .type = c("segment", "allele", "family"),
                      .norm = FALSE) {
  .type <- .type[1]
  .ambig <- .ambig[1]
  .quant <- .quant[1]
  .gene <- .gene[1]

  if (has_class(.data, "list")) {
    if (length(.data) == 1) {
      return(geneUsage(.data[[1]], .gene = .gene, .quant = .quant, .ambig = .ambig, .norm = .norm))
    } else {
      genedata <- lapply(.data, geneUsage,
        .gene = .gene, .quant = .quant,
        .ambig = .ambig, .type = .type, .norm = .norm
      )

      # Add check if list has 1 element
      for (i in 2:(length(genedata))) {
        if (i > 2) {
          allres <- full_join(allres, genedata[[i]], by = "Names")
          names(allres)[i + 1] <- names(genedata)[i]
        } else {
          allres <- full_join(genedata[[1]], genedata[[2]], by = "Names")
          names(allres)[2] <- names(genedata)[1]
          names(allres)[3] <- names(genedata)[2]
        }
      }
      allres <- allres[order(allres$Names), ]

      return(add_class(allres, "immunr_gene_usage"))
    }
  }

  "%!in%" <- function(x, y) !("%in%"(x, y))

  which_gene <- strsplit(.gene, ".", TRUE)[[1]][2]
  gene_col <- tolower(substr(which_gene, 4, 4))
  if (gene_col == "v") {
    gene_col <- IMMCOL$v
  } else if (gene_col == "d") {
    gene_col <- IMMCOL$d
  } else if (gene_col == "j") {
    gene_col <- IMMCOL$j
  } else {
    stop("The entered gene_col name is invalid")
  }

  gene_col_new_name <- paste0(gene_col, ".no.alleles")
  .data %<>% add_column_without_alleles(gene_col, gene_col_new_name)
  gene_col <- gene_col_new_name

  .gene <- get_genes(.gene, .type)

  # Bad Cases:
  # X.name is segments, .type == "alleles"
  # X.name is families, .type == "segments"
  # X.name is families, .type == "alleles"

  if (.ambig == "inc") {
    if (.type == "segment") {
      .data[[gene_col]] <- return_segments(.data[[gene_col]])
    } else if (.type == "family") {
      .data[[gene_col]] <- return_families(.data[[gene_col]])
    }
  }

  if (has_class(.data, "data.table")) {
    .data <- .data %>% lazy_dt()
  }
  dataset <- .data %>%
    select(Gene = gene_col, IMMCOL$count) %>%
    group_by(Gene) %>%
    rename(Quant = IMMCOL$count)

  if (is.na(.quant)) {
    dataset <- dataset %>%
      summarise(count = n())
  } else {
    dataset <- dataset %>%
      summarise(count = sum(Quant))
  }

  dataset <- dataset %>%
    collect(n = Inf)

  if (.ambig == "inc") {
    res <- dataset
    names(res) <- c("Names", IMMCOL$count)
  } else {
    gene_col_split <- strsplit(dataset[["Gene"]], ", ", fixed = TRUE, useBytes = TRUE)
    counts <- dataset$count
    names(counts) <- dataset[["Gene"]]

    if (.ambig == "exc") {
      included_names <- names(counts) %in% gene_col_split[sapply(gene_col_split, length) == 1]
      chosen <- counts[included_names]
      chosen <- chosen[names(chosen) %in% .gene]
      other <- counts[names(chosen) %!in% .gene]
      res <- data.frame(
        Names = names(chosen),
        Count = unlist(chosen), stringsAsFactors = FALSE
      )
    } else if (.ambig == "maj") {
      maj_counts <- data.frame(Names = sapply(gene_col_split, "[[", 1), Counts = counts, stringsAsFactors = FALSE)
      row.names(maj_counts) <- NULL
      res <- maj_counts %>%
        group_by(Names) %>%
        summarise(Counts = sum(Counts))
    }
  }
  if (length(.gene %in% res[[1]]) < length(.gene)) {
    missing <- data.frame(.gene[!.gene %in% res[[1]]],
      0,
      stringsAsFactors = FALSE
    )
    names(missing) <- names(res)
    res <- rbind(res, missing)
  }
  row.names(res) <- NULL
  res <- res[order(res[[1]]), ]
  if (.norm) {
    res[, 2] <- res[, 2] / sum(res[, 2])
  }
  res <- as_tibble(as.data.frame(res))
  add_class(res, "immunr_gene_usage")
}

#' WIP
#'
#' @concept gene_usage
#'
#' @importFrom dplyr n
#'
#' @aliases gene_stats
#'
#' @return
#' \code{gene_stats} returns all segment gene statistics
#'
#' @examples
#' gene_stats()
#' get_genes("hs.trbv", "segment")
#' @export gene_stats
gene_stats <- function() {
  res <- GENE_SEGMENTS %>%
    group_by(alias, species, gene) %>%
    summarise(n = n()) %>%
    reshape2::dcast(alias + species ~ gene, value.var = "n")
  res[is.na(res)] <- 0
  res
}

get_genes <- function(.gene = c("hs.trbv", "HomoSapiens.trbj", "macmul.ighv"), .type = c("segment", "allele", "family")) {
  .gene <- .gene[1]
  .type <- .type[1]

  # Species
  .gene <- tolower(.gene)
  which_species <- strsplit(.gene, ".", TRUE)[[1]][1]
  if (which_species %in% GENE_SEGMENTS$alias) {
    which_species_col <- "alias"
  } else {
    if (which_species %in% tolower(GENE_SEGMENTS$species)) {
      which_species_col <- "species"
    } else {
      stop("Unknown species name / alias")
    }
  }

  # Genes
  which_gene <- strsplit(.gene, ".", TRUE)[[1]][2]
  if (!(which_gene %in% GENE_SEGMENTS[["gene"]][tolower(GENE_SEGMENTS[[which_species_col]]) == which_species])) {
    stop("Unknown gene name")
  }

  # Id's type
  which_type <- paste0(.type, "_id")

  sort(unique(GENE_SEGMENTS[[which_type]][tolower(GENE_SEGMENTS[[which_species_col]]) == which_species & GENE_SEGMENTS[["gene"]] == which_gene]))
}
