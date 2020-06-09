data(genesegments, envir = environment())
data(aaproperties, envir = environment())



#' Gene segments table
#'
#' @concept data
#'
#' @docType data
#'
#' @name gene_segments
#'
#' @aliases gene_segments genes segments GENE_SEGMENTS
NULL


#' Tables with amino acid properties
#'
#' @concept data
#'
#' @docType data
#'
#' @name aa_properties
#'
#' @aliases aa_properties properties kidera KIDERA aa_prop AA_PROP atchley ATCHLEY
NULL


#' Amino acid / codon table
#'
#' @concept data
#'
#' @docType data
#'
#' @name aa_table
#'
#' @aliases AA_TABLE AA_TABLE_REVERSED
AA_TABLE <- table(c(
  "TCA", "TCG", "TCC", "TCT", "TTT", "TTC", "TTA", "TTG", "TAT", "TAC", "TAA", "TAG", "TGT",
  "TGC", "TGA", "TGG", "CTA", "CTG", "CTC", "CTT", "CCA", "CCG", "CCC", "CCT", "CAT", "CAC",
  "CAA", "CAG", "CGA", "CGG", "CGC", "CGT", "ATT", "ATC", "ATA", "ATG", "ACA", "ACG", "ACC",
  "ACT", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTA", "GTG", "GTC", "GTT",
  "GCA", "GCG", "GCC", "GCT", "GAT", "GAC", "GAA", "GAG", "GGA", "GGG", "GGC", "GGT", "NNN"
))
AA_TABLE[c(
  "TCA", "TCG", "TCC", "TCT", "TTT", "TTC", "TTA", "TTG", "TAT", "TAC", "TAA", "TAG", "TGT",
  "TGC", "TGA", "TGG", "CTA", "CTG", "CTC", "CTT", "CCA", "CCG", "CCC", "CCT", "CAT", "CAC",
  "CAA", "CAG", "CGA", "CGG", "CGC", "CGT", "ATT", "ATC", "ATA", "ATG", "ACA", "ACG", "ACC",
  "ACT", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTA", "GTG", "GTC", "GTT",
  "GCA", "GCG", "GCC", "GCT", "GAT", "GAC", "GAA", "GAG", "GGA", "GGG", "GGC", "GGT", "NNN"
)] <- c(
  "S", "S", "S", "S", "F", "F", "L", "L", "Y", "Y", "*", "*", "C", "C", "*",
  "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R",
  "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S",
  "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E",
  "G", "G", "G", "G", "~"
)
AA_TABLE_REVERSED <- sapply(unique(AA_TABLE), function(aa) {
  names(AA_TABLE)[AA_TABLE == aa]
})
AA_TABLE_REVERSED <- AA_TABLE_REVERSED[order(names(AA_TABLE_REVERSED))]


#' Immune repertoire dataset
#'
#' @concept data
#'
#' @description A dataset with single chain data for testing and examplatory purposes.
#'
#' @format A list of two elements. First element ("data") is a list with data frames with clonotype tables.
#' Second element ("meta") is a metadata table.
#' \describe{
#'   \item{data}{List of immune repertoire data frames.}
#'   \item{meta}{Metadata}
#'   ...
#' }
"immdata"
