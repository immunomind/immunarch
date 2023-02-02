#' Specification of the data format for Single Cell
#'
#' @name immunr_sc_data_format
#'
#' @concept data
#'
#' @aliases immunr_sc_data_format immunarch_sc_data_format
#'
#' @docType data
#'
#' @description
#'
#' "Reads" - number of reads aligned to this contig;
#'
#' "UMIs" - number of distinct UMIs aligned to this contig;
#'
#' "FR1.nt" - predicted FR1 nucleotide sequence;
#'
#' "FR1.aa" - predicted FR1 amino acid sequence;
#'
#' "CDR1.nt" - predicted CDR1 nucleotide sequence;
#'
#' "CDR1.aa" - predicted CDR1 amino acid sequence;
#'
#' "FR2.nt" - predicted FR2 nucleotide sequence;
#'
#' "FR2.aa" - predicted FR2 amino acid sequence;
#'
#' "CDR2.nt" - predicted CDR2 nucleotide sequence;
#'
#' "CDR2.aa" - predicted CDR2 amino acid sequence;
#'
#' "FR3.nt" - predicted FR3 nucleotide sequence;
#'
#' "FR3.aa" - predicted FR3 amino acid sequence;
#'
#' "CDR3.nt" - predicted CDR3 nucleotide sequence;
#'
#' "CDR3.aa" - predicted CDR3 amino acid sequence;
#'
#' "FR4.nt" - predicted FR4 nucleotide sequence;
#'
#' "FR4.aa" - predicted FR4 amino acid sequence;
#'
#' "V.name" - name of aligned Variable gene segment;
#'
#' "D.name" - name of aligned Diversity gene segment;
#'
#' "J.name" - name of aligned Joining gene segment;
#'
#' "C.name" - name of aligned Constant gene segment;
#'
#' "Full.length" - TRUE or FALSE statement about whether the contig is full length;
#'
#' "Is.cell" - TRUE or FALSE value indicating whether the barcode was called as a cell;
#'
#' "Contig.ID" - Unique identifier for this contig;
#'
#' "Productive" - TRUE or FALSE value indicating if the contig was declared as productive;
#'
#' "Raw.consensus.ID" - ID of the clonotype to which this cell barcode was assigned;
#'
#' "Raw.clonotype.ID" - ID of the consensus sequence to which this contig was assigned.
#'
NULL

# column names and types for Immunarch Single Cell format
IMMCOL_SC <- list(
  Cell.ID = "character",
  Reads = "integer",
  UMIs = "integer",
  FR1.nt = "character",
  FR1.aa = "character",
  CDR1.nt = "character",
  CDR1.aa = "character",
  FR2.nt = "character",
  FR2.aa = "character",
  CDR2.nt = "character",
  CDR2.aa = "character",
  FR3.nt = "character",
  FR3.aa = "character",
  CDR3.nt = "character",
  CDR3.aa = "character",
  FR4.nt = "character",
  FR4.aa = "character",
  V.name = "character",
  D.name = "character",
  J.name = "character",
  C.name = "character",
  Full.length = "logical",
  Is.cell = "logical",
  Contig.ID = "character",
  Productive = "logical",
  Raw.consensus.ID = "character",
  Raw.clonotype.ID = "character"
)
