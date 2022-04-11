#' Specification of the data format used by immunarch dataframes
#'
#' @name immunr_data_format
#'
#' @concept data
#'
#' @aliases immunr_data_format immunarch_data_format
#'
#' @docType data
#'
#' @description
#'
#' - "Clones" - number of barcodes (events, UMIs) or reads;
#'
#' - "Proportion" - proportion of barcodes (events, UMIs) or reads;
#'
#' - "CDR3.nt" - CDR3 nucleotide sequence;
#'
#' - "CDR3.aa" - CDR3 amino acid sequence;
#'
#' - "V.name" - names of aligned Variable gene segments;
#'
#' - "D.name" - names of aligned Diversity gene segments or NA;
#'
#' - "J.name" - names of aligned Joining gene segments;
#'
#' - "V.end" - last positions of aligned V gene segments (1-based);
#'
#' - "D.start" - positions of D'5 end of aligned D gene segments (1-based);
#'
#' - "D.end" - positions of D'3 end of aligned D gene segments (1-based);
#'
#' - "J.start" - first positions of aligned J gene segments (1-based);
#'
#' - "VJ.ins" - number of inserted nucleotides (N-nucleotides) at V-J junction (-1 for receptors with VDJ recombination);
#'
#' - "VD.ins" - number of inserted nucleotides (N-nucleotides) at V-D junction (-1 for receptors with VJ recombination);
#'
#' - "DJ.ins" - number of inserted nucleotides (N-nucleotides) at D-J junction (-1 for receptors with VJ recombination);
#'
#' - "Sequence" - full nucleotide sequence.
#'
NULL

# Immunarch basic column names
IMMCOL <- new.env()

IMMCOL$count <- "Clones"
IMMCOL$prop <- "Proportion"
IMMCOL$cdr3nt <- "CDR3.nt"
IMMCOL$cdr3aa <- "CDR3.aa"
IMMCOL$v <- "V.name"
IMMCOL$d <- "D.name"
IMMCOL$j <- "J.name"
IMMCOL$ve <- "V.end"
IMMCOL$ds <- "D.start"
IMMCOL$de <- "D.end"
IMMCOL$js <- "J.start"
IMMCOL$vnj <- "VJ.ins"
IMMCOL$vnd <- "VD.ins"
IMMCOL$dnj <- "DJ.ins"
IMMCOL$seq <- "Sequence"
IMMCOL$order <- c(
  IMMCOL$count, IMMCOL$prop, IMMCOL$cdr3nt, IMMCOL$cdr3aa,
  IMMCOL$v, IMMCOL$d, IMMCOL$j,
  IMMCOL$ve, IMMCOL$ds, IMMCOL$de, IMMCOL$js,
  IMMCOL$vnj, IMMCOL$vnd, IMMCOL$dnj, IMMCOL$seq
)
IMMCOL$type <- c(
  "numeric", "numeric", "character", "character",
  "character", "character", "character",
  "integer", "integer", "integer", "integer",
  "integer", "integer", "integer", "character"
)

# TODO: move V/D/J coordinates here
# Immunarch extended column names
IMMCOL_EXT <- new.env()

IMMCOL_EXT$bestv <- "Best.V"
IMMCOL_EXT$bestj <- "Best.J"
IMMCOL_EXT$cdr3s <- "CDR3.start"
IMMCOL_EXT$cdr3e <- "CDR3.end"
IMMCOL_EXT$c <- "C.name"
IMMCOL_EXT$cs <- "C.start"
IMMCOL_EXT$ce <- "C.end"
IMMCOL_EXT$cdr1nt <- "CDR1.nt"
IMMCOL_EXT$cdr1aa <- "CDR1.aa"
IMMCOL_EXT$cdr2nt <- "CDR2.nt"
IMMCOL_EXT$cdr2aa <- "CDR2.aa"
IMMCOL_EXT$fr1nt <- "FR1.nt"
IMMCOL_EXT$fr1aa <- "FR1.aa"
IMMCOL_EXT$fr2nt <- "FR2.nt"
IMMCOL_EXT$fr2aa <- "FR2.aa"
IMMCOL_EXT$fr3nt <- "FR3.nt"
IMMCOL_EXT$fr3aa <- "FR3.aa"
IMMCOL_EXT$fr4nt <- "FR4.nt"
IMMCOL_EXT$fr4aa <- "FR4.aa"
IMMCOL_EXT$v3del <- "V3.Deletions"
IMMCOL_EXT$j3del <- "J3.Deletions"
IMMCOL_EXT$order <- c(
  IMMCOL_EXT$bestv, IMMCOL_EXT$bestj, IMMCOL_EXT$cdr3s, IMMCOL_EXT$cdr3e,
  IMMCOL_EXT$c, IMMCOL_EXT$cs, IMMCOL_EXT$ce,
  IMMCOL_EXT$cdr1nt, IMMCOL_EXT$cdr1aa, IMMCOL_EXT$cdr2nt, IMMCOL_EXT$cdr2aa,
  IMMCOL_EXT$fr1nt, IMMCOL_EXT$fr1aa, IMMCOL_EXT$fr2nt, IMMCOL_EXT$fr2aa,
  IMMCOL_EXT$fr3nt, IMMCOL_EXT$fr3aa, IMMCOL_EXT$fr4nt, IMMCOL_EXT$fr4aa,
  IMMCOL_EXT$v3del, IMMCOL_EXT$j3del
)
IMMCOL_EXT$type <- c(
  "character", "character", "integer", "integer",
  "character", "integer", "integer",
  "character", "character", "character", "character",
  "character", "character", "character", "character",
  "character", "character", "character", "character",
  "integer", "integer"
)

# Additional information
IMMCOL_ADD <- new.env()
# Separator for values inside columns, e.g., V assignments.
IMMCOL_ADD$valsep <- ","
# Separator for paired chain data
IMMCOL_ADD$scsep <- ";"

IMMUNR_ERROR_NOT_IMPL <- "ERROR: not implemented yet. Please contact us via support@immunomind.io if you need it in your research."
