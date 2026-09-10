write_imgt_fasta <- function(lines) {
  path <- tempfile(fileext = ".fasta")
  writeLines(lines, path)
  path
}


testthat::test_that("IMGT reference reader keeps nucleotide IG V/J records and their original metadata", {
  path <- write_imgt_fasta(c(
    ">M99641|IGHV1-18*01|Homo_sapiens|F|V-REGION|188..483|296 nt|1| | | | |",
    "caggtt...", "acgt",
    ">J00256|IGHJ1*01|Homo_sapiens|F|J-REGION|723..774|52 nt|1| | | | |",
    "gctgaata",
    ">X97051|IGHD1-1*01|Homo_sapiens|F|D-REGION|33714..33730|17 nt|1| | | | |",
    "ggtacaactggaacgac",
    ">AE000658|TRAV1-1*01|Homo_sapiens|F|V-REGION|128090..128364|275 nt|1| | | | |",
    "ggacaaag",
    ">M99641|IGHV1-18*01|Homo_sapiens|F|V-REGION|188..483|296 aa|1| | | | |",
    "CAGGTT"
  ))

  actual <- read_germline_references(path)

  testthat::expect_s3_class(actual, "tbl_df")
  testthat::expect_named(actual, c("species", "locus", "segment", "allele_id", "sequence"))
  testthat::expect_equal(
    actual,
    tibble::tribble(
      ~species, ~locus, ~segment, ~allele_id, ~sequence,
      "Homo_sapiens", "IGH", "V", "IGHV1-18*01", "CAGGTT...ACGT",
      "Homo_sapiens", "IGH", "J", "IGHJ1*01", "GCTGAATA"
    )
  )
})


testthat::test_that("IMGT reference reader rejects duplicate relevant records", {
  path <- write_imgt_fasta(c(
    ">A1|IGHV1*01|Homo_sapiens|F|V-REGION|1..4|4 nt|1|", "ACGT",
    ">A2|IGHV1*01|Homo_sapiens|F|V-REGION|1..4|4 nt|1|", "ACGT"
  ))

  testthat::expect_error(read_germline_references(path), "(?i)duplicate")
})


testthat::test_that("IMGT reference reader rejects malformed relevant records and files without IG V/J", {
  malformed <- write_imgt_fasta(c(
    ">A1|IGHV1*01|Homo_sapiens|F|V-REGION", "ACGT"
  ))
  no_usable <- write_imgt_fasta(c(
    ">A1|IGHD1*01|Homo_sapiens|F|D-REGION|1..4|4 nt|1|", "ACGT",
    ">A2|TRAV1*01|Homo_sapiens|F|V-REGION|1..4|4 nt|1|", "ACGT"
  ))

  testthat::expect_error(read_germline_references(malformed), "(?i)malformed")
  testthat::expect_error(read_germline_references(no_usable), "(?i)usable")
})
