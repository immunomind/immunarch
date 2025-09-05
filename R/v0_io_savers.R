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
