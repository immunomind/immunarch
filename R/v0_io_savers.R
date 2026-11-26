save_immunarch <- function(.data, .path, .compress = TRUE) {
  if (.compress) {
    filepath <- gzfile(paste0(.path, ".tsv.gz"), open = "wb", compression = 9)
    on.exit(close(filepath), add = TRUE)
  } else {
    filepath <- paste0(.path, ".tsv")
  }
  readr::write_lines(paste0("# Exported from immunarch ", packageVersion("immunarch"), " https://immunarch.com"),
    file = filepath
  )
  readr::write_tsv(x = .data, file = filepath, append = TRUE, col_names = TRUE)
}

save_vdjtools <- function(.data, .path, .compress = TRUE) {
  if (.compress) {
    filepath <- gzfile(paste0(.path, ".tsv.gz"), open = "wb", compression = 9)
    on.exit(close(filepath), add = TRUE)
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

  matching_names <- match(names(.data), old)
  names(.data)[!is.na(matching_names)] <- new[matching_names[!is.na(matching_names)]]

  readr::write_tsv(x = .data, file = filepath)
}
