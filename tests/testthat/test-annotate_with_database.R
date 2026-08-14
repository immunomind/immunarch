test_that("annotate_with_database joins duckplyr database and applies filters", {

  idata <- make_test_idata()
  ann <- dplyr::collect(idata$annotations)

  keys <- unique(ann$cdr3_aa)[1:6]
  db <- tibble::tibble(
    cdr3_db = keys,
    species = c("A", "B", "C", "D", "E", "F"),
    gene = c("TRB", "TRA", "TRB", "TRA", "TRB", "TRA")
  )
  db <- duckplyr::as_duckdb_tibble(db, prudence = "lavish")

  out <- immunarch:::annotate_with_database(
    idata,
    db = db,
    by = c(cdr3_aa = "cdr3_db"),
    label_col = "species",
    gene == "TRB"
  )

  out_ann <- dplyr::collect(out$annotations)
  expect_true("species" %in% names(out_ann))
  expect_setequal(unique(stats::na.omit(out_ann$species)), c("A", "C", "E"))
})

test_that("annotate_with_database reads database from CSV file", {

  idata <- make_test_idata()
  ann <- dplyr::collect(idata$annotations)

  keys <- unique(ann$cdr3_aa)[1:4]
  db <- tibble::tibble(
    cdr3_db = keys,
    species = c("S1", "S2", "S3", "S4")
  )
  db <- duckplyr::as_duckdb_tibble(db, prudence = "lavish")

  db_path <- tempfile(fileext = ".csv")
  on.exit(unlink(db_path), add = TRUE)
  utils::write.csv(dplyr::collect(db), db_path, row.names = FALSE)

  out <- immunarch:::annotate_with_database(
    idata,
    db = db_path,
    by = c(cdr3_aa = "cdr3_db"),
    label_col = "species"
  )

  out_ann <- dplyr::collect(out$annotations)
  expect_true("species" %in% names(out_ann))
  expect_gt(sum(!is.na(out_ann$species)), 0)
})

test_that("annotate_with_database errors when join columns are missing", {

  idata <- make_test_idata()
  ann <- dplyr::collect(idata$annotations)
  db <- tibble::tibble(
    cdr3_db = unique(ann$cdr3_aa)[1:3],
    species = c("A", "B", "C")
  )
  db <- duckplyr::as_duckdb_tibble(db, prudence = "lavish")

  expect_error(
    immunarch:::annotate_with_database(
      idata,
      db = db,
      by = c(not_a_column = "cdr3_db"),
      label_col = "species"
    ),
    "are not found in ImmunData"
  )

  expect_error(
    immunarch:::annotate_with_database(
      idata,
      db = db,
      by = c(cdr3_aa = "missing_in_db"),
      label_col = "species"
    ),
    "are not found in the database"
  )
})

test_that("annotate_with_database errors when label_col is missing", {

  idata <- make_test_idata()
  ann <- dplyr::collect(idata$annotations)
  db <- tibble::tibble(
    cdr3_db = unique(ann$cdr3_aa)[1:3],
    other_label = c("x", "y", "z")
  )
  db <- duckplyr::as_duckdb_tibble(db, prudence = "lavish")

  expect_error(
    immunarch:::annotate_with_database(
      idata,
      db = db,
      by = c(cdr3_aa = "cdr3_db"),
      label_col = "species"
    ),
    "is not found in the database"
  )
})

test_that("annotate_with_database errors for missing database file", {

  idata <- make_test_idata()
  missing_path <- file.path(tempdir(), "definitely_missing_db.csv")

  expect_error(
    immunarch:::annotate_with_database(
      idata,
      db = missing_path,
      by = c(cdr3_aa = "cdr3_db"),
      label_col = "species"
    ),
    "does not exist"
  )
})
