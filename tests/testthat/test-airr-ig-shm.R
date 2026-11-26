airr_shm_rows <- function() {
  values <- c(1, 3, 10, NA, 5, NA, NA, 9)
  data.frame(
    sequence_id = paste0("s", seq_along(values)),
    sample = c(rep("S1", 7), "S2"), locus = "IGH",
    clone_id = c("A", "A", "B", "C", NA, "D", "D", "A"),
    mu_count = values, mu_count_r = values,
    mu_count_s = ifelse(is.na(values), NA, 0),
    mu_freq = values / 100, mu_freq_r = values / 100,
    mu_freq_s = ifelse(is.na(values), NA, 0)
  )
}


make_airr_shm_idata <- function(rows = airr_shm_rows(), repertoire_schema = "sample") {
  input <- tempfile(fileext = ".tsv")
  readr::write_tsv(rows, input, na = "")
  immundata::read_repertoires(
    input, schema = "sequence_id", repertoire_schema = repertoire_schema,
    preprocess = NULL, postprocess = NULL, rename_columns = NULL,
    output_folder = tempfile("airr-shm-"), verbose = FALSE
  )
}


test_that("SHM repertoire summaries distinguish sequence and clone weighting", {
  idata <- make_airr_shm_idata()
  expect_s3_class(idata$annotations, "duckplyr_df")
  result <- airr_ig_shm(idata)
  expect_s3_class(result, "immunarch_res_airr_ig_shm")
  expect_equal(nrow(result), 12L)
  actual <- result |>
    dplyr::filter(sample == "S1", mutation == "total", measure == "freq")
  expect_equal(actual$n_sequences, 7)
  expect_equal(actual$n_valid, 4)
  expect_equal(actual$n_missing, 3)
  expect_equal(actual$n_valid_without_clone, 1)
  expect_equal(actual$n_clones, 4)
  expect_equal(actual$n_valid_clones, 2)
  expect_equal(actual$n_valid_with_clone, 3)
  expect_equal(actual$sequence_mean, .0475)
  expect_equal(actual$sequence_median, .04)
  expect_equal(actual$sequence_q25, .025)
  expect_equal(actual$sequence_q75, .0625)
  expect_equal(actual$clone_mean_weighted_mean, .14 / 3)
  expect_equal(actual$clone_mean_mean, .06)
  expect_equal(actual$clone_mean_median, .06)
  expect_equal(actual$clone_mean_q25, .04)
  expect_equal(actual$clone_mean_q75, .08)
  expect_equal(actual$clone_median_mean, .06)
  expect_equal(actual$clone_median_q25, .04)
  expect_equal(actual$clone_median_q75, .08)
  # Reusing clone A in a different repertoire must not combine the clones.
  second <- result |>
    dplyr::filter(sample == "S2", mutation == "total", measure == "freq")
  expect_equal(second$clone_mean_mean, .09)
  expect_equal(second$n_clones, 1)
  zeros <- result |>
    dplyr::filter(sample == "S1", mutation == "silent", measure == "freq")
  expect_equal(zeros$n_zero, 4)
  expect_equal(zeros$zero_prop, 1)
  expect_equal(zeros$n_zero_clone_means, 2)
  expect_equal(zeros$n_zero_clone_medians, 2)
  expect_s3_class(idata$annotations, "duckplyr_df")
})


test_that("missing SHM values, clone IDs and grouping keys keep their meaning", {
  rows <- airr_shm_rows()
  rows$locus[rows$sample == "S2"] <- NA_character_
  rows$clone_id[rows$sample == "S2"] <- NA_character_
  rows[rows$sample == "S2", grepl("^mu_", names(rows))] <- NA_real_
  result <- airr_ig_shm(make_airr_shm_idata(rows)) |>
    dplyr::filter(sample == "S2")
  expect_equal(nrow(result), 6L)
  expect_true(all(is.na(result$locus)))
  expect_true(all(result$n_sequences == 1))
  expect_true(all(result$n_valid == 0))
  expect_true(all(result$n_missing == 1))
  expect_true(all(result$n_clones == 0))
  expect_true(all(result$n_valid_clones == 0))
  expect_true(all(result$n_valid_with_clone == 0))
  expect_true(all(result$n_zero == 0))
  expect_true(all(result$n_zero_clone_means == 0))
  expect_true(all(is.na(result$zero_prop)))
  expect_true(all(is.na(result$sequence_mean)))
  expect_true(all(is.na(result$clone_mean_mean)))
  expect_true(all(is.na(result$clone_mean_weighted_mean)))
})


test_that("clone means and clone medians preserve different distributions", {
  rows <- airr_shm_rows()[1:4, ]
  rows$clone_id <- c("A", "A", "A", "B")
  rows$mu_count <- rows$mu_count_r <- c(0, 0, 9, 6)
  rows$mu_count_s <- 0
  rows$mu_freq <- rows$mu_freq_r <- rows$mu_count / 100
  rows$mu_freq_s <- 0
  actual <- airr_ig_shm(make_airr_shm_idata(rows)) |>
    dplyr::filter(mutation == "total", measure == "count")
  expect_equal(actual$sequence_mean, 3.75)
  expect_equal(actual$clone_mean_weighted_mean, 3.75)
  expect_equal(actual$clone_mean_mean, 4.5)
  expect_equal(actual$clone_median_mean, 3)
  expect_equal(actual$clone_median_q25, 1.5)
  expect_equal(actual$clone_median_q75, 4.5)
  expect_equal(actual$n_zero, 2)
  expect_equal(actual$zero_prop, .5)
  expect_equal(actual$n_zero_clone_means, 0)
  expect_equal(actual$n_zero_clone_medians, 1)
})


test_that("SHM grouping, metadata and long-only output follow AIRR conventions", {
  rows <- airr_shm_rows()
  rows$locus[2] <- "IGK"
  idata <- make_airr_shm_idata(rows)
  separated <- airr_ig_shm(idata, autojoin = FALSE)
  expect_equal(nrow(separated), 18L)
  expect_false("sample" %in% names(separated))
  pooled <- airr_ig_shm(idata, by = NULL)
  expect_equal(nrow(pooled), 12L)
  expect_false("locus" %in% names(pooled))
  # A repertoire metadata column used for grouping must not get .x/.y suffixes.
  grouped <- airr_ig_shm(idata, by = "sample")
  expect_true("sample" %in% names(grouped))
  expect_false(any(grepl("\\.[xy]$", names(grouped))))
  grouped <- grouped |> dplyr::arrange(sample, region, mutation, measure)
  pooled <- pooled |> dplyr::arrange(sample, region, mutation, measure)
  expect_equal(grouped$sequence_mean, pooled$sequence_mean)
  no_locus <- make_airr_shm_idata(rows[, names(rows) != "locus"])
  expect_false("locus" %in% colnames(no_locus))
  expect_equal(nrow(airr_ig_shm(no_locus)), 12L)
  expect_error(airr_ig_shm(no_locus, by = "locus"), "Grouping column")
  long <- airr_ig_shm(idata, format = "wide", autojoin = FALSE)
  expect_true(all(c("region", "mutation", "measure") %in% names(long)))
  expect_equal(nrow(long), nrow(separated))
})


test_that("SHM summaries validate inputs and quote custom column names", {
  rows <- airr_shm_rows()
  names(rows)[names(rows) == "clone_id"] <- "clone ID"
  rows[["group label"]] <- rows$sample
  idata <- make_airr_shm_idata(rows)
  expect_equal(nrow(airr_ig_shm(idata, clone_col = "clone ID", by = "group label")), 12L)
  expect_error(airr_ig_shm(idata), "not in the input ImmunData")
  expect_error(airr_ig_shm(idata, clone_col = "clone ID", by = "clone ID"), "must not include")
  expect_error(airr_ig_shm(idata, clone_col = "clone ID", by = TRUE), "character")
  expect_error(airr_ig_shm(make_airr_shm_idata(repertoire_schema = NULL)), "Repertoire aggregation")
  expect_error(airr_ig_shm(data.frame()), "R6")
  rows$mu_freq <- NULL
  expect_error(airr_ig_shm(make_airr_shm_idata(rows), clone_col = "clone ID"), "Missing SHM column")
})


test_that("regional totals use existing R/S values without combining frequencies", {
  rows <- airr_shm_rows()
  rows$mu_count_fwr1_r <- 1
  rows$mu_count_fwr1_s <- 2
  rows$mu_freq_fwr1_r <- .1
  rows$mu_freq_fwr1_s <- .2
  rows$mu_count_v <- 8
  rows$mu_count_v_r <- 5
  rows$mu_count_v_s <- 3
  rows$mu_freq_v <- .04
  rows$mu_freq_v_r <- .025
  rows$mu_freq_v_s <- .015
  result <- airr_ig_shm(make_airr_shm_idata(rows)) |>
    dplyr::filter(sample == "S1", mutation == "total")
  expect_setequal(result$region, c("overall", "v", "fwr1"))
  expect_equal(result$sequence_mean[result$region == "fwr1" & result$measure == "count"], 3)
  expect_equal(result$sequence_mean[result$region == "fwr1" & result$measure == "freq"], .3)
  expect_equal(result$sequence_mean[result$region == "v" & result$measure == "freq"], .04)
})
