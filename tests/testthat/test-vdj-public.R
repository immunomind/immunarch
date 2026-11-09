test_that("default publicness selects shared identities and ordinary observation rows", {
  idata <- make_public_idata()
  before <- dplyr::collect(idata) |> dplyr::arrange(row_id)
  selected <- vdj_public(idata)
  actual <- dplyr::collect(selected) |> dplyr::arrange(row_id)
  expect_s3_class(selected, "ImmunData")
  expect_s3_class(selected$annotations, "duckplyr_df")
  expect_equal(actual$row_id, c(1L, 2L, 4L, 5L, 6L))
  expect_equal(actual$imd_orig_count, before$imd_count[actual$row_id])
  expect_equal(actual$imd_orig_proportion, before$imd_proportion[actual$row_id])
  expect_equal(dplyr::collect(idata) |> dplyr::arrange(row_id), before)
  expect_false(any(grepl("^vdj_public_incidence", names(actual))))
  expect_null(attr(selected, "vdj_context", exact = TRUE))
  collected <- vdj_public(idata, output = "df") |> dplyr::arrange(row_id)
  expect_s3_class(collected, "data.frame")
  expect_s3_class(collected, "immunarch_res_vdj_public")
  expect_equal(tibble::as_tibble(collected), actual)
})


test_that("rules select identities by IDs and retain occurrences outside the population", {
  idata <- make_public_idata()
  reps <- idata$repertoires
  ids <- reps$imd_repertoire_id[reps$sample %in% c("S1", "S2")]
  actual <- vdj_public(idata, make_public_filter(repertoires = ids, min_incidence = 2),
                       output = "df")
  expect_setequal(actual$row_id, c(1L, 2L, 4L, 5L, 6L))
  expect_true(any(actual$sample == "S3"))
})


test_that("strata labels and IDs select the same population and intersect repertoire IDs", {
  idata <- make_public_idata()
  strata <- idata$strata
  idata <- immundata::rename_strata(
    idata, stats::setNames(paste0("label-", strata$group), as.character(strata$imd_strata_id))
  )
  case_id <- strata$imd_strata_id[strata$group == "Case"]
  by_label <- vdj_public(idata, make_public_filter(strata = "label-Case", min_incidence = 2),
                         output = "df")
  by_id <- vdj_public(idata, make_public_filter(strata = case_id, min_incidence = 2),
                      output = "df")
  expect_setequal(by_label$row_id, c(1L, 2L, 4L, 5L, 6L))
  expect_equal(dplyr::arrange(by_id, row_id), dplyr::arrange(by_label, row_id))
  sid <- idata$repertoires$imd_repertoire_id[idata$repertoires$sample %in% c("S2", "S3")]
  intersected <- vdj_public(idata, make_public_filter(repertoires = sid, strata = case_id),
                            output = "df")
  expect_setequal(intersected$row_id, c(1L, 2L, 4L, 5L, 6L))
})


test_that("incidence distinguishes distinct repertoires from distinct strata", {
  idata <- make_public_idata()
  repertoires <- vdj_public(idata, make_public_filter(min_incidence = 2), output = "df")
  strata <- vdj_public(idata, make_public_filter(incidence_by = "strata", min_incidence = 2),
                       output = "df")
  expect_setequal(repertoires$junction, c("AAA", "CCC"))
  expect_setequal(strata$row_id, c(1L, 4L, 6L))
})


test_that("all rules use original input and presence and absence are combined with AND", {
  idata <- make_public_idata()
  groups <- idata$strata
  case <- groups$imd_strata_id[groups$group == "Case"]
  control <- groups$imd_strata_id[groups$group == "Control"]
  present <- make_public_filter(strata = case, min_incidence = 2,
                                min_proportion = 1 / 7, max_proportion = 2 / 7)
  absent <- make_public_filter(strata = control, max_incidence = 0)
  forward <- vdj_public(idata, present, absent, output = "df")
  backward <- vdj_public(idata, absent, present, output = "df")
  expect_setequal(forward$row_id, c(2L, 5L))
  expect_equal(dplyr::arrange(forward, row_id), dplyr::arrange(backward, row_id))
})


test_that("inclusive abundance bounds qualify occurrences before incidence", {
  idata <- make_public_idata()
  capped <- vdj_public(idata, make_public_filter(max_count = 10, min_incidence = 2),
                       output = "df")
  expect_setequal(capped$row_id, c(1L, 2L, 4L, 5L, 6L))
  expect_true(any(capped$imd_orig_count == 20))
  bounded <- vdj_public(idata, make_public_filter(min_count = 5, max_count = 10,
                                                min_incidence = 2), output = "df")
  expect_setequal(bounded$row_id, c(1L, 4L, 6L))
  proportions <- vdj_public(idata, make_public_filter(min_proportion = 1 / 7,
                                                     max_proportion = 2 / 7,
                                                     min_incidence = 2), output = "df")
  expect_setequal(proportions$row_id, c(2L, 5L))
})


test_that("implicit minimum one and explicit maximum alone have different zero behavior", {
  idata <- make_public_idata()
  threshold <- vdj_public(idata, make_public_filter(min_count = 5), output = "df")
  expect_setequal(threshold$row_id, c(1L, 4L, 6L, 7L))
  below <- vdj_public(idata, make_public_filter(min_count = 5, max_incidence = 0),
                      output = "df")
  expect_setequal(below$row_id, c(2L, 3L, 5L, 8L))
  at_most_one <- vdj_public(idata, make_public_filter(min_count = 5, max_incidence = 1),
                            output = "df")
  expect_setequal(at_most_one$row_id, c(2L, 3L, 5L, 7L, 8L))
})


test_that("paired chains and repeated cells are deduplicated without splitting receptors", {
  rows <- tibble::tibble(
    sample = rep(c("S1", "S1", "S1", "S2", "S2", "S3"), each = 2),
    cell = rep(c("a1", "a2", "a3", "b1", "b2", "c1"), each = 2),
    junction = c("AAA", "CCC", "AAA", "CCC", "AAA", "GGG",
                 "AAA", "CCC", "TTT", "GGG", "AAA", "GGG"),
    v_call = c("HV1", "LV1", "HV1", "LV1", "HV1", "LV2",
               "HV1", "LV1", "HV2", "LV2", "HV1", "LV2"),
    locus = rep(c("IGH", "IGL"), 6), umi = 10L, row_id = seq_len(12L)
  )
  idata <- make_track_idata(rows, paired = TRUE)
  shared <- vdj_public(idata, output = "df")
  expect_setequal(shared$cell, c("a1", "a2", "a3", "b1", "c1"))
  expect_equal(nrow(shared), 10L)
  expect_setequal(shared$locus, c("IGH", "IGL"))
  bounded <- vdj_public(idata, make_public_filter(min_count = 2, max_count = 2),
                        output = "df")
  expect_setequal(bounded$cell, c("a1", "a2", "b1"))
  expect_true(all(bounded$imd_orig_count[bounded$sample == "S1"] == 2))
  empty <- vdj_public(idata, make_public_filter(min_incidence = 3), output = "df")
  expect_equal(nrow(empty), 0L)
})


test_that("repeated calls replace saved abundances and normal proportions are recalculated", {
  idata <- make_public_idata() |>
    immundata::mutate_immundata(imd_orig_count = -1, imd_orig_proportion = -1)
  first <- vdj_public(idata)
  second <- vdj_public(first, make_public_filter(min_count = 5), output = "df")
  expect_equal(second$imd_orig_count[second$sample == "S1"], 5)
  expect_equal(second$imd_orig_proportion[second$sample == "S1"], 5 / 6)
  expect_true(all(second$imd_proportion == 1))
  expect_true(all(dplyr::collect(idata)$imd_orig_count == -1))
})


test_that("empty selections return ordinary empty outputs", {
  idata <- make_public_idata()
  empty <- vdj_public(idata, make_public_filter(min_incidence = 5))
  expect_s3_class(empty, "ImmunData")
  expect_s3_class(empty$annotations, "duckplyr_df")
  expect_equal(nrow(dplyr::collect(empty)), 0L)
  collected <- vdj_public(idata, make_public_filter(min_incidence = 5), output = "df")
  expect_s3_class(collected, "data.frame")
  expect_equal(nrow(collected), 0L)
})


test_that("unknown selectors and empty comparison populations error", {
  idata <- make_public_idata()
  expect_error(vdj_public(idata, make_public_filter(repertoires = 999L)), "Unknown repertoire")
  expect_error(vdj_public(idata, make_public_filter(strata = "unknown")), "Unknown strata")
  case <- idata$strata$imd_strata_id[idata$strata$group == "Case"]
  sid <- idata$repertoires$imd_repertoire_id[idata$repertoires$sample == "S3"]
  expect_error(vdj_public(idata, make_public_filter(repertoires = sid, strata = case,
                                                  max_incidence = 0)), "[Ee]mpty")
})


test_that("repertoire rules work without strata but strata rules require them", {
  idata <- make_track_idata()
  actual <- vdj_public(idata, output = "df")
  expect_setequal(actual$row_id, setdiff(seq_len(15L), c(4L, 7L, 14L)))
  expect_error(vdj_public(idata, make_public_filter(incidence_by = "strata")),
               "Strata must be defined")
  expect_error(vdj_public(idata, make_public_filter(strata = "Case")),
               "Strata must be defined")
})
