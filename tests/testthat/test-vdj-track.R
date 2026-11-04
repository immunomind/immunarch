test_that("default tracking preserves every observation and input abundance", {
  idata <- make_track_idata()
  before <- dplyr::collect(idata) |> dplyr::arrange(row_id)
  before_reps <- idata$repertoires
  tracked <- vdj_track(idata)
  actual <- dplyr::collect(tracked) |> dplyr::arrange(row_id)

  expect_s3_class(tracked, "ImmunData")
  expect_s3_class(tracked$annotations, "duckplyr_df")
  expect_false(inherits(tracked, "immunarch_res"))
  expect_null(attr(tracked, "vdj_context", exact = TRUE))
  expect_equal(actual$row_id, before$row_id)
  expect_equal(actual$imd_orig_count, before$imd_count)
  expect_equal(actual$imd_orig_proportion, before$imd_proportion)
  expect_equal(actual$abundance, before$abundance)
  expect_equal(tracked$schema_receptor, idata$schema_receptor)
  expect_equal(dplyr::collect(idata) |> dplyr::arrange(row_id), before)
  expect_equal(idata$repertoires, before_reps)
})


test_that("repertoire IDs select identities that are followed throughout the target", {
  idata <- make_track_idata()
  reps <- idata$repertoires
  sid <- reps$imd_repertoire_id[reps$sample == "S1"]
  one <- vdj_track(idata, from = sid) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(one$row_id, c(1L, 2L, 3L, 4L, 5L, 6L, 8L, 9L, 10L, 12L))
  expect_setequal(one$sample, c("S1", "S2", "S3", "S4"))

  ids <- reps$imd_repertoire_id[reps$sample %in% c("S4", "S5")]
  multiple <- vdj_track(idata, from = ids) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_setequal(multiple$junction, c("CCC", "ACG", "GGA"))
  repeated <- vdj_track(idata, from = rep(ids, 2L)) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(repeated, multiple)
})


test_that("top-N shortcuts select per repertoire and follow the union everywhere", {
  idata <- make_track_idata()
  reps <- idata$repertoires
  sid <- reps$imd_repertoire_id[reps$sample == "S1"]
  from_one <- vdj_track(idata, from = sid, n = 2L) |> dplyr::collect()
  expect_setequal(from_one$junction, c("AAA", "CCC"))
  expect_true(all(from_one$v_call[from_one$junction == "AAA"] == "V1"))
  expect_true(any(from_one$sample == "S2" & from_one$junction == "AAA"))

  source <- immundata::filter_immundata(idata, imd_repertoire_id == !!sid) |>
    immundata::top_n(n = 3L)
  shortcut <- vdj_track(idata, from = sid, n = 3L) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  upstream <- vdj_track(idata, from = source) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(shortcut, upstream)

  each <- vdj_track(idata, n = 1L)
  identities <- each$receptors |> dplyr::collect()
  expect_equal(nrow(identities), 4L)
  expect_setequal(identities$junction, c("AAA", "CCC", "ACG", "GGA"))
  low_rank <- dplyr::collect(each) |>
    dplyr::filter(sample == "S2", junction == "AAA")
  expect_equal(low_rank$imd_orig_count, 1)
  expect_equal(low_rank$imd_orig_proportion, 1 / 15)

  fifth_id <- reps$imd_repertoire_id[reps$sample == "S5"]
  fifth <- vdj_track(idata, from = fifth_id, n = 3L) |> dplyr::collect()
  expect_setequal(fifth$junction, c("ACG", "GGA"))
  expect_true(any(fifth$sample == "S3"))
  oversized <- vdj_track(idata, n = 100L) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  all_receptors <- vdj_track(idata) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(oversized, all_receptors)
})


test_that("collected output has identical rows and original target abundances", {
  idata <- make_track_idata()
  query <- tibble::tibble(junction = "CCC")
  tracked <- vdj_track(idata, from = query)
  actual <- vdj_track(idata, from = query, output = "data.frame") |>
    dplyr::arrange(row_id)
  expect_s3_class(actual, "data.frame")
  expect_equal(actual, dplyr::collect(tracked) |> dplyr::arrange(row_id))
  expect_null(attr(actual, "vdj_context", exact = TRUE))
  expect_equal(actual$imd_orig_count, c(2, 9, 6))
  expect_equal(actual$imd_orig_proportion, c(2 / 12, 9 / 15, 6 / 7))
  expect_true(all(actual$imd_proportion == 1))
  expect_equal(actual$imd_count, actual$imd_orig_count)
})


test_that("ImmunData sources match features rather than independent numeric IDs", {
  idata <- make_track_idata()
  dictionary <- idata$receptors |> dplyr::collect()
  chosen <- dictionary[which.max(dictionary$imd_receptor_id), ]
  rows <- dplyr::collect(idata) |>
    dplyr::filter(junction == chosen$junction, v_call == chosen$v_call) |>
    dplyr::select(sample, junction, v_call, abundance, row_id)
  source <- make_track_idata(rows, repertoires = NULL)
  expect_false(identical(source$receptors |> dplyr::pull(imd_receptor_id),
                         chosen$imd_receptor_id))
  actual <- vdj_track(idata, from = source) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(actual$row_id, sort(rows$row_id))
  expect_true(all(actual$junction == chosen$junction))
  expect_true(all(actual$v_call == chosen$v_call))

  unmatched_rows <- rows[1, ]
  unmatched_rows$junction <- "ATAT"
  unmatched <- vdj_track(idata, from = make_track_idata(unmatched_rows))
  expect_equal(nrow(dplyr::collect(unmatched)), 0L)
  incompatible <- make_track_idata(features = "junction")
  expect_error(vdj_track(idata, from = incompatible), "same receptor features")
})


test_that("query tables match exact tuples while partial queries retain separate receptors", {
  idata <- make_track_idata()
  queries <- tibble::tibble(junction = c("AAA", "CCC"), v_call = c("V1", "V2"))
  actual <- vdj_track(idata, from = queries) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_false(any(actual$junction == "AAA" & actual$v_call == "V2"))
  expect_equal(nrow(actual), 6L)
  duplicated_query <- vdj_track(idata, from = queries[c(1, 2, 1), ]) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  lazy_query <- vdj_track(idata, from = duckplyr::as_duckdb_tibble(queries)) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(duplicated_query, actual)
  expect_equal(lazy_query, actual)

  partial <- vdj_track(idata, from = tibble::tibble(junction = "AAA"))
  expect_equal(nrow(partial$receptors |> dplyr::collect()), 2L)
  expect_setequal(dplyr::collect(partial)$v_call, c("V1", "V2"))
  hit <- vdj_track(idata, from = tibble::tibble(trait = "hit")) |>
    dplyr::collect() |> dplyr::arrange(row_id)
  expect_equal(hit$row_id, c(2L, 6L, 12L))
  expect_setequal(hit$trait, c("hit", "other"))
  literal <- vdj_track(idata, from = tibble::tibble(junction = "aaa"))
  expect_equal(nrow(dplyr::collect(literal)), 0L)

  receptor_id <- idata$receptors |>
    dplyr::filter(junction == "TTT") |> dplyr::pull(imd_receptor_id)
  by_id <- vdj_track(idata, from = tibble::tibble(imd_receptor_id = receptor_id))
  expect_equal(dplyr::collect(by_id)$row_id, 7L)
})


test_that("paired matching preserves whole receptors and keeps joint identities distinct", {
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
  before <- dplyr::collect(idata) |> dplyr::arrange(row_id)
  source <- make_track_idata(rows[rows$cell == "a3", ], paired = TRUE,
                             repertoires = NULL)
  paired <- vdj_track(idata, from = source) |> dplyr::collect()
  expect_setequal(paired$cell, c("a3", "c1"))
  expect_equal(nrow(paired), 4L)
  expect_setequal(paired$locus, c("IGH", "IGL"))

  heavy <- vdj_track(idata, from = tibble::tibble(locus = "IGH", junction = "AAA"))
  actual <- dplyr::collect(heavy)
  expect_equal(nrow(heavy$receptors |> dplyr::collect()), 2L)
  expect_equal(nrow(actual), 10L)
  expect_setequal(actual$cell, c("a1", "a2", "a3", "b1", "c1"))
  expect_setequal(actual$locus, c("IGH", "IGL"))
  expect_true(all(actual$imd_orig_count[actual$cell %in% c("a1", "a2")] == 2))
  expect_true(all(actual$imd_orig_proportion[actual$cell %in% c("a1", "a2")] == 2 / 3))

  crossed <- vdj_track(idata, from = tibble::tibble(junction = "AAA", v_call = "LV1"))
  expect_equal(nrow(dplyr::collect(crossed)), 0L)
  expect_equal(dplyr::collect(idata) |> dplyr::arrange(row_id), before)
})


test_that("original abundance columns are replaced relative to each call's input", {
  idata <- make_track_idata() |>
    immundata::mutate_immundata(imd_orig_count = -1, imd_orig_proportion = -1)
  first <- vdj_track(idata, from = tibble::tibble(junction = "AAA"))
  second <- vdj_track(first, from = tibble::tibble(v_call = "V1")) |>
    dplyr::collect()
  s1 <- second |> dplyr::filter(sample == "S1")
  expect_equal(s1$imd_orig_count, 8)
  expect_equal(s1$imd_orig_proportion, 8 / 9)
  expect_equal(s1$imd_proportion, 1)
  expect_true(all(dplyr::collect(idata)$imd_orig_count == -1))
})


test_that("empty selections are ordinary empty objects and collected data frames", {
  idata <- make_track_idata()
  empty_queries <- list(tibble::tibble(junction = "not-observed"),
                        tibble::tibble(junction = character()))
  for (query in empty_queries) {
    empty <- vdj_track(idata, from = query)
    expect_s3_class(empty, "ImmunData")
    expect_s3_class(empty$annotations, "duckplyr_df")
    expect_equal(nrow(dplyr::collect(empty)), 0L)
    collected <- vdj_track(idata, from = query, output = "data.frame")
    expect_s3_class(collected, "data.frame")
    expect_equal(nrow(collected), 0L)
    expect_true(all(c("imd_orig_count", "imd_orig_proportion") %in% names(collected)))
  }
  zero <- vdj_track(idata, n = 0L)
  expect_equal(nrow(dplyr::collect(zero)), 0L)
  zero_from_one <- vdj_track(idata, from = idata$repertoires$imd_repertoire_id[1], n = 0L)
  expect_equal(nrow(dplyr::collect(zero_from_one)), 0L)
})


test_that("tracking preserves ordinary strata and target annotations", {
  idata <- make_track_idata(repertoires = c("sample", "group")) |>
    immundata::agg_strata("group")
  labels <- idata$strata
  idata <- immundata::rename_strata(
    idata, stats::setNames(paste0("label-", labels$group), as.character(labels$imd_strata_id))
  )
  actual <- vdj_track(idata, from = tibble::tibble(junction = "GGA"))
  expect_equal(actual$schema_strata, "group")
  expect_equal(actual$strata$strata_name, "label-Case")
  expect_equal(dplyr::collect(actual)$row_id, 14L)
  expect_equal(dplyr::collect(actual)$group, "Case")
})


test_that("tracking validates source forms, queries, ranking arguments and output", {
  idata <- make_track_idata()
  expect_error(vdj_track(idata, from = "S1"), "must be NULL")
  expect_error(vdj_track(idata, from = list(1)), "must be NULL")
  expect_error(vdj_track(idata, from = numeric()), "length")
  expect_error(vdj_track(idata, from = NA_real_))
  expect_error(vdj_track(idata, from = 1.5))
  expect_error(vdj_track(idata, from = 999L), "Unknown repertoire IDs")
  expect_error(vdj_track(idata, from = c(1L, 999L)), "Unknown repertoire IDs")
  expect_error(vdj_track(idata, from = tibble::tibble()), "at least one column")
  expect_error(vdj_track(idata, from = tibble::tibble(unknown = "AAA")), "Unknown query columns")
  expect_error(vdj_track(idata, from = tibble::tibble(junction = NA_character_)), "Missing query values")
  expect_error(vdj_track(idata, from = tibble::tibble(junction = c("AAA", NA_character_))), "Missing query values")
  expect_error(vdj_track(idata, from = idata, n = 1L), "applies only to repertoire sources")
  expect_error(vdj_track(idata, from = tibble::tibble(junction = "AAA"), n = 1L), "applies only")
  expect_error(vdj_track(idata, n = -1L))
  expect_error(vdj_track(idata, n = 0.5))
  expect_error(vdj_track(idata, n = NA_real_))
  expect_error(vdj_track(idata, output = "matrix"), "arg")
  expect_error(vdj_track(make_track_idata(repertoires = NULL)), "Repertoire aggregation")
  expect_named(formals(vdj_track), c("idata", "from", "n", "output"))
})
