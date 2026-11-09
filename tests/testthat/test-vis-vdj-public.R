test_that("public plots preserve the original repertoire universe and abundance baseline", {
  idata <- make_public_idata()
  selected <- vdj_public(idata)
  before <- dplyr::collect(selected) |> dplyr::arrange(row_id)
  p <- vis(selected, plot = "public", reference = idata, xval = "sample")
  actual <- p$data |> dplyr::arrange(sample)
  expect_s3_class(p, "ggplot")
  expect_equal(as.character(actual$sample), paste0("S", 1:4))
  expect_equal(actual$n_selected, rep(2, 4))
  expect_equal(actual$n_detected, c(2, 2, 1, 0))
  expect_equal(actual$coverage, c(1, 1, .5, 0))
  expect_equal(actual$abundance, c(6 / 7, 1, 20 / 28, 0))
  expect_equal(actual$coverage_label, c("2/2", "2/2", "1/2", "0/2"))
  expect_equal(actual$abundance_label[4], "0.00%")
  expect_equal(actual$imd_repertoire_id,
    idata$repertoires$imd_repertoire_id[match(as.character(actual$sample), idata$repertoires$sample)])
  expect_equal(dplyr::collect(selected) |> dplyr::arrange(row_id), before)
  collected <- vdj_public(idata, output = "df")
  expect_equal(vis(collected, reference = idata, xval = "sample")$data, p$data)
  expect_error(ggplot2::ggplotGrob(p), NA)
})

test_that("public plotting matches rebuilt repertoire IDs by their original metadata", {
  idata <- make_public_idata()
  strata <- idata$strata
  selected <- vdj_public(idata,
    make_public_filter(strata = strata$imd_strata_id[strata$group == "Case"], min_incidence = 2),
    make_public_filter(strata = strata$imd_strata_id[strata$group == "Control"], max_incidence = 0))
  p <- vis(selected, plot = "public", reference = idata, xval = "sample", fill = "group")
  actual <- p$data |> dplyr::arrange(sample)
  expect_equal(actual$n_detected, c(1, 1, 0, 0))
  expect_equal(actual$coverage, c(1, 1, 0, 0))
  expect_equal(actual$abundance, c(1 / 7, 4 / 14, 0, 0))
  expect_setequal(actual$group, c("Case", "Control"))
  expect_error(vis(selected, plot = "public", reference = idata, xval = "group"),
               "one repertoire")
})

test_that("paired chains and repeated cells do not inflate public coverage or abundance", {
  rows <- tibble::tibble(
    sample = rep(c("S1", "S1", "S1", "S2", "S2", "S3"), each = 2),
    cell = rep(c("a1", "a2", "a3", "b1", "b2", "c1"), each = 2),
    junction = c("AAA", "CCC", "AAA", "CCC", "AAA", "GGG",
                 "AAA", "CCC", "TTT", "GGG", "AAA", "GGG"),
    v_call = c("HV1", "LV1", "HV1", "LV1", "HV1", "LV2",
               "HV1", "LV1", "HV2", "LV2", "HV1", "LV2"),
    locus = rep(c("IGH", "IGL"), 6), umi = 10L)
  idata <- make_track_idata(rows, paired = TRUE)
  p <- vis(vdj_public(idata, output = "df"), reference = idata, xval = "sample")
  actual <- p$data |> dplyr::arrange(sample)
  expect_equal(actual$n_selected, rep(2, 3))
  expect_equal(actual$n_detected, c(2, 1, 1))
  expect_equal(actual$coverage, c(1, .5, .5))
  expect_equal(actual$abundance, c(1, .5, 1))
})

test_that("ten strata receive defined colours", {
  rows <- tibble::tibble(
    sample = rep(paste0("S", 1:10), each = 2),
    group = rep(paste0("G", 1:10), each = 2),
    junction = as.vector(rbind(rep("AAA", 10), paste0("private", 1:10))),
    v_call = "V1", abundance = 1L)
  idata <- make_track_idata(rows, repertoires = c("sample", "group")) |>
    immundata::agg_strata("group")
  p <- vis(vdj_public(idata, output = "df"), reference = idata, xval = "sample")
  expect_equal(nlevels(p$data$public_group), 10L)
  expect_warning(built <- ggplot2::ggplot_build(p), NA)
  expect_false(anyNA(built$data[[2]]$fill))
  expect_equal(length(unique(built$data[[2]]$fill)), 10L)
  expect_error(ggplot2::ggplotGrob(p), NA)
})

test_that("public selection plots require a reference and saved original proportions", {
  idata <- make_track_idata()
  selected <- vdj_public(idata)
  expect_error(vis(selected, plot = "public"), "reference")
  expect_error(vis(idata, plot = "public", reference = idata), "imd_orig_proportion")
  empty <- vdj_public(idata, make_public_filter(min_incidence = 99), output = "df")
  expect_error(vis(empty, reference = idata), "No selected receptors")
  # The existing default pairwise ImmunData plot keeps its meaning.
  p <- vis(idata, density = FALSE, lm = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("x", "y") %in% names(p$data)))
})
