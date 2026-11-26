test_that("unstratified incidence counts each identity once across all repertoires", {
  idata <- make_track_idata()
  annotated <- suppressMessages(annotate_public(idata))
  before <- dplyr::collect(annotated)
  p <- vis(annotated, plot = "incidence")
  counts <- p$data |> dplyr::arrange(incidence)
  expect_s3_class(p, "ggplot")
  expect_equal(counts$incidence, c(1, 3))
  expect_equal(counts$receptors, c(3L, 4L))
  expect_equal(sum(counts$receptors), 7)
  expect_equal(nlevels(counts$stratum), 1L)
  expect_equal(p$guides$guides$fill, "none")
  expect_equal(p$theme$legend.position, "right")
  expect_equal(p$labels$title, "Immune receptor sharing across repertoires")
  expect_equal(dplyr::collect(annotated), before)
  expect_error(ggplot2::ggplotGrob(p), NA)
})

test_that("unstratified selected receptor plots retain only the measure legend", {
  idata <- make_track_idata()
  p <- vis(vdj_public(idata), plot = "public", reference = idata)
  expect_equal(nlevels(p$data$public_group), 1L)
  expect_equal(p$guides$guides$fill, "none")
  expect_setequal(p$scales$get_scales("shape")$get_labels(), c("Abundance", "Coverage"))
  expect_equal(p$theme$legend.position, "right")
})

test_that("stratified incidence stacks within-stratum identities without chain duplication", {
  idata <- make_public_idata()
  annotated <- annotate_public(idata)
  p <- vis(annotated, plot = "incidence")
  actual <- p$data |> dplyr::arrange(stratum, incidence)
  expect_equal(actual$incidence, c(1, 2, 1))
  expect_equal(actual$receptors, c(1L, 2L, 3L))
  # Five global identities produce six receptor-stratum occurrences: the
  # shared AAA identity is observed in both Case and Control.
  expect_equal(sum(actual$receptors), 6)
  expect_equal(p$labels$x, "Number of repertoires within stratum")
  expect_equal(p$labels$title, "Immune receptor sharing within strata")
  expect_equal(p$theme$legend.position, "right")
  expect_error(ggplot2::ggplotGrob(p), NA)
  expect_null(vis(annotated, plot = "incidence", title = NULL)$labels$title)
  expect_error(vis(idata, plot = "incidence"), "annotate_public")
})

test_that("long stratum names wrap only for display in both publicness plots", {
  idata <- make_public_idata()
  strata <- dplyr::collect(idata$strata)
  labels <- paste("A descriptive cohort with a deliberately long name for", strata$group)
  idata <- immundata::rename_strata(idata,
    stats::setNames(labels, as.character(strata$imd_strata_id)))
  plots <- list(
    vis(annotate_public(idata), plot = "incidence"),
    vis(vdj_public(idata), plot = "public", reference = idata))
  for (p in plots) {
    scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")
    expect_true(all(grepl("\n", scale$get_labels())))
    expect_setequal(gsub("\n", " ", scale$get_labels()), labels)
    expect_equal(p$theme$legend.position, "right")
    expect_error(ggplot2::ggplotGrob(p), NA)
  }
  expect_setequal(dplyr::collect(idata$strata)$strata_name, labels)
})

test_that("paired observations do not inflate annotated incidence histograms", {
  rows <- tibble::tibble(
    sample = rep(c("S1", "S1", "S2"), each = 2),
    cell = rep(c("a1", "a2", "b1"), each = 2),
    junction = rep(c("AAA", "CCC"), 3),
    v_call = rep(c("HV1", "LV1"), 3),
    locus = rep(c("IGH", "IGL"), 3), umi = 10L)
  idata <- make_track_idata(rows, paired = TRUE)
  p <- vis(suppressMessages(annotate_public(idata)), plot = "incidence")
  expect_equal(p$data$incidence, 2)
  expect_equal(p$data$receptors, 1L)
})
