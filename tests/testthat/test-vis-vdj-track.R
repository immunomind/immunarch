test_that("tracking plots dispatch on summaries and preserve their rows and abundances", {
  idata <- make_track_idata()
  tracked <- vdj_track(idata, from = tibble::tibble(junction = c("AAA", "CCC")),
                       output = "df")
  original <- tracked
  p <- vis(tracked, xval = "sample", plot = "line")
  expect_s3_class(tracked, "immunarch_res_vdj_track")
  expect_equal(p$data$imd_orig_proportion, tracked$imd_orig_proportion)
  expect_equal(as.numeric(as.character(p$data$imd_receptor_id)), tracked$imd_receptor_id)
  expect_equal(p$data$imd_repertoire_id, tracked$imd_repertoire_id)
  expect_true(any(p$data$imd_orig_proportion == 0))
  expect_equal(nrow(p$data), length(unique(tracked$imd_receptor_id)) *
                 nrow(idata$repertoires))
  expect_equal(tracked, original)
  expect_null(attr(tracked, "repertoires", exact = TRUE))
  expect_error(vis(dplyr::bind_rows(tracked, tracked)), "one row per receptor")
  incomplete <- tracked
  incomplete$imd_orig_proportion[1] <- NA_real_
  expect_error(vis(incomplete), "missing")
  # The renderer uses only rows supplied by the upstream summary.
  shortened <- tracked[tracked$sample != "S5", ]
  expect_equal(nrow(vis(shortened, plot = "line")$data), nrow(shortened))
  expect_false(inherits(vdj_track(idata), "immunarch_res"))
  expect_null(utils::getS3method("vis", "data.frame", optional = TRUE))
})

test_that("all three tracking geometries build on proportion and count scales", {
  skip_if_not_installed("ggalluvial")
  tracked <- vdj_track(make_track_idata(repertoires = c("sample", "group")),
                       n = 1, output = "df")
  tracked$sample <- factor(tracked$sample, levels = paste0("S", 5:1))
  for (kind in c("alluvial", "area", "line")) {
    for (measure in c("imd_orig_count", "imd_orig_proportion")) {
      p <- vis(tracked, xval = "sample", yval = measure, plot = kind)
      expect_s3_class(p, "ggplot")
      expect_error(ggplot2::ggplot_build(p), NA)
      expect_equal(levels(p$data$sample), paste0("S", 5:1))
    }
  }
  p <- vis(tracked)
  expect_true(is.factor(p$data$imd_repertoire_id))
  tracked$day <- as.integer(sub("S", "", tracked$sample)) * 7
  p <- vis(tracked, xval = "day", plot = "line")
  expect_equal(sort(unique(ggplot2::layer_data(p)$x)), sort(unique(tracked$day)))
})

test_that("faceting does not pool repertoires or invent cross-panel samples", {
  tracked <- vdj_track(make_track_idata(repertoires = c("sample", "group")),
                       n = 1, output = "df")
  tracked$visit <- ifelse(tracked$sample %in% c("S1", "S2"), "Before", tracked$sample)
  expect_error(vis(tracked, xval = "visit"), "one repertoire")
  p <- vis(tracked, xval = "visit", facet = "group", plot = "line")
  expect_error(ggplot2::ggplot_build(p), NA)
  expect_setequal(unique(p$data$imd_repertoire_id[p$data$group == "Case"]),
                  unique(tracked$imd_repertoire_id[tracked$group == "Case"]))
  expect_setequal(unique(p$data$imd_repertoire_id[p$data$group == "Control"]),
                  unique(tracked$imd_repertoire_id[tracked$group == "Control"]))
  expect_error(vis(tracked, yval = "imd_proportion"), "arg")
  expect_error(vis(tracked[0, ]), "rows")
})

test_that("samples and whole panels with no tracked receptors remain at zero", {
  rows <- tibble::tibble(
    sample = c("T0", "T0", "T7", "T14", "T14", "Other"),
    patient = c(rep("P1", 5), "P2"), day = c(0, 0, 7, 14, 14, 0),
    junction = c("AAA", "CCC", "CCC", "AAA", "CCC", "CCC"),
    v_call = "V1", abundance = c(10L, 90L, 100L, 20L, 80L, 100L)
  )
  idata <- make_track_idata(rows, repertoires = c("sample", "patient", "day"))
  tracked <- vdj_track(idata, from = tibble::tibble(junction = "AAA"), output = "df")
  expect_equal(nrow(tracked), 4L)
  expect_setequal(tracked$sample, rows$sample)
  expect_equal(tracked$imd_orig_count[tracked$sample %in% c("T7", "Other")], c(0, 0))
  expect_equal(tracked$imd_orig_proportion[tracked$sample %in% c("T7", "Other")], c(0, 0))
  expect_null(attr(tracked, "repertoires", exact = TRUE))
  expect_equal(tracked$imd_repertoire_id, idata$repertoires$imd_repertoire_id[
    match(tracked$sample, idata$repertoires$sample)
  ])

  for (kind in c("line", "area", "alluvial")) {
    if (kind == "alluvial" && !requireNamespace("ggalluvial", quietly = TRUE)) next
    for (measure in c("imd_orig_count", "imd_orig_proportion")) {
      p <- vis(tracked, xval = "day", facet = "patient", plot = kind, yval = measure)
      expect_setequal(p$data$imd_repertoire_id, idata$repertoires$imd_repertoire_id)
      zero_ids <- idata$repertoires$imd_repertoire_id[idata$repertoires$sample %in% c("T7", "Other")]
      expect_equal(p$data[[measure]][p$data$imd_repertoire_id %in% zero_ids], c(0, 0))
      expect_equal(p$data[[measure]], tracked[[measure]])
      built <- ggplot2::ggplot_build(p)
      expect_setequal(built$layout$layout$patient, c("P1", "P2"))
      expect_equal(sort(unique(built$data[[1]]$x[built$data[[1]]$PANEL == 1])), c(0, 7, 14))
      expect_error(ggplot2::ggplotGrob(p), NA)
    }
  }

  tracked$sample <- factor(tracked$sample, levels = c("T0", "T7", "T14", "Other"))
  p <- vis(tracked, xval = "sample", plot = "line")
  expect_setequal(p$data$sample, rows$sample)
  expect_equal(nrow(p$data), nrow(idata$repertoires))
})
