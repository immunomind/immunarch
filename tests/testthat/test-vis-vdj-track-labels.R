test_that("feature labels are used by all tracking plot types without changing identity", {
  tracked <- vdj_track(make_track_idata(), from = tibble::tibble(junction = c("AAA", "CCC")),
                       output = "df")
  expect_setequal(tracked$imd_receptor_label, c("AAA V1", "AAA V2", "CCC V2"))
  for (kind in c("alluvial", "area", "line")) {
    if (kind == "alluvial" && !requireNamespace("ggalluvial", quietly = TRUE)) next
    p <- vis(tracked, plot = kind)
    scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales(if (kind == "line") "colour" else "fill")
    expect_setequal(unname(scale$get_labels()), unique(tracked$imd_receptor_label))
  }
  tracked$imd_receptor_label <- "Same display label"
  built <- ggplot2::ggplot_build(vis(tracked, plot = "line"))
  expect_equal(length(unique(built$data[[1]]$group)), length(unique(tracked$imd_receptor_id)))
})

test_that("paired receptors get complete feature labels on the default plot", {
  rows <- tibble::tibble(
    sample = rep(c("S1", "S1", "S2"), each = 2),
    cell = rep(c("a1", "a2", "b1"), each = 2),
    junction = c("AAA", "CCC", "AAA", "GGG", "AAA", "CCC"),
    v_call = c("HV1", "LV1", "HV1", "LV2", "HV1", "LV1"),
    locus = rep(c("IGH", "IGL"), 3), umi = 10L
  )
  idata <- make_track_idata(rows, paired = TRUE)
  tracked <- vdj_track(idata, output = "df")
  expected <- c("IGH AAA HV1 + IGL CCC LV1", "IGH AAA HV1 + IGL GGG LV2")
  expect_setequal(tracked$imd_receptor_label, expected)
  kinds <- c("line", "area", if (requireNamespace("ggalluvial", quietly = TRUE)) "alluvial")
  for (kind in kinds) {
    p <- vis(tracked, plot = kind)
    built <- ggplot2::ggplot_build(p)
    scale <- built$plot$scales$get_scales(if (kind == "line") "colour" else "fill")
    expect_setequal(unname(scale$get_labels()), expected)
    expect_equal(length(unique(built$data[[1]]$group)), 2L)
  }
})
