test_that("repertoire similarity plots omit the diagonal and name the metric", {
  idata <- make_repsim_idata(
    tibble::tribble(
      ~receptor, ~repertoire, ~count,
      "r1", "R1", 2,
      "r2", "R1", 1,
      "r1", "R2", 1,
      "r3", "R2", 2,
      "r2", "R3", 1,
      "r3", "R3", 2
    ),
    normalize = TRUE
  )
  methods <- list(
    "Shared receptors" = repsim_intersection,
    "Jaccard similarity" = repsim_jaccard,
    "Chao-Jaccard similarity" = repsim_chao_jaccard,
    "Morisita-Horn similarity" = repsim_morisita_horn,
    "Bray-Curtis similarity" = repsim_bray
  )

  for (title in names(methods)) {
    result <- methods[[title]](idata)
    plot <- vis(result)

    expect_false(anyNA(diag(result)), info = title)
    expect_equal(sum(is.na(plot$data$value)), nrow(result), info = title)
    expect_null(plot$labels$caption)
    expect_identical(ggplot2::get_labs(plot)$fill, title)
    fill_scale <- plot$scales$get_scales("fill")
    limits <- range(plot$data$value, na.rm = TRUE)
    expect_equal(fill_scale$breaks, unique(c(limits[1], mean(limits), limits[2])))
    expect_false(any(grepl("Min|Max", fill_scale$labels)))
    expect_silent(ggplot2::ggplot_build(plot))
  }
})


test_that("a single-repertoire plot can omit its only cell", {
  one <- matrix(1, 1, 1, dimnames = list("R1", "R1"))
  expect_silent(ggplot2::ggplot_build(vis_repsim_jaccard_impl(one)))
})
