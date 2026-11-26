test_that("airr_clonality_line returns proportional top ranks per repertoire", {
  repertoire_col <- immundata::imd_schema("repertoire")
  count_col <- immundata::imd_schema("count")
  proportion_col <- immundata::imd_schema("proportion")

  counts <- rep(12:1, times = 2)
  repertoires <- rep(c("R1", "R2"), each = 12)
  receptors <- paste0(repertoires, "_r", rep(seq_len(12), times = 2))

  idata <- make_test_repertoire_idata(
    receptors = receptors,
    repertoires = repertoires,
    counts = counts,
    proportions = counts / sum(12:1)
  )

  out <- airr_clonality_line(idata, limit = 10, autojoin = FALSE)

  expect_true(proportion_col %in% names(out))
  expect_equal(
    as.integer(table(out[[repertoire_col]])),
    c(10L, 10L)
  )

  for (repertoire in c("R1", "R2")) {
    repertoire_out <- out[out[[repertoire_col]] == repertoire, ]

    expect_equal(repertoire_out$index, seq_len(10))
    expect_equal(repertoire_out[[count_col]], 12:3)
    expect_equal(repertoire_out[[proportion_col]], (12:3) / sum(12:1))
  }

  count_plot <- ggplot2::ggplot(
    out,
    ggplot2::aes(
      x = index,
      y = imd_count,
      colour = imd_repertoire_id,
      group = imd_repertoire_id
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10()

  proportion_plot <- ggplot2::ggplot(
    out,
    ggplot2::aes(
      x = index,
      y = imd_proportion,
      colour = imd_repertoire_id,
      group = imd_repertoire_id
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10()

  expect_silent(ggplot2::ggplot_build(count_plot))
  expect_silent(ggplot2::ggplot_build(proportion_plot))

  expect_error(
    airr_clonality_line(idata, limit = 9, autojoin = FALSE)
  )
  expect_error(
    airr_clonality_line(idata, limit = 10.5, autojoin = FALSE)
  )
})
