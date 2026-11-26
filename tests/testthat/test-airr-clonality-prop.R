test_that("proportion summaries and legends follow the supplied bin order", {
  rows <- tibble::tibble(
    sample = rep(c("S1", "S2"), each = 6),
    junction = rep(paste0("r", 1:6), times = 2),
    v_call = "V1",
    abundance = rep(c(1L, 9L, 40L, 100L, 200L, 650L), times = 2)
  )
  idata <- make_track_idata(rows)
  bins <- c(Medium = .1, Hyperexpanded = .5, Small = .01, Large = .2, Rare = .005)
  occupied <- c(Medium = .1, Hyperexpanded = .65, Small = .04,
                Large = .2, Rare = .009, "Ultra-rare" = .001)

  for (thresholds in list(bins, rev(bins))) {
    expected <- c(names(thresholds), "Ultra-rare")
    for (autojoin in c(FALSE, TRUE)) {
      out <- airr_clonality_prop(idata, bins = thresholds, autojoin = autojoin)

      expect_equal(out$clonal_prop_bin, rep(expected, each = 2))
      expect_equal(out$occupied_prop, rep(unname(occupied[expected]), each = 2))

      scale <- ggplot2::ggplot_build(vis(out))$plot$scales$get_scales("fill")
      expect_equal(scale$get_limits(), expected)
    }
  }

  out <- airr_clonality_prop(idata)
  expect_equal(unique(out$clonal_prop_bin), c("Hyperexpanded", "Large"))
})
