test_that("vis() for airr_desc_chains shows repertoire richness", {
  idata <- make_aggregated_test_idata("Therapy")
  chains <- airr_desc_chains(idata)
  p <- vis(chains)

  expect_s3_class(p, "ggplot")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p)), NA)
  expect_equal(p$labels$title, "Receptor richness by repertoire")
  expect_equal(nrow(p$data), length(unique(chains$imd_repertoire_id)))

  point_layer <- which(vapply(
    p$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)
  ))
  points <- ggplot2::layer_data(p, point_layer)
  expect_equal(sort(10^points$x), sort(chains$n_barcodes), tolerance = 1e-8)
  expect_equal(sort(10^points$y), sort(chains$n_receptors), tolerance = 1e-8)
  expect_setequal(
    p$layers[[2]]$data$label,
    c("100% (1:1)", "50% (1:2)", "25% (1:4)")
  )
})

test_that("vis() labels and colours one point per repertoire across loci", {
  idata <- make_aggregated_test_idata("Therapy")
  chains <- airr_desc_chains(idata)
  chains_by_locus <- dplyr::bind_rows(
    chains, dplyr::mutate(chains, locus = "TRA")
  )

  p <- vis_airr_desc_chains_impl(
    chains_by_locus, label = "Therapy", fill = "Therapy"
  )

  expect_equal(nrow(p$data), length(unique(chains$imd_repertoire_id)))
  expect_setequal(p$data$plot_label, unique(chains$Therapy))
  expect_equal(rlang::as_label(p$mapping$colour), "Therapy")
  expect_error(suppressWarnings(ggplot2::ggplot_build(p)), NA)
})
