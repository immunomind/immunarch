#' @keywords internal
vis_airr_desc_chains_impl <- function(
  .data,
  label = "imd_repertoire_id",
  colour = NULL,
  fill = NULL,
  title = "Receptor richness by repertoire"
) {
  checkmate::assert_data_frame(.data)
  checkmate::assert_character(label, min.len = 1, any.missing = FALSE)
  checkmate::assert_string(colour, null.ok = TRUE)
  checkmate::assert_string(fill, null.ok = TRUE)
  if (is.null(colour)) colour <- fill
  checkmate::assert_subset(
    c("imd_repertoire_id", "n_barcodes", "n_receptors", label, colour),
    choices = names(.data)
  )

  # The repertoire totals are repeated for each locus in airr_desc_chains().
  plot_data <- .data[!duplicated(.data$imd_repertoire_id), , drop = FALSE]
  plot_data$plot_label <- apply(
    plot_data[, label, drop = FALSE], 1, paste, collapse = " / "
  )
  guide_data <- data.frame(
    x = max(plot_data$n_barcodes) * 0.7,
    ratio = c(1, 0.5, 0.25),
    label = c("100% (1:1)", "50% (1:2)", "25% (1:4)")
  )
  guide_data$y <- guide_data$x * guide_data$ratio

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n_barcodes, y = n_receptors)
  ) +
    ggplot2::geom_abline(
      intercept = log10(c(1, 0.5, 0.25)), slope = 1,
      colour = "grey80", linetype = "dashed"
    ) +
    ggplot2::geom_label(
      data = guide_data,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, hjust = 1,
      size = 3, colour = "grey45", fill = "white"
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_comma(), expand = ggplot2::expansion(mult = 0.12)
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_comma(), expand = ggplot2::expansion(mult = 0.12)
    ) +
    ggplot2::labs(
      x = "Repertoire size (n_barcodes)", y = "Distinct receptors",
      title = title,
      subtitle = "Dashed lines show the overall ratio of distinct receptors to barcodes"
    ) +
    ggthemes::theme_few()

  if (!is.null(colour)) {
    p <- p + ggplot2::aes(colour = !!rlang::sym(colour))
  }

  p +
    ggplot2::geom_point(size = 2.5) +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = plot_label),
      size = 3, show.legend = FALSE,
      max.overlaps = Inf, seed = 1
    )
}

register_immunarch_visualisation(vis_airr_desc_chains_impl, "airr_desc", "chains")


#' @keywords internal
vis_airr_desc_lengths_impl <- make_fixed_col_plot(
  x_col = "seq_len",
  y_col = "prop",
  title = "CDR3 length distribution", xlab = "CDR3 length", ylab = "Proportion"
)

register_immunarch_visualisation(vis_airr_desc_lengths_impl, "airr_desc", "lengths")


#' @keywords internal
vis_airr_desc_genes_impl <- make_dotplot(
  title_default = "Gene usage",
  size_default = "Proportion (all data)",
  fill_default = "No. Receptors"
)

register_immunarch_visualisation(vis_airr_desc_genes_impl, "airr_desc", "genes")
