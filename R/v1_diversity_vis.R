#' @keywords internal
vis_airr_diversity_dxx_impl <- make_dynam_col_plot(
  y_default     = "dxx",
  title_default = "Coverage diversity (Dxx)",
  position      = "dodge"
)
register_immunarch_visualisation(
  vis_airr_diversity_dxx_impl,
  "airr_diversity",
  "dxx"
)


#' @keywords internal
vis_airr_diversity_chao1_impl <- make_dynam_col_plot(
  y_default     = "Estimator",
  title_default = "Chao1 richness estimator",
  position      = "dodge"
)
register_immunarch_visualisation(
  vis_airr_diversity_chao1_impl,
  "airr_diversity",
  "chao1"
)


#' @keywords internal
vis_airr_diversity_shannon_impl <- make_dynam_col_plot(
  y_default     = "shannon",
  title_default = "Shannon entropy (bits)",
  position      = "dodge"
)
register_immunarch_visualisation(
  vis_airr_diversity_shannon_impl,
  "airr_diversity",
  "shannon"
)


#' @keywords internal
vis_airr_diversity_pielou_impl <- make_dynam_col_plot(
  y_default     = "pielou",
  title_default = "Pielou evenness",
  position      = "dodge"
)
register_immunarch_visualisation(
  vis_airr_diversity_pielou_impl,
  "airr_diversity",
  "pielou"
)


#' @keywords internal
vis_airr_diversity_index_impl <- make_dynam_col_plot(
  y_default     = "hill_number",
  title_default = "Hill diversity index (q = 1)",
  position      = "dodge"
)
register_immunarch_visualisation(
  vis_airr_diversity_index_impl,
  "airr_diversity",
  "index"
)


#' @keywords internal
vis_airr_diversity_rarefaction_impl <- function(.data,
                                                color = immundata::imd_schema("repertoire"),
                                                show_ci = TRUE,
                                                log = FALSE,
                                                ...) {
  checkmate::assert_data_frame(.data)
  checkmate::assert_string(color)
  checkmate::assert_choice(color, names(.data))
  checkmate::assert_logical(show_ci, len = 1)
  checkmate::assert_logical(log, len = 1)
  checkmate::assert_subset(c("size", "mean"), names(.data))

  has_ci <- all(c("q_low", "q_high") %in% names(.data))
  has_type <- "type" %in% names(.data)
  plot_data <- .data
  plot_data[[color]] <- as.character(plot_data[[color]])

  if (has_type) {
    plot_data$type <- factor(plot_data$type, levels = c("interpolation", "extrapolation"))
  }

  is_norm <- max(plot_data$size, na.rm = TRUE) <= 1 + sqrt(.Machine$double.eps) &&
    max(plot_data$mean, na.rm = TRUE) <= 1 + sqrt(.Machine$double.eps)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$size,
      y = .data$mean,
      colour = .data[[color]],
      group = if (has_type) interaction(.data[[color]], .data$type, drop = TRUE) else .data[[color]]
    )
  )

  if (show_ci && has_ci) {
    ribbon_data <- plot_data
    if (has_type) {
      ribbon_data <- ribbon_data |> dplyr::filter(.data$type == "interpolation")
    }

    p <- p + ggplot2::geom_ribbon(
      data = ribbon_data,
      ggplot2::aes(
        x = .data$size,
        ymin = .data$q_low,
        ymax = .data$q_high,
        fill = .data[[color]],
        group = .data[[color]]
      ),
      inherit.aes = FALSE,
      alpha = 0.15,
      colour = NA
    )
  }

  if (has_type) {
    p <- p + ggplot2::geom_line(ggplot2::aes(linetype = .data$type), linewidth = 0.8, na.rm = TRUE)
  } else {
    p <- p + ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE)
  }

  p <- p +
    ggplot2::ggtitle("Rarefaction analysis") +
    ggplot2::xlab(if (is_norm) "Sample size (proportion of clones)" else "Sample size (clones)") +
    ggplot2::ylab(if (is_norm) "Estimated richness (relative to observed)" else "Estimated richness (unique receptors)")

  if (isTRUE(log)) {
    p <- p + ggplot2::scale_x_log10()
  }

  k <- length(unique(plot_data[[color]]))
  if (requireNamespace("ggsci", quietly = TRUE) && k <= 11) {
    p <- p +
      ggsci::scale_color_locuszoom(name = color) +
      ggsci::scale_fill_locuszoom(name = color)
  } else {
    p <- p +
      ggplot2::scale_colour_viridis_d(option = "H", name = color) +
      ggplot2::scale_fill_viridis_d(option = "H", name = color)
  }

  if (requireNamespace("ggthemes", quietly = TRUE)) {
    p <- p + ggthemes::theme_few()
  }

  p
}

register_immunarch_visualisation(
  vis_airr_diversity_rarefaction_impl,
  "airr_diversity",
  "rarefaction"
)
