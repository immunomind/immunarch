#' @keywords internal
vis_vdj_track_impl <- function(
    .data,
    plot = c("alluvial", "area", "line"),
    xval = "imd_repertoire_id",
    yval = c("imd_orig_proportion", "imd_orig_count"),
    facet = NULL,
    title = "Receptor tracking") {
  checkmate::assert_data_frame(.data, min.rows = 1)
  plot <- match.arg(plot)
  yval <- match.arg(yval)
  checkmate::assert_string(xval)
  checkmate::assert_string(facet, null.ok = TRUE)
  checkmate::assert_string(title)
  columns <- unique(c("imd_receptor_id", "imd_repertoire_id", xval, yval, facet))
  checkmate::assert_subset(c(columns, "imd_receptor_label"), names(.data))
  labels <- .data[!duplicated(.data$imd_receptor_id), c("imd_receptor_id", "imd_receptor_label")]
  receptor_labels <- stats::setNames(labels$imd_receptor_label, labels$imd_receptor_id)
  plot_data <- .data[columns]
  checkmate::assert_numeric(plot_data[[yval]], lower = 0, finite = TRUE,
                            any.missing = FALSE)
  if (anyDuplicated(plot_data[c("imd_repertoire_id", "imd_receptor_id")])) {
    cli::cli_abort("The tracking summary must contain one row per receptor and repertoire.")
  }
  axis_columns <- unique(c("imd_repertoire_id", xval, facet))
  axes <- unique(plot_data[axis_columns])
  if (anyNA(axes[c(xval, facet)])) {
    cli::cli_abort("Plotting metadata must not contain missing values.")
  }
  if (anyDuplicated(axes["imd_repertoire_id"]) ||
      anyDuplicated(axes[unique(c(xval, facet))])) {
    cli::cli_abort("Each x value within a facet must identify one repertoire. Define repertoires upstream or use {.arg facet}.")
  }
  if (xval == "imd_repertoire_id") {
    plot_data[[xval]] <- factor(plot_data[[xval]], levels = sort(unique(plot_data[[xval]])))
  }
  plot_data$imd_receptor_id <- factor(plot_data$imd_receptor_id,
    levels = sort(unique(plot_data$imd_receptor_id)))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = !!rlang::sym(xval), y = !!rlang::sym(yval),
    group = imd_receptor_id, fill = imd_receptor_id,
    colour = imd_receptor_id
  ))
  if (plot == "alluvial") {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggalluvial} is required for alluvial plots.")
    }
    p <- p + ggplot2::aes(stratum = imd_receptor_id, alluvium = imd_receptor_id) +
      ggalluvial::geom_alluvium(width = 0.18, alpha = 0.65, colour = NA) +
      ggalluvial::geom_stratum(
        data = plot_data[plot_data[[yval]] > 0, , drop = FALSE],
        width = 0.18, colour = "white", linewidth = 0.3
      )
  } else if (plot == "area") {
    p <- p + ggplot2::geom_area(colour = "white", linewidth = 0.3)
  } else {
    p <- p + ggplot2::geom_line(linewidth = 0.85) +
      ggplot2::geom_point(size = 2)
  }
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet)), scales = "free_x")
  }
  p + ggplot2::scale_y_continuous(
    labels = if (yval == "imd_orig_proportion") scales::label_percent() else scales::label_number(),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
    ggplot2::scale_fill_viridis_d(option = "H", name = "Receptor", labels = receptor_labels) +
    ggplot2::scale_colour_viridis_d(option = "H", name = "Receptor", labels = receptor_labels) +
    ggplot2::guides(fill = if (plot == "line") "none" else "legend",
                    colour = if (plot == "line") "legend" else "none") +
    ggplot2::labs(title = title, x = xval,
      y = if (yval == "imd_orig_proportion") "Original repertoire proportion" else "Original receptor count",
      caption = "Unobserved combinations shown at zero; not evidence of biological absence.") +
    ggthemes::theme_few() +
    ggplot2::theme(legend.position = "right", plot.caption = ggplot2::element_text(colour = "grey45"))
}

register_immunarch_visualisation(vis_vdj_track_impl, "vdj", "track")
