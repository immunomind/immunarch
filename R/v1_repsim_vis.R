#' @keywords internal
make_repsim_dotplot <- function(title) {
  dotplot <- make_dotplot(title, title, title)
  function(.data, ...) {
    plot_data <- .data
    diag(plot_data) <- NA_real_
    plot <- dotplot(plot_data, ...) +
      ggplot2::guides(size = "none")

    values <- plot_data[is.finite(plot_data)]
    if (length(values)) {
      limits <- range(values)
      breaks <- unique(c(limits[1], mean(limits), limits[2]))
      labels <- vapply(breaks, format, character(1), digits = 5, trim = TRUE)
      fill_scale <- plot$scales$get_scales("fill")
      fill_scale$breaks <- breaks
      fill_scale$labels <- labels
    }

    plot
  }
}


#' @keywords internal
vis_repsim_intersection_impl <- make_repsim_dotplot("Shared receptors")

register_immunarch_visualisation(vis_repsim_intersection_impl, "repsim", "intersection")


#' @keywords internal
vis_repsim_jaccard_impl <- make_repsim_dotplot("Jaccard similarity")

register_immunarch_visualisation(vis_repsim_jaccard_impl, "repsim", "jaccard")


#' @keywords internal
vis_repsim_chao_jaccard_impl <- make_repsim_dotplot("Chao-Jaccard similarity")

register_immunarch_visualisation(vis_repsim_chao_jaccard_impl, "repsim", "chao_jaccard")


#' @keywords internal
vis_repsim_morisita_horn_impl <- make_repsim_dotplot("Morisita-Horn similarity")

register_immunarch_visualisation(vis_repsim_morisita_horn_impl, "repsim", "morisita_horn")


#' @keywords internal
vis_repsim_bray_impl <- make_repsim_dotplot("Bray-Curtis similarity")

register_immunarch_visualisation(vis_repsim_bray_impl, "repsim", "bray")
