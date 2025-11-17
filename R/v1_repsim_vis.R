#' @keywords internal
vis_repsim_intersection_impl <- make_dotplot(
  title_default = "No. of public receptors",
  size_default = "No. receptors",
  fill_default = "No. receptors"
)

register_immunarch_visualisation(vis_repsim_intersection_impl, "repsim", "intersection")


#' @keywords internal
vis_repsim_jaccard_impl <- make_dotplot(
  title_default = "Jaccard similarity index",
  size_default = "Value",
  fill_default = "Value"
)

register_immunarch_visualisation(vis_repsim_jaccard_impl, "repsim", "jaccard")
