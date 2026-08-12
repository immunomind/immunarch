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


#' @keywords internal
vis_repsim_chao_jaccard_impl <- make_dotplot(
  title_default = "Chao-Jaccard similarity index",
  size_default = "Value",
  fill_default = "Value"
)

register_immunarch_visualisation(vis_repsim_chao_jaccard_impl, "repsim", "chao_jaccard")


#' @keywords internal
vis_repsim_morisita_horn_impl <- make_dotplot(
  title_default = "Morisita-Horn similarity index",
  size_default = "Value",
  fill_default = "Value"
)

register_immunarch_visualisation(vis_repsim_morisita_horn_impl, "repsim", "morisita_horn")


#' @keywords internal
vis_repsim_bray_impl <- make_dotplot(
  title_default = "Bray-Curtis dissimilarity",
  size_default = "Distance",
  fill_default = "Distance"
)

register_immunarch_visualisation(vis_repsim_bray_impl, "repsim", "bray")
