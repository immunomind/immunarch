#' @keywords internal
vis_airr_public_intersection_impl <- make_dotplot(
  title_default = "No. of public receptors",
  size_default = "No. receptors",
  fill_default = "No. receptors"
)

register_immunarch_visualisation(vis_airr_public_intersection_impl, "airr_public", "intersection")


#' @keywords internal
vis_airr_public_jaccard_impl <- make_dotplot(
  title_default = "Jaccard similarity index",
  size_default = "Value",
  fill_default = "Value"
)

register_immunarch_visualisation(vis_airr_public_jaccard_impl, "airr_public", "jaccard")
