#' @keywords internal
vis_airr_desc_chains_impl <- make_dynam_col_plot(
  y_default = "n_chains",
  title_default = "No. chains per sample",
  position = "dodge"
)

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
