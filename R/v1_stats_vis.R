#' @keywords internal
vis_airr_stats_chains_impl <- make_dynam_col_plot(
  y_default = "n_receptors",
  title_default = "No. receptors per sample",
  position = "dodge"
)

register_immunarch_visualisation(vis_airr_stats_chains_impl, "airr_stats", "chains")


#' @keywords internal
vis_airr_stats_lengths_impl <- make_fixed_col_plot(
  x_col = "seq_len",
  y_col = "prop",
  title = "CDR3 length distribution", xlab = "CDR3 length", ylab = "Proportion"
)

register_immunarch_visualisation(vis_airr_stats_lengths_impl, "airr_stats", "lengths")


#' @keywords internal
vis_airr_stats_genes_impl <- make_dotplot(
  title_default = "Gene usage",
  size_default = "Proportion (all data)",
  fill_default = "No. Receptors"
)

register_immunarch_visualisation(vis_airr_stats_genes_impl, "airr_stats", "genes")
