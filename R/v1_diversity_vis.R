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
