#' Get a list of package updates
#' @concept migration_utility
immunarch_v1_updates <- c(
  "Oct 2025" = immunarch_v1_update_oct_2025
)

.onAttach <- function(libname, pkgname) {
  lines <- c(
    paste0("Hi, this is Vadim Nazarov - author of ", cli::col_green("immunarch"), "."),
    paste0(
      cli::col_green("immunarch"),
      " is moving toward its 1.0 release. It will be faster, more user-friendly, and ready for publication. ",
      "Some functions will be deprecated or replaced with newer, more powerful methods."
    ),
    "",
    "- Click {.run [get_immunarch_news()](immunarch::get_immunarch_news())} to read the latest update (what changed, what's new, how to migrate, what's next).",
    "",
    "- Click {.run [list_immunarch_news()](immunarch::list_immunarch_news())} to list all updates (latest: {cli::col_yellow('#1, Oct 2025')}).",
    "",
    "- Migration guide: {.url https://immunomind.github.io/docs/tutorials/migration}",
    "",
    "To load the package without this message: {.code suppressPackageStartupMessages(library('immunarch'))}",
    "",
    "- Vadim I. Nazarov"
  )

  msg <- paste(lines, collapse = "\n")
  cli::cli_inform(cli::format_inline(msg), class = "packageStartupMessage")
}


.onLoad <- function(libname, pkgname) {
  # register_immunarch_visualisation(vis_airr_stats_lengths_impl, "airr_stats", "lengths")
  # register_immunarch_visualisation(vis_airr_stats_chains_impl, "airr_stats", "chains")
  # register_immunarch_visualisation(vis_airr_stats_genes_impl, "airr_stats", "genes")
  #
  # register_immunarch_visualisation(
  #   vis_airr_diversity_dxx_impl,
  #   "airr_diversity",
  #   "dxx"
  # )
  # register_immunarch_visualisation(
  #   vis_airr_diversity_chao1_impl,
  #   "airr_diversity",
  #   "chao1"
  # )
  # register_immunarch_visualisation(
  #   vis_airr_diversity_shannon_impl,
  #   "airr_diversity",
  #   "shannon"
  # )
  # register_immunarch_visualisation(
  #   vis_airr_diversity_pielou_impl,
  #   "airr_diversity",
  #   "pielou"
  # )
  # register_immunarch_visualisation(
  #   vis_airr_diversity_index_impl,
  #   "airr_diversity",
  #   "index"
  # )
  #
  #
  # register_immunarch_visualisation(vis_airr_public_intersection_impl, "airr_public", "intersection")
  # register_immunarch_visualisation(vis_airr_public_jaccard_impl, "airr_public", "jaccard")

  op <- options()
  op.immunarch <- list(
    immunarch.autojoin = FALSE # default
  )
  toset <- !(names(op.immunarch) %in% names(op))
  if (any(toset)) options(op.immunarch[toset])
  invisible()
}
