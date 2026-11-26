#' Get a list of package updates
#' @concept migration_utility
immunarch_v1_updates <- c(
  "Oct 2025" = immunarch_v1_update_oct_2025,
  "Nov 2026" = immunarch_v1_update_nov_2026
)

.onAttach <- function(libname, pkgname) {
  lines <- c(
    paste0("Hi, this is Vadim Nazarov - author of ", cli::col_green("immunarch"), "."),
    "The January 2027 update will remove legacy functions that currently raise deprecation warnings. ",
    "",
    "- Run {.run [get_immunarch_news()](immunarch::get_immunarch_news())} for migration help.",
    "",
    "- Run {.run [list_immunarch_news()](immunarch::list_immunarch_news())} to list all updates.",
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
  op <- options()
  op.immunarch <- list(IMMUNARCH_AUTOJOIN_DEFAULT)
  names(op.immunarch) <- IMMUNARCH_AUTOJOIN_OPTION
  toset <- !(names(op.immunarch) %in% names(op))
  if (any(toset)) options(op.immunarch[toset])
  invisible()
}
