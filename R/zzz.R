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
  op <- options()
  op.immunarch <- list(
    immunarch.autojoin = FALSE # default
  )
  toset <- !(names(op.immunarch) %in% names(op))
  if (any(toset)) options(op.immunarch[toset])
  invisible()
}
