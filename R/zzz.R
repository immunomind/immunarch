#' Get a list of package updates
#' @concept migration_utility
immunarch_v1_updates <- c(
  "Sep 2025" = immunarch_v1_update_sep_2025
)

.onAttach <- function(libname, pkgname) {
  msg <- paste0(
    "Hi, this is Vadim Nazarov speaking -- author of ",
    cli::col_green("immunarch"),
    ".\n",
    cli::col_green("immunarch"),
    " is evolving towards its 1.0 release. Soon it will be faster, more user-friendly, and ready for its long-awaited publication. Some functions will no longer be supported or will be replaced with new, more powerful methods.\n",
    "\n -- Please click on ",
    cli::col_cyan("{.run [get_immunarch_news()](immunarch::get_immunarch_news())}"),
    " or run it in your R console to read the latest update and learn what has changed, what's new, how to migrate your code, and what changes are planned for the next update.\n",
    "\n -- Click on ",
    cli::col_cyan("{.run [list_immunarch_news()](immunarch::list_immunarch_news())}"),
    " or run it to list all available updates and catch up on any you may have missed. Latest update: ",
    cli::col_yellow("#1, Sep 2025"),
    "\n",
    "\n -- To import the package without this message, run ",
    cli::col_cyan("suppressPackageStartupMessages(library(\"immunarch\"))"),
    "\n",
    "\nMigration guide is available online:\n\n-- in R: ",
    "\n",
    " {.url https://immunomind.github.io/docs/tutorials/migration}",
    "\n\nThank you.\n",
    "\n- Vadim I. Nazarov"
  )

  cli::cli_inform(msg, class = "packageStartupMessage")

  # Show registered methods?
}
