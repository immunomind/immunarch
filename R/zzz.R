#' @importFrom cli cli_h1 cli_alert_warning cli_par cli_text cli_end cli_alert_info cli_bullets spark_line cli_alert_success start_app cli_inform
#' @importFrom stats runif
immunarch_v1_update_apr_2025 <- function() {
  cli::cli_h1("{cli::col_green('immunarch')} {cli::col_yellow('0.9.x')} – Critical Pre-release Notice")

  cli::cli_alert_warning("Update #1 [Apr 2025] -- Major changes are coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0.0')}!")
  cli::cli_text(cli::col_yellow(cli::spark_line(runif(110, 0, 1))))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "Hi, this is Vadim Nazarov speaking -- author of {cli::col_green('immunarch')}. ",
    "{cli::col_green('immunarch')} is finally graduating from out of the {cli::col_yellow('0.x.y')} development cycle. ",
    "I'm preparing our {cli::col_yellow('1.0.0')} release, which will remain stable and free of sudden changes until we approach {cli::col_yellow('2.0.0')}, along with ",
    "a scientific publication for proper citations. ",
    "Significant changes are coming, and I want to ensure you have everything you need to migrate to the new version."
  )

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("Here’s a preview of what's coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0.0')}:")
  cli::cli_bullets(c(
    "i" = "Some computationally intensive or advanced features (e.g., distance computations, graph-based analyses, dimensionality reduction techniques) will move to separate packages, making {cli::col_green('immunarch')} much more lightweight to install and manage;",
    "i" = "New functions will be introduced instead of the left old ones to make code more readable and maintainable. Legacy functions will remain temporarily, but they won't be updated and will be removed by {cli::col_yellow('~2027')};",
    "i" = "We will discontinue support for most custom file formats because the AIRR ecosystem is now mature enough that the majority of tools adhere to the AIRR standard;",
    "i" = "The package will transition from data frames to the new {cli::col_blue('ImmunData')} structure -- better suited for handling modern larger, more complex, and multi-modal datasets (e.g., single-cell, spatial);",
    "i" = "{cli::col_blue('ImmunData')} is available in the separate {cli::col_blue('immundata')} package, which you can already install via {cli::col_cyan('install.packages(\"immundata\")')};",
    "i" = "The {cli::col_blue('ImmunData')}-based computations will be significantly faster, support datasets larger than RAM, and fully adhere to AIRR Community standards."
  ))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "See the dedicated migration guide for migration on what you can do now and how to prepare for the future:"
  )
  cli::cli_text(">> run {cli::col_cyan('vignette(\"immunarch_v1_migration\")')}, or")
  cli::cli_text(">> visit {cli::col_cyan('https://immunarch.com/articles/immunarch_v1_migration.html')}")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_warning("Keep an eye on the update number and date at the beginning. I’ll share updates, tips, and important dates for our major transformation. If you happen to miss some update, call {cli::col_cyan('immunarch::immunarch_updates()')} to read previous updates.")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_success("Thank you for supporting {cli::col_green('immunarch')} from its early days. Your feedback, contributions, and trust have driven its evolution, and I deeply appreciate it.")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_info("Questions, comments, ideas? I'm available via:")
  cli::cli_text(">> Support email: {cli::col_cyan('support@immunomind.com')}")
  cli::cli_text(">> GitHub tickets: {cli::col_cyan('https://github.com/immunomind/immunarch')}")
  cli::cli_text(">> LinkedIn: {cli::col_cyan('https://www.linkedin.com/in/vdnaz/')}")
}

#' Some description
immunarch_v1_updates <- c(
  "Apr 2025" = immunarch_v1_update_apr_2025
)

#' Get the Latest immunarch Update
#'
#' Retrieves an update message for immunarch.
#'
#' If `datepoint` is set to `"latest"`, the function returns the most recent update.
#' Otherwise, specify the update date key (e.g., `"Apr 2025"`) to retrieve that particular update.
#' If no matching update is found, a warning is issued along with available update keys.
#'
#' @param datepoint A string specifying the update date. Use `"latest"` for the most recent update
#'   or supply a valid date key (e.g., `"Apr 2025"`).
#'
#' @return A character string with the update details or a warning if the key is not found.
#'
#' @seealso [list_immunarch_news()]
#'
#' @export
get_immunarch_news <- function(datepoint = "latest") {
  if (datepoint == "latest") {
    immunarch_v1_updates[[length(immunarch_v1_updates)]]()
  } else if (datepoint %in% names(immunarch_v1_updates)) {
    immunarch_v1_updates[[datepoint]]()
  } else {
    cli::cli_alert_warning("No {datepoint} date in the list of {cli::col_green('immunarch')} updates. Available update names are: {immunarch:::list_immunarch_news()}")
  }
}

#' List Available immunarch Updates
#'
#' Returns the list of available update keys for immunarch v1.
#'
#' @return A character vector containing all the date keys for the available updates.
#'
#' @seealso [get_immunarch_news()]
#'
#' @export
list_immunarch_news <- function() {
  names(immunarch_v1_updates)
}

.onAttach <- function(libname, pkgname) {
  msg <- paste0(
    "Hi, this is Vadim Nazarov speaking -- author of ",
    cli::col_green("immunarch"),
    ".\n",
    cli::col_green("immunarch"),
    " is evolving towards its 1.0.0 release. Soon it will be faster, more user-friendly, and ready for its long-awaited publication. Some functions will no longer be supported or will be replaced with new, more powerful methods.\n",
    "\n -- Please run ",
    cli::col_cyan("get_immunarch_news()"),
    " in your R console to read the latest update and learn what has changed, what's new, how to migrate your code, and what changes are planned for the next update.\n",
    "\n -- Run ",
    cli::col_cyan("list_immunarch_news()"),
    " to list all available updates and catch up on any you may have missed. Latest update: ", cli::col_yellow("#1, Apr 2025"), "\n",
    "\n -- To import the package without this message, run ",
    cli::col_cyan("suppressPackageStartupMessages(library(\"immunarch\"))"),
    "\n",
    "Thank you."
  )

  cli::cli_inform(msg, class = "packageStartupMessage")
}
