#' Get the Latest immunarch Update
#'
#' Retrieves an update message for immunarch.
#'
#' If `datepoint` is set to `"latest"`, the function returns the most recent update.
#' Otherwise, specify the update date key (e.g., `"Apr 2025"`) to retrieve that particular update.
#' If no matching update is found, a warning is issued along with available update keys.
#'
#' @param datepoint A single character string specifying the update date. Use
#'   `"latest"` for the most recent update or supply a valid date key (e.g.,
#'   `"Apr 2025"`).
#'
#' @return A character string with the update details or a warning if the key is not found.
#'
#' @seealso [list_immunarch_news()]
#'
#' @concept migration_utility
#'
#' @export
get_immunarch_news <- function(datepoint = "latest") {
  if (datepoint == "latest") {
    immunarch_v1_updates[[length(immunarch_v1_updates)]]()
  } else if (datepoint %in% names(immunarch_v1_updates)) {
    immunarch_v1_updates[[datepoint]]()
  } else {
    cli::cli_alert_warning("No immunarch update for {.val {datepoint}}. Available dates are {.val {names(immunarch_v1_updates)}}.")
  }
}

#' List Available immunarch Updates
#'
#' Prints each available update date and its command. Click the command to
#' run it in consoles that support clickable R commands.
#'
#' @return Invisibly returns `NULL`.
#'
#' @seealso [get_immunarch_news()]
#'
#' @concept migration_utility
#'
#' @export
list_immunarch_news <- function() {
  for (datepoint in names(immunarch_v1_updates)) {
    command <- sprintf('immunarch::get_immunarch_news("%s")', datepoint)
    line <- cli::format_inline("{datepoint} -> run {.run [{command}]({command})}")
    cat(line, "\n", sep = "")
  }
}
