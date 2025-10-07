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
#' @concept migration_utility
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
#' @concept migration_utility
#'
#' @export
list_immunarch_news <- function() {
  for (i in seq_along(names(immunarch_v1_updates))) {
    cat(names(immunarch_v1_updates), " -> ", "run immunarch::get_immunarch_news(", '"', names(immunarch_v1_updates), '"', ")", sep = "")
  }
}
