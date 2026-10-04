#' Show Warnings as App Notifications
#'
#' Runs `expr` and, for every warning it raises, shows a warning-style notification
#' to the user (for example "1 restriction site(s) could not be removed ..."). The
#' warning is not swallowed, so it still appears in the server log as well. Errors
#' are not touched.
#'
#' @param expr The code to run.
#' @param notify A function taking the warning text; defaults to a Shiny warning
#'   notification that stays until the user dismisses it.
#' @return The value of `expr`.
#' @noRd
with_user_notifications <- function(expr,
                                    notify = function(msg) {
                                      shiny::showNotification(paste("Warning:", msg), type = "warning", duration = NULL)
                                    }) {
  withCallingHandlers(expr, warning = function(w) notify(conditionMessage(w)))
}
