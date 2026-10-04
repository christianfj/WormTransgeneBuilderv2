#' Run the Shiny Application
#'
#' Initializes and launches the WormTransgeneBuilder Shiny application.
#' This function acts as the primary entry point for the package.
#'
#' @param onStart A function that will be called before the app is actually run.
#'   This is only needed for `shinyAppObj`, since in the `shinyAppDir` case, a
#'   `global.R` file can be used for this purpose.
#' @param options Named options that should be passed to the `runApp` call
#'   (these can be any of the following: "port", "launch.browser", "host", "quiet",
#'   "display.mode" and "test.mode").
#' @param enableBookmarking Can be one of "url", "server", or "disable".
#' @param uiPattern A regular expression that will be applied to each `GET` request
#'   to determine whether the `ui` should be used to handle the request.
#' @param ... Additional arguments to pass to `golem_opts`.
#'   See `?golem::get_golem_options` for more details.
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
    onStart = NULL,
    options = list(launch.browser = TRUE),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
