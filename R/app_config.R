#' Access files in the current app
#'
#' Helper function to robustly resolve file paths within the package installation.
#' This ensures the app can find its data files whether running locally or deployed.
#'
#' @param ... Character vectors, specifying subdirectory and file(s)
#'   within your package. The default, none, returns the root of the app.
#'
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "WormTransgeneBuilder")
}

#' Read App Configuration
#'
#' Retrieves variables from the `golem-config.yml` file, allowing the application
#' to dynamically adjust its behavior based on the active environment (e.g., dev, prod).
#'
#' @param value Character string; the specific config value to retrieve.
#' @param config Character string; the active configuration environment. Defaults to GOLEM_CONFIG_ACTIVE or R_CONFIG_ACTIVE.
#' @param use_parent Logical; whether to scan the parent directory for the config file.
#' @param file Character string; the explicit location of the config file.
#'
#' @noRd
get_golem_config <- function(
    value,
    config = Sys.getenv(
      "GOLEM_CONFIG_ACTIVE",
      Sys.getenv(
        "R_CONFIG_ACTIVE",
        "default"
      )
    ),
    use_parent = TRUE,
    file = app_sys("golem-config.yml")
) {
  config::get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
}
