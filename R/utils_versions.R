#' Versions of the Software Behind a Result
#'
#' Returns the version of this app, of R and of ViennaRNA (`RNAfold`), so that a result can be tied to the
#' exact software that produced it. It is shown in the About tab and printed when the app starts.
#'
#' @param rnafold Path to `RNAfold`, or `""` if it is not installed (the default looks for it as the
#'   ribosome-binding-site step does).
#' @return A named character vector with elements `app`, `R` and `ViennaRNA`; ViennaRNA is
#'   "not installed" when `RNAfold` cannot be found or run.
#' @noRd
software_versions <- function(rnafold = find_rnafold()) {
  vienna <- "not installed"
  if (nzchar(rnafold)) {
    out <- tryCatch(suppressWarnings(system2(rnafold, "--version", stdout = TRUE, stderr = FALSE)),
                    error = function(e) character(0))
    version <- sub("^RNAfold\\s+", "", out[1])
    if (length(out) > 0 && !is.na(version) && nzchar(version)) vienna <- version
  }
  c(app = as.character(utils::packageVersion("WormTransgeneBuilder")),
    R = paste(R.version$major, R.version$minor, sep = "."),
    ViennaRNA = vienna)
}
