#' Run Code with a Temporary Random Seed
#'
#' Evaluates `code` with `set.seed(seed)` in effect and then puts R's random number
#' generator back exactly as it was, even if `code` fails. Unlike calling
#' `set.seed()` directly, this cannot affect any random choice made afterwards, by
#' the same user or by anyone else using the same R process (a Shiny server shares
#' one R process between all sessions).
#'
#' @param seed Integer seed, or `NULL` to run `code` without touching the seed.
#' @param code The code to run. It is evaluated after the seed is set.
#' @return The value of `code`.
#' @noRd
with_local_seed <- function(seed, code) {
  if (is.null(seed)) return(code)

  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  set.seed(seed)
  code
}
