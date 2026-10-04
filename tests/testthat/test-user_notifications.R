# with_user_notifications() turns the warnings an optimization raises into pop-up
# notifications in the app, so a user sees "1 restriction site could not be removed"
# instead of it appearing only in the server's console log.

test_that("each warning is passed to the notify function", {
  msgs <- character(0)
  suppressWarnings(with_user_notifications({
    warning("first problem", call. = FALSE)
    warning("second problem", call. = FALSE)
  }, notify = function(m) msgs <<- c(msgs, m)))
  expect_equal(msgs, c("first problem", "second problem"))
})

test_that("the value of the code is returned", {
  expect_equal(with_user_notifications(6 * 7, notify = function(m) NULL), 42)
  res <- suppressWarnings(with_user_notifications({ warning("x"); "done" }, notify = function(m) NULL))
  expect_equal(res, "done")
})

test_that("nothing is notified when there is no warning", {
  msgs <- character(0)
  with_user_notifications(1 + 1, notify = function(m) msgs <<- c(msgs, m))
  expect_length(msgs, 0)
})

test_that("the warning still reaches R, so it also appears in the server log", {
  expect_warning(with_user_notifications(warning("visible"), notify = function(m) NULL), "visible")
})

test_that("errors are not swallowed", {
  expect_error(with_user_notifications(stop("boom"), notify = function(m) NULL), "boom")
})

test_that("the default notifier shows a warning-style notification", {
  shown <- list()
  testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "message", duration = 5) {
    shown[[length(shown) + 1]] <<- list(ui = ui, type = type, duration = duration)
  }, .package = "shiny")
  suppressWarnings(with_user_notifications(warning("look here")))
  expect_length(shown, 1)
  expect_match(shown[[1]]$ui, "look here")
  expect_equal(shown[[1]]$type, "warning")
  expect_null(shown[[1]]$duration)  # stays until dismissed
})

test_that("the app passes an optimization warning through and still finishes", {
  # An enzyme name that is not in the table cannot be chosen in the app, but the server
  # must cope with it: it warns and still produces a sequence.
  cds <- paste0("ATG", strrep("GCTAAG", 30), "TAA")
  expect_warning(seq <- run_server_once("1", cds, enzymes = "NotAnEnzyme", introns = NULL), "NotAnEnzyme")
  expect_gt(nchar(seq), 0)
})
