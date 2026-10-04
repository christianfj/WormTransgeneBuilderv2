# Small housekeeping checks: the version numbers agree, and the search for the
# external RNAfold program contains no path that only exists on one person's computer.

test_that("the golem config version matches the package version in DESCRIPTION", {
  expect_equal(get_golem_config("golem_version"),
               as.character(utils::packageVersion("WormTransgeneBuilder")))
})

# ---- find_rnafold() ------------------------------------------------------------

make_fake_program <- function(dir, name = "RNAfold") {
  path <- file.path(dir, name)
  writeLines("#!/bin/sh", path)
  path
}

test_that("find_rnafold returns the program found on the PATH", {
  expect_equal(find_rnafold(on_path = "/usr/bin/RNAfold", fallbacks = character(0)), "/usr/bin/RNAfold")
})

test_that("the PATH result wins over the fallback locations", {
  dir <- withr::local_tempdir()
  fake <- make_fake_program(dir)
  expect_equal(find_rnafold(on_path = "/usr/bin/RNAfold", fallbacks = fake), "/usr/bin/RNAfold")
})

test_that("find_rnafold uses the first fallback that exists", {
  d1 <- withr::local_tempdir(); d2 <- withr::local_tempdir()
  second <- make_fake_program(d2)
  expect_equal(find_rnafold(on_path = "", fallbacks = c(file.path(d1, "RNAfold"), second)), second)
  first <- make_fake_program(d1)
  expect_equal(find_rnafold(on_path = "", fallbacks = c(first, second)), first)
})

test_that("find_rnafold expands ~ in fallback locations", {
  home <- withr::local_tempdir()
  fake <- make_fake_program(home)
  withr::local_envvar(HOME = home)
  expect_equal(find_rnafold(on_path = "", fallbacks = "~/RNAfold"), path.expand("~/RNAfold"))
})

test_that("find_rnafold returns an empty string when nothing is found", {
  expect_identical(find_rnafold(on_path = "", fallbacks = "/no/such/place/RNAfold"), "")
  expect_identical(find_rnafold(on_path = "", fallbacks = character(0)), "")
})

test_that("find_rnafold returns a plain string, not a named vector", {
  expect_null(names(find_rnafold(on_path = c(RNAfold = "/usr/bin/RNAfold"), fallbacks = character(0))))
})

test_that("the default search list holds no personal or absolute home-directory path", {
  defaults <- eval(formals(find_rnafold)$fallbacks)
  expect_gt(length(defaults), 0)
  expect_false(any(grepl("^/Users/|^/home/", defaults)))
})

test_that("optimize_rbs warns and returns the input when RNAfold cannot be found", {
  testthat::local_mocked_bindings(find_rnafold = function(...) "")
  start <- substr(glo_cds(), 1, 39)
  expect_warning(out <- optimize_rbs(start, "Ubiq"), "RNAfold binary not found")
  expect_identical(out, start)
})
