# The restriction-enzyme menu is grouped: Commonly used, Golden Gate, All other enzymes. Every enzyme
# is in exactly one group (BsaI and BsmBI, which are both common and Golden Gate, are listed only
# under Golden Gate), each group is alphabetical, and there is no "select all".

menu_data <- function() utils::read.csv(system.file("extdata", "Enzymes.csv", package = "WormTransgeneBuilder"),
                                        stringsAsFactors = FALSE)

# ---- the grouping rule ---------------------------------------------------------------------------

test_that("the menu has the three groups in the order Commonly used, Golden Gate, All other enzymes", {
  expect_identical(names(enzyme_menu_choices(menu_data())), c("Commonly used", "Golden Gate", "All other enzymes"))
})

test_that("every enzyme is in exactly one group", {
  groups <- enzyme_menu_choices(menu_data())
  all_listed <- unlist(groups, use.names = FALSE)
  expect_equal(anyDuplicated(all_listed), 0)
  expect_setequal(all_listed, menu_data()$Enzyme)
  expect_length(all_listed, nrow(menu_data()))
})

test_that("Golden Gate holds every Golden Gate enzyme, including BsaI and BsmBI", {
  golden <- enzyme_menu_choices(menu_data())[["Golden Gate"]]
  expect_setequal(golden, c("BbsI", "BsaI", "BsmBI", "BspQI", "BtgZI", "Esp3I", "PaqCI", "SapI"))
})

test_that("BsaI and BsmBI are listed only under Golden Gate, not also under Commonly used", {
  groups <- enzyme_menu_choices(menu_data())
  expect_false(any(c("BsaI", "BsmBI") %in% groups[["Commonly used"]]))
  expect_false(any(c("BsaI", "BsmBI") %in% groups[["All other enzymes"]]))
})

test_that("Commonly used is the common set without the Golden Gate enzymes", {
  common <- enzyme_menu_choices(menu_data())[["Commonly used"]]
  expect_setequal(common, setdiff(common_enzymes, c("BsaI", "BsmBI")))
  expect_length(common, 22)
})

test_that("All other enzymes holds only enzymes that are neither common nor Golden Gate", {
  other <- enzyme_menu_choices(menu_data())[["All other enzymes"]]
  expect_false(any(other %in% common_enzymes))
  expect_false(any(other %in% enzyme_menu_choices(menu_data())[["Golden Gate"]]))
  expect_length(other, nrow(menu_data()) - 22 - 8)
})

test_that("each group is in alphabetical order, ignoring case", {
  for (g in enzyme_menu_choices(menu_data())) expect_identical(g, g[order(tolower(g))])
})

test_that("a Golden Gate enzyme marked common is still listed once, under Golden Gate", {
  e <- data.frame(Enzyme = c("A", "B", "C"), GG = c(TRUE, FALSE, FALSE), Common = c(TRUE, TRUE, FALSE))
  expect_identical(enzyme_menu_choices(e), list("Commonly used" = "B", "Golden Gate" = "A", "All other enzymes" = "C"))
})

test_that("empty groups are left out, and missing flag columns count as FALSE", {
  expect_identical(enzyme_menu_choices(data.frame(Enzyme = c("X", "Y"))), list("All other enzymes" = c("X", "Y")))
  e <- data.frame(Enzyme = c("X", "Y"), GG = c(FALSE, TRUE))
  expect_identical(names(enzyme_menu_choices(e)), c("Golden Gate", "All other enzymes"))
})

test_that("no group offers to select everything", {
  expect_false(any(grepl("all", names(enzyme_menu_choices(menu_data())[c("Commonly used", "Golden Gate")]), ignore.case = TRUE)))
  # The only group with "all" in its name is the catch-all for the rest, which is a list, not an action.
  expect_true(is.character(enzyme_menu_choices(menu_data())[["All other enzymes"]]))
})

# ---- the table the menu is built from ------------------------------------------------------------------

test_that("the Common column in the enzyme table matches the fixture list of common enzymes", {
  d <- menu_data()
  expect_setequal(d$Enzyme[d$Common], common_enzymes)
  expect_type(d$Common, "logical")
  expect_false(anyNA(d$Common))
})

# ---- the menu in the real app, in a real browser -----------------------------------------------------------
# Starts the app and opens it in headless Chrome, so this checks what a user actually sees. It is skipped
# when Chrome or the packages that drive it are not available, or when the tests are not run from the
# source tree.

start_app_in_background <- function() {
  root <- normalizePath(testthat::test_path("..", ".."))
  skip_if(!file.exists(file.path(root, "DESCRIPTION")), "not running from the source tree")
  port <- httpuv::randomPort()
  code <- sprintf("pkgload::load_all('%s', quiet = TRUE); run_app(options = list(port = %d, host = '127.0.0.1', launch.browser = FALSE))", root, port)
  app <- processx::process$new(file.path(R.home("bin"), "Rscript"), c("-e", code), stdout = tempfile(), stderr = tempfile())
  withr::defer(app$kill(), envir = parent.frame())
  ready <- FALSE
  for (i in 1:120) {
    ready <- tryCatch({ con <- socketConnection("127.0.0.1", port, blocking = TRUE, timeout = 1); close(con); TRUE },
                      error = function(e) FALSE, warning = function(w) FALSE)
    if (ready) break
    Sys.sleep(0.5)
  }
  skip_if(!ready, "the app did not start in time")
  port
}

test_that("in the real app the enzyme menu is grouped, empty at first, and has no 'select all'", {
  skip_on_cran()
  skip_if_not_installed("chromote")
  skip_if_not_installed("processx")
  skip_if_not_installed("pkgload")
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NA_character_)
  skip_if(is.na(chrome) || !nzchar(chrome), "Chrome is not available")

  # Chrome's own sandbox does not start inside containers or on CI runners that restrict user namespaces
  # (Ubuntu 24.04). This test only opens the app on 127.0.0.1, so it runs Chrome without the sandbox.
  chromote::set_chrome_args(c(chromote::default_chrome_args(), "--no-sandbox"))
  port <- start_app_in_background()
  b <- tryCatch(chromote::ChromoteSession$new(), error = function(e) NULL)
  skip_if(is.null(b), "Chrome could not be started here")
  withr::defer(b$close())
  b$Page$navigate(sprintf("http://127.0.0.1:%d", port))
  b$Page$loadEventFired()
  js <- function(expr) b$Runtime$evaluate(expr, returnByValue = TRUE)$result$value
  select <- "document.getElementById('Oenzymes').selectize"

  # Wait for the server to fill the menu.
  n <- -1
  for (i in 1:80) {
    n <- js("(function(){var e = document.getElementById('Oenzymes'); return e && e.selectize ? Object.keys(e.selectize.options).length : -1})()")
    if (isTRUE(n > 200)) break
    Sys.sleep(0.5)
  }
  expect_equal(n, nrow(menu_data()))

  groups <- jsonlite::fromJSON(js(sprintf("JSON.stringify(Object.keys(%s.optgroups).map(function(k){return %s.optgroups[k].label}))", select, select)))
  expect_identical(groups, c("Commonly used", "Golden Gate", "All other enzymes"))

  sizes <- jsonlite::fromJSON(js(sprintf("(function(){var c = {}; Object.values(%s.options).forEach(function(o){[].concat(o.optgroup).forEach(function(g){c[g] = (c[g] || 0) + 1})}); return JSON.stringify(c)})()", select)))
  expect_equal(unlist(sizes[groups]), c("Commonly used" = 22, "Golden Gate" = 8, "All other enzymes" = nrow(menu_data()) - 30))

  # Every enzyme is in exactly one group, and BsaI and BsmBI are under Golden Gate only.
  bsai <- jsonlite::fromJSON(js(sprintf("JSON.stringify([].concat(%s.options['BsaI'].optgroup))", select)))
  bsmbi <- jsonlite::fromJSON(js(sprintf("JSON.stringify([].concat(%s.options['BsmBI'].optgroup))", select)))
  expect_identical(bsai, "Golden Gate")
  expect_identical(bsmbi, "Golden Gate")

  # Nothing is pre-selected.
  expect_equal(length(jsonlite::fromJSON(js(sprintf("JSON.stringify(%s.items)", select)))), 0)

  # Open the menu (the sidebar panel has to be expanded and the checkbox ticked first).
  js("Array.from(document.querySelectorAll('.accordion-button')).filter(function(x){return /Biological Shields/.test(x.innerText)})[0].click()")
  Sys.sleep(1)
  js("document.getElementById('checkEnzySites').click()")
  Sys.sleep(1)
  js(sprintf("%s.focus(); %s.open()", select, select))
  Sys.sleep(1)
  headers <- jsonlite::fromJSON(js("JSON.stringify(Array.from(document.querySelectorAll('#Oenzymes ~ .selectize-control .optgroup-header')).map(function(e){return e.innerText}))"))
  expect_identical(headers, c("Commonly used", "Golden Gate", "All other enzymes"))

  # There is no way to select everything at once.
  menu_text <- js("document.querySelector('#Oenzymes ~ .selectize-control').parentElement.innerText")
  expect_false(grepl("select\\s*all", menu_text, ignore.case = TRUE))
})

test_that("choosing no enzyme while the box is ticked runs without error and removes nothing", {
  # The menu starts empty, so a user can tick "Remove Restriction Sites" and choose nothing.
  cds <- "ATGCTCGAGGGCGCCTAA"
  d <- run_server_downloads("100", cds, extra = list(checkEnzySites = TRUE, Oenzymes = NULL), enzymes = NULL)
  expect_identical(d$viewer_seq, cds)                       # method "No optimization": nothing changed
  comments <- parse_genbank(d$opt)$comments
  expect_true(any(grepl("Remove RE sites: No", comments)))   # the log does not claim enzymes were removed
  expect_false(any(grepl("Restriction sites removed", comments)))
})
