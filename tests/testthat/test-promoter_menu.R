# The promoter menu is grouped by the heading rows of Promoters.csv (rows with an empty Sequence). It shows
# "Name Gene", but the value it returns, and so the label used in the viewer and the GenBank file, is the
# Name alone. Heading rows can never be chosen.

promoters <- function() utils::read.csv(system.file("extdata", "Promoters.csv", package = "WormTransgeneBuilder"),
                                        stringsAsFactors = FALSE)

test_that("the promoter file has two groups and no heading is offered as a promoter", {
  d <- promoters()
  headings <- d$Name[d$Sequence == ""]
  expect_length(headings, 2)
  menu <- promoter_menu_choices(d)
  expect_identical(names(menu), headings)
  expect_false(any(headings %in% unlist(menu, use.names = FALSE)))
  expect_false(any(headings %in% names(unlist(menu))))
})

test_that("every promoter with a sequence is offered exactly once, in file order", {
  d <- promoters()
  real <- d$Name[d$Sequence != ""]
  expect_identical(unname(unlist(promoter_menu_choices(d))), real)
})

test_that("the menu shows Name and Gene, and returns the Name", {
  menu <- promoter_menu_choices(promoters())
  adf <- menu[[1]][["pADF (WBGene00005359, srh-142)"]]
  expect_identical(adf, "pADF")
  afd <- names(menu[[1]])[menu[[1]] == "pAFD"]
  expect_identical(afd, "pAFD (WBGene00001535, gcy-8)")
})

test_that("promoter names contain no gene text and are unique", {
  d <- promoters()
  real <- d[d$Sequence != "", ]
  expect_equal(anyDuplicated(real$Name), 0)
  expect_false(any(grepl("WBGene|\\(", real$Name)))
  expect_true(all(nchar(real$Sequence) %in% c(296:300)))
})

test_that("grouping rule on a small example, including rows before any heading and a missing Gene column", {
  d <- data.frame(Name = c("p1", "Group A", "p2", "p3", "Group B", "Group C", "p4"),
                  Gene = c("(g1)", "", "(g2)", "", "", "", "(g4)"),
                  Sequence = c("acgt", "", "acgt", "acgt", "", "", "acgt"))
  expect_identical(promoter_menu_choices(d),
                   list("p1 (g1)" = "p1", "Group A" = c("p2 (g2)" = "p2", "p3" = "p3"), "Group C" = c("p4 (g4)" = "p4")))
  expect_identical(names(promoter_menu_choices(d[, c("Name", "Sequence")])), c("p1", "Group A", "Group C"))
})

test_that("choosing a group heading name gives no promoter in the output", {
  cds <- paste0("ATG", strrep("GCTAAG", 10), "TAA")
  heading <- promoters()$Name[1]
  d <- run_server_downloads("100", cds, list(checkPromoters = TRUE, oppromo = heading), enzymes = NULL)
  p <- parse_genbank(d$opt)
  expect_false(any(p$labels == heading))
  expect_false(any(p$labels == "promoter"))
})
