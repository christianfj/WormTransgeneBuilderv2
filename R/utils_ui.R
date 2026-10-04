#' Grouped Choices for the Restriction-Enzyme Menu
#'
#' With 234 enzymes in the list, a plain alphabetical menu is hard to use, so the menu is split
#' into three groups, shown in this order:
#'
#' * **Commonly used**: the enzymes most often avoided in a transgene (`Common` is TRUE in
#'   `inst/extdata/Enzymes.csv`), except the Golden Gate ones;
#' * **Golden Gate**: the enzymes used for modular Golden Gate cloning (`GG` is TRUE), including
#'   BsaI and BsmBI even though they are also commonly avoided;
#' * **All other enzymes**: everything else.
#'
#' Every enzyme is in exactly one group, and each group is in alphabetical order. There is no
#' "select all": removing many enzymes forces many codon changes and lowers how far the sequence
#' can be optimized, so the menu does not encourage it.
#'
#' @param enzymes A data frame with an `Enzyme` column and, optionally, logical `GG` and `Common`
#'   columns (a missing column counts as FALSE).
#' @return A named list of character vectors, one per non-empty group, in the form
#'   `shiny::selectizeInput()` accepts for grouped choices.
#' @noRd
enzyme_menu_choices <- function(enzymes) {
  flag <- function(column) {
    if (column %in% names(enzymes)) as.logical(enzymes[[column]]) %in% TRUE else rep(FALSE, nrow(enzymes))
  }
  golden_gate <- flag("GG")
  common <- flag("Common") & !golden_gate
  other <- !golden_gate & !common

  alphabetical <- function(x) x[order(tolower(x))]
  groups <- list(
    "Commonly used" = alphabetical(enzymes$Enzyme[common]),
    "Golden Gate" = alphabetical(enzymes$Enzyme[golden_gate]),
    "All other enzymes" = alphabetical(enzymes$Enzyme[other])
  )
  groups[lengths(groups) > 0]
}

#' Grouped Choices for the Promoter Menu
#'
#' `Promoters.csv` has the columns Name, Gene and Sequence. A row with an empty Sequence is a group
#' heading (its text is in Name), not a promoter, and cannot be chosen: it becomes the title of a
#' group in the menu. Each promoter is shown as "Name Gene" (for example "pADF (WBGene00005359,
#' srh-142)"), but the value it returns is the Name alone, which is what the sequence viewer and the
#' GenBank file use as the label. Promoters keep the order of the file.
#'
#' @param promoters A data frame with `Name`, `Gene` and `Sequence` columns.
#' @return A named list of named character vectors (names are the shown text, values are the promoter
#'   names), in the form `shiny::selectInput()` accepts for grouped choices. Promoters that come before
#'   any heading are returned as a plain named vector in front of the groups.
#' @noRd
promoter_menu_choices <- function(promoters) {
  text <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; trimws(x) }
  name <- text(promoters$Name)
  gene <- if ("Gene" %in% names(promoters)) text(promoters$Gene) else rep("", length(name))
  sequence <- text(promoters$Sequence)

  is_heading <- sequence == "" & name != ""
  is_promoter <- sequence != "" & name != ""
  shown <- trimws(paste(name, gene))
  # Which heading each row falls under (0 = none before it).
  group <- cumsum(is_heading)

  choices <- list()
  before_headings <- is_promoter & group == 0
  if (any(before_headings)) choices <- c(choices, as.list(stats::setNames(name[before_headings], shown[before_headings])))
  for (g in which(is_heading)) {
    members <- is_promoter & group == group[g]
    if (any(members)) choices[[name[g]]] <- stats::setNames(name[members], shown[members])
  }
  choices
}
