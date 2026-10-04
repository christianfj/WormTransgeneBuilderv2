#' GenBank Date
#'
#' Formats a date the way the GenBank LOCUS line expects it: a two-digit day, a
#' three-letter English month and the year (for example `09-Sep-2026`). The month
#' comes from `month.abb`, so it does not depend on the computer's language settings.
#'
#' @param date A `Date`; defaults to today.
#' @return A character string.
#' @noRd
genbank_date <- function(date = Sys.Date()) {
  paste(format(date, "%d"), month.abb[as.integer(format(date, "%m"))], format(date, "%Y"), sep = "-")
}

#' Generate GenBank (.gb) File for Export
#'
#' Compiles a sequence, its optimization metrics, and all mapped features into a
#' standard GenBank file format. It injects specific `/ApEinfo` tags so that downstream
#' software like "A plasmid Editor" (ApE) correctly parses colors and tooltips.
#'
#' @param locus_name Character string; the identifier used in the LOCUS field.
#' @param sequence Character string; the core DNA sequence to export.
#' @param starts Numeric vector; 1-based start coordinates of features.
#' @param ends Numeric vector; 1-based end coordinates of features.
#' @param colors Character vector; hex color codes assigned to features.
#' @param tooltips Character vector; textual labels for the features.
#' @param feat_types Character vector; GenBank feature keys (e.g., "misc_feature", "piRNA").
#' @param cai_score Numeric; Codon Adaptation Index score to append to COMMENTS.
#' @param gc_score Character string; Formatted GC content percentage to append.
#' @param glo_score Character string; Formatted GLO score to append.
#' @param input_seq Character string; The original, un-optimized input sequence (optional).
#' @param extracomments Character vector; Additional workflow parameters to log in COMMENTS.
#' @return A character vector where each element is a line of the GenBank file.
#' @export
export_genbank <- function(locus_name, sequence, starts=c(), ends=c(), colors=c(), tooltips=c(), feat_types=c(),
                           cai_score=NULL, gc_score=NULL, glo_score=NULL, input_seq=NULL, extracomments=c()) {

  if (is.null(sequence) || !is.character(sequence) || nchar(sequence) == 0) return(c())

  n_features <- length(starts)
  if (!all(c(length(ends), length(colors), length(tooltips), length(feat_types)) == n_features)) {
    stop("starts, ends, colors, tooltips and feat_types must all have the same length.")
  }

  FileLines <- c()
  FileLines <- append(FileLines, paste("LOCUS", paste(locus_name, sep="", collapse=""), paste(nchar(sequence),"bp ds-DNA", sep=""), "linear", genbank_date(), sep="     "))
  FileLines <- append(FileLines, c("DEFINITION     .", "ACCESSION     .", "VERSION     .", "SOURCE     .", "ORGANISM     C.elegans"))
  FileLines <- append(FileLines, paste("COMMENT", paste(locus_name, " (Coding sequence is annoted in uppercase)"), sep="     "))

  if (!is.null(cai_score)) FileLines <- append(FileLines, paste("COMMENT     Codon Adaptation Index", cai_score))
  if (!is.null(gc_score)) FileLines <- append(FileLines, paste("COMMENT     GC", gc_score))
  if (!is.null(glo_score)) FileLines <- append(FileLines, paste("COMMENT     GLO Score", glo_score))
  if (!is.null(input_seq)) FileLines <- append(FileLines, paste("COMMENT     Input sequence:", input_seq))

  if (length(extracomments) > 0) {
    FileLines <- append(FileLines, paste("COMMENT    ", extracomments))
  }

  FileLines <- append(FileLines, "FEATURES             Location/Qualifiers")

  if (length(starts) > 0) {
    for (i in seq_along(starts)) {
      start_pos <- starts[i]
      end_pos <- ends[i]
      ftype <- feat_types[i]

      if (ftype == "piRNA") FileLines <- append(FileLines, paste0("     ", ftype, "           complement(", start_pos, "..", end_pos, ")"))
      else FileLines <- append(FileLines, paste0("     ", ftype, "           ", start_pos, "..", end_pos))

      FileLines <- append(FileLines, paste0("                     /locus_tag=\"", tooltips[i], "\""))
      FileLines <- append(FileLines, paste0("                     /ApEinfo_label=\"", tooltips[i], "\""))
      FileLines <- append(FileLines, paste0("                     /ApEinfo_fwdcolor=\"", colors[i], "\""))
      FileLines <- append(FileLines, paste0("                     /ApEinfo_revcolor=\"", colors[i], "\""))
      FileLines <- append(FileLines, "                     /ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"")
    }
  }

  FileLines <- append(FileLines, "ORIGIN")
  Compseq <- unlist(strsplit(sequence, ""))
  partseq <- c()

  for (i in seq(1, length(Compseq), 10)) {
    endseq <- min(i + 9, length(Compseq))
    partseq <- append(partseq, paste(Compseq[i:endseq], collapse=""))
  }

  idx <- 1
  for (num in seq(1, length(Compseq), 60)) {
    endpart <- min(idx + 5, length(partseq))
    FileLines <- append(FileLines, paste(sprintf("%9s", num), paste(partseq[idx:endpart], collapse=" ")))
    idx <- idx + 6
  }

  FileLines <- append(FileLines, "//")
  return(FileLines)
}
