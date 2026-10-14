# Asks the C. elegans Codon Adapter web tool (Redemann et al. 2011; https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py) for the codon weights it
# uses, and compares them with inst/extdata/Codon_frequencies.csv. Run from the project root:
#   Rscript data-raw/codon_adapter/query_codon_adapter.R
#
# How it works: the page is an HTML form (method POST, fields name, sequence, introns, score, rnai_resistance, optimize_rbs, remove_splice_sites). The
# result page contains, in plain HTML, the original and optimized codons, and, in JavaScript arrays (arrSeriesOrig / arrSeriesOpt / arrTicks)
# that draw its "Codon Weights Details" bar chart, the weight of every codon of the input. So one input that contains every codon once gives every
# codon's weight, and the optimized sequence at CAI 1.0 gives the best codon of each amino acid. In a browser the same data is visible with
# right click > "View page source" after pressing Optimize.
url <- "https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py"
codons <- c(outer(c("A", "C", "G", "T"), outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0), paste0))
stops <- c("TAA", "TAG", "TGA")
sense <- setdiff(codons, stops)
sense <- c("ATG", setdiff(sense, "ATG"))                                   # starts with ATG, as a coding sequence should
query <- function(dna, name) {
  tmp <- tempfile(fileext = ".html")
  status <- system2("curl", c("-sS", "-m", "120", "-o", tmp, "--data-urlencode", shQuote(paste0("name=", name)), "--data-urlencode", shQuote(paste0("sequence=", dna)),
                              "--data-urlencode", "introns=0", "--data-urlencode", "score=1.0", shQuote(url)))
  stopifnot(status == 0)
  readLines(tmp, warn = FALSE)
}
parse <- function(html, n_codons) {
  txt <- paste(html, collapse = "\n")
  orig_w <- as.numeric(sub(".*arrSeriesOrig.push\\( \\[[^,]*, ([0-9.]+)\\] \\);.*", "\\1", regmatches(txt, gregexpr("arrSeriesOrig.push\\( \\[[^]]*\\] \\);", txt))[[1]]))
  opt_w <- as.numeric(sub(".*arrSeriesOpt.push\\( \\[[^,]*, ([0-9.]+)\\] \\);.*", "\\1", regmatches(txt, gregexpr("arrSeriesOpt.push\\( \\[[^]]*\\] \\);", txt))[[1]]))
  strip <- function(x) trimws(gsub("&nbsp;", " ", gsub("<[^>]*>", "", x)))
  tokens <- function(label_pattern) { x <- unlist(strsplit(strip(grep(label_pattern, html, value = TRUE)), " +")); x[grepl("^[ACGT]{3}$", x)] }   # the result is printed in rows of 30 codons
  orig_c <- tokens("ORIGINAL&nbsp;&nbsp;&nbsp;&nbsp;[ACGT]{3}")
  opt_c <- tokens("OPTIMIZED&nbsp;&nbsp;&nbsp;(<font[^>]*>)?[ACGT]{3}")
  cai_lines <- html[grep("Codon Adaptation Index", html) + 1]
  list(orig_w = orig_w, opt_w = opt_w, orig_c = orig_c, opt_c = opt_c, cai = trimws(sub("&nbsp;.*", "", cai_lines)))
}
res <- list()
for (stop in stops) {
  dna <- paste0(paste(sense, collapse = ""), stop)
  p <- parse(query(dna, paste0("all_codons_", stop)), length(sense) + 1)
  cat(sprintf("stop %s: %d codons sent; weights returned %d; codons echoed %d; optimized codons %d; CAI shown: %s\n", stop, length(sense) + 1, length(p$orig_w), length(p$orig_c), length(p$opt_c), paste(p$cai, collapse = " | ")))
  res[[stop]] <- p
}
sent <- c(sense, "TAA")
p <- res[["TAA"]]
stopifnot(identical(p$orig_c, sent), length(p$orig_w) == length(sent))
tool <- data.frame(codon = c(sense, stops), weight = NA_real_, optimized = NA_character_, stringsAsFactors = FALSE)
tool$weight[match(sent, tool$codon)] <- p$orig_w
tool$optimized[match(sent, tool$codon)] <- p$opt_c
for (stop in c("TAG", "TGA")) { q <- res[[stop]]; tool$weight[tool$codon == stop] <- q$orig_w[length(q$orig_w)]; tool$optimized[tool$codon == stop] <- q$opt_c[length(q$opt_c)] }
utils::write.csv(tool, "data-raw/codon_adapter/codon_adapter_weights.csv", row.names = FALSE)
cat("saved data-raw/codon_adapter/codon_adapter_weights.csv\n")
