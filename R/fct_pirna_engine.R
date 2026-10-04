#' Fast Exact Hamming Distances Between 20-mers and the piRNA Database
#'
#' Every function that screens a sequence for piRNA matches needs the Hamming distance
#' (the number of differing bases) between each 20-base window of the sequence and each
#' of the ~18,000 piRNAs. Comparing them character by character is the slowest step of
#' the pipeline, so this does it with whole-number arithmetic instead, giving exactly the
#' same distances.
#'
#' Each base is stored in 2 bits (A = 0, C = 1, G = 2, T = 3), so a 20-mer becomes two
#' whole numbers holding 10 bases (20 bits) each. The bitwise XOR of two such numbers is
#' non-zero in exactly the base positions where the 20-mers differ, and a lookup table
#' (`mismatch_count_table()`) says how many of the 10 positions those are. The distance
#' is the count for the first half plus the count for the second half. One window is
#' compared with every piRNA at once, as whole vectors.
#'
#' A window that is not plain A/C/G/T (for example one containing `N`) cannot be stored
#' this way, so it is compared with `stringdist` exactly as before. So is a piRNA that is
#' 20 long but contains `N`.
#'
#' Long sequences are split into blocks of consecutive windows that are scanned on
#' several CPU cores at once, and the pieces are put back together in a fixed order, so
#' the result does not depend on the number of cores used.
#' @noRd
NULL

# For every 20-bit code (10 bases), how many of the 10 two-bit groups are non-zero.
mismatch_count_table <- local({
  tbl <- NULL
  function() {
    if (is.null(tbl)) {
      x <- 0:(4^10 - 1)
      n <- integer(length(x))
      for (g in 0:9) n <- n + (((x %/% 4^g) %% 4) != 0)
      tbl <<- n
    }
    tbl
  }
})

# Split 20-mers into two 10-base whole numbers. `valid` marks the ones that are plain
# A/C/G/T and exactly 20 long; the numbers for the others are NA.
encode_20mers <- function(seqs) {
  valid <- !is.na(seqs) & grepl("^[ACGT]{20}$", seqs)
  codes <- chartr("ACGT", "0123", seqs)
  lo <- strtoi(substr(codes, 1, 10), base = 4L)
  hi <- strtoi(substr(codes, 11, 20), base = 4L)
  valid <- valid & !is.na(lo) & !is.na(hi)
  lo[!valid] <- NA_integer_
  hi[!valid] <- NA_integer_
  list(lo = lo, hi = hi, valid = valid)
}

#' Prepare a Set of piRNA Sequences for Fast Comparison
#'
#' @param seqs Character vector of piRNA sequences (may contain `NA` or odd entries).
#' @return A list with `source` (the sequences as given), their encoded halves `lo` and
#'   `hi`, `valid` (which ones could be encoded), `odd` (exactly 20 long but not plain
#'   A/C/G/T, for example containing `N`: compared the slow, general way) and `unusable`
#'   (missing or the wrong length: they cannot be compared at all).
#' @noRd
make_pirna_reference <- function(seqs) {
  enc <- encode_20mers(seqs)
  odd <- which(!enc$valid & !is.na(seqs) & nchar(seqs) == 20)
  list(source = seqs, lo = enc$lo, hi = enc$hi, valid = enc$valid,
       odd = odd, unusable = setdiff(which(!enc$valid), odd))
}

# The full piRNA database, encoded once and kept for later calls. It is encoded again
# only if the dataset itself has changed.
pirna_reference <- local({
  cache <- NULL
  function() {
    p <- WormTransgeneBuilder::pirna_sequences
    if (is.null(cache) || !identical(cache$source, p)) cache <<- make_pirna_reference(p)
    cache
  }
})

#' How Many CPU Cores to Scan With
#'
#' Small scans are done on one core, because starting the workers would take longer than
#' the scan. Larger ones use all cores but one (as the previous implementation did), or the
#' number set in `options(wormtransgenebuilder.cores = )`. Windows cannot fork, so it
#' always uses one.
#'
#' @param n_windows Number of 20-base windows to scan.
#' @return An integer, at least 1.
#' @noRd
pirna_scan_cores <- function(n_windows) {
  if (.Platform$OS.type == "windows") return(1L)
  if (n_windows < 300) return(1L)  # checked first, so small scans never ask about cores
  cores <- getOption("wormtransgenebuilder.cores", available_cores())
  as.integer(max(1L, min(cores, n_windows %/% 150)))
}

#' Number of Cores Available for Scanning
#'
#' All cores but one, as `parallel::detectCores()` reports them. That call asks the operating
#' system and takes several milliseconds, and the piRNA shield scans about 160 times per run,
#' so the answer is looked up once and kept.
#'
#' @param refresh If `TRUE`, ask again.
#' @return An integer, at least 1.
#' @noRd
available_cores <- local({
  n <- NULL
  function(refresh = FALSE) {
    if (refresh || is.null(n)) n <<- max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
    n
  }
})

# Scan one block of windows against the piRNAs. Returns the smallest distance per piRNA
# and per window (as integers, with 1000 standing for "could not be compared") and, if a
# threshold is given, the piRNAs closer than it for each window. Window numbers are
# relative to the block.
scan_windows <- function(kmers, enc, reference, threshold) {
  no_match <- 1000L
  tbl <- mismatch_count_table()
  n_ref <- length(reference$lo)
  n_win <- length(kmers)
  unusable <- reference$unusable
  odd <- reference$odd
  col_min <- rep(no_match, n_ref)
  row_min <- integer(n_win)
  pair_window <- vector("list", n_win)
  pair_pirna <- vector("list", n_win)

  for (w in seq_len(n_win)) {
    if (enc$valid[w]) {
      d <- tbl[bitwXor(reference$lo, enc$lo[w]) + 1L] + tbl[bitwXor(reference$hi, enc$hi[w]) + 1L]
      if (length(odd) > 0) {
        d[odd] <- suppressWarnings(as.integer(stringdist::stringdist(kmers[w], reference$source[odd], method = "hamming")))
      }
    } else {
      # Not plain A/C/G/T: compare the slow, general way. NA and Inf (different length)
      # both mean "cannot be compared".
      d <- suppressWarnings(as.integer(stringdist::stringdist(kmers[w], reference$source, method = "hamming")))
    }
    d[is.na(d)] <- no_match
    if (length(unusable) > 0) d[unusable] <- no_match

    col_min <- pmin.int(col_min, d)
    row_min[w] <- min(d)
    if (!is.null(threshold)) {
      hit <- which(d < threshold)
      if (length(hit) > 0) {
        pair_window[[w]] <- rep.int(w, length(hit))
        pair_pirna[[w]] <- hit
      }
    }
  }
  list(col_min = col_min, row_min = row_min, pair_window = pair_window, pair_pirna = pair_pirna)
}

#' Compare 20-mers with the piRNAs
#'
#' @param kmers Character vector of 20-base windows (upper case).
#' @param reference Encoded piRNAs from `make_pirna_reference()`; defaults to the full database.
#' @param threshold If given, also return every (window, piRNA) pair closer than this.
#' @param cores Number of CPU cores to scan with; by default chosen from the number of windows.
#' @return A list with
#'   * `col_min`: for each piRNA, the smallest distance to any window (`NA` if it could not be
#'     compared), the same as the result of `get_pirna_distances()`;
#'   * `row_min`: for each window, the smallest distance to any piRNA (`Inf` if none could be
#'     compared);
#'   * `pairs`: only if `threshold` is given; a data frame of `window` and `pirna` indices
#'     closer than the threshold, ordered by piRNA and then by window.
#' @noRd
scan_kmers_against_pirnas <- function(kmers, reference = pirna_reference(), threshold = NULL,
                                      cores = pirna_scan_cores(length(kmers))) {
  n_ref <- length(reference$lo)
  n_win <- length(kmers)
  no_match <- 1000L
  empty_pairs <- data.frame(window = integer(0), pirna = integer(0))

  if (n_win == 0 || n_ref == 0) {
    return(list(col_min = rep(NA_real_, n_ref), row_min = rep(Inf, n_win),
                pairs = if (is.null(threshold)) NULL else empty_pairs))
  }

  enc <- encode_20mers(kmers)
  cores <- max(1L, min(as.integer(cores), n_win))
  blocks <- if (cores > 1) split(seq_len(n_win), cut(seq_len(n_win), cores, labels = FALSE)) else list(seq_len(n_win))
  scan_block <- function(i) {
    scan_windows(kmers[i], list(lo = enc$lo[i], hi = enc$hi[i], valid = enc$valid[i]), reference, threshold)
  }

  parts <- if (length(blocks) > 1) parallel::mclapply(blocks, scan_block, mc.cores = length(blocks)) else list(scan_block(blocks[[1]]))
  if (any(vapply(parts, function(p) !is.list(p) || is.null(p$col_min), TRUE))) {
    parts <- lapply(blocks, scan_block)  # a worker failed: do it again on this core
  }

  col_min <- Reduce(pmin.int, lapply(parts, `[[`, "col_min"))
  row_min <- unlist(lapply(parts, `[[`, "row_min"), use.names = FALSE)
  col_min <- as.numeric(col_min)
  col_min[col_min >= no_match] <- NA_real_
  row_min <- as.numeric(row_min)
  row_min[row_min >= no_match] <- Inf

  pairs <- NULL
  if (!is.null(threshold)) {
    offsets <- c(0L, cumsum(lengths(blocks))[-length(blocks)])
    window <- unlist(Map(function(p, off) unlist(p$pair_window, use.names = FALSE) + off, parts, offsets), use.names = FALSE)
    pirna <- unlist(lapply(parts, function(p) unlist(p$pair_pirna, use.names = FALSE)), use.names = FALSE)
    if (length(window) == 0) {
      pairs <- empty_pairs
    } else {
      # The same order as which(distance_matrix < threshold, arr.ind = TRUE): by piRNA
      # (column), then by window (row).
      o <- order(pirna, window)
      pairs <- data.frame(window = as.integer(window[o]), pirna = as.integer(pirna[o]))
    }
  }
  list(col_min = col_min, row_min = row_min, pairs = pairs)
}
