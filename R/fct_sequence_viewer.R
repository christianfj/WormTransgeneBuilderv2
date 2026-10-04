#' Generate Sequence Viewer HTML/JS
#'
#' Constructs the HTML container and injects the corresponding JavaScript
#' required to render the calipho-sib sequence-viewer widget. It maps
#' annotation coordinates (starts/ends) to specific colors and tooltips.
#'
#' @param div_id Character string; the unique HTML ID for the target rendering `div`.
#' @param sequence Character string; the primary DNA or Protein sequence to display.
#' @param starts Numeric vector; the start coordinates for all annotations (1-based).
#' @param ends Numeric vector; the end coordinates for all annotations (1-based).
#' @param colors Character vector; hex color codes corresponding to each annotation.
#' @param tooltips Character vector; hover text strings for each annotation.
#' @param seq_title Character string; the title displayed above the viewer (e.g., "DNA Sequence").
#' @param underline_starts,underline_ends Numeric vectors; 1-based start and end coordinates of
#'   regions to underline (for example codons that are not the best codon). An underline is drawn
#'   on top of any colour annotation in the same place, so both remain visible.
#' @param underline_label Character string; legend and hover text for the underlined regions.
#' @return A `shiny::HTML` object containing the `div` and `<script>` tags.
#' @export
generate_sequence_viewer <- function(div_id, sequence, starts = c(), ends = c(), colors = c(), tooltips = c(), seq_title = "DNA Sequence",
                                     underline_starts = c(), underline_ends = c(), underline_label = "Non-optimal codon") {
  if (is.null(sequence) || sequence == "") return(NULL)

  sequence_length <- nchar(sequence)
  coverage_str <- ""
  legend_str <- ""

  if (length(starts) > 0 || length(underline_starts) > 0) {
    color_map <- rep(NA, sequence_length)
    tooltip_map <- rep(NA, sequence_length)
    underline_map <- rep(FALSE, sequence_length)

    for (i in seq_along(starts)) {
      for (k in starts[i]:ends[i]) {
        if (k > 0 && k <= sequence_length) {
          if (is.na(color_map[k])) {
            color_map[k] <- colors[i]
            tooltip_map[k] <- tooltips[i]
          } else {
            if (!grepl(tooltips[i], tooltip_map[k], fixed = TRUE)) {
              tooltip_map[k] <- paste0(tooltip_map[k], " + ", tooltips[i])
              color_map[k] <- "#ccd1d1"
            }
          }
        }
      }
    }

    # Underlines are stored separately from the colour so that an underlined base
    # keeps the colour and hover text of any annotation that covers it.
    for (i in seq_along(underline_starts)) {
      ks <- underline_starts[i]:underline_ends[i]
      underline_map[ks[ks > 0 & ks <= sequence_length]] <- TRUE
    }
    tooltip_map[underline_map & is.na(tooltip_map)] <- underline_label

    matches_df <- data.frame(start=integer(), end=integer(), color=character(), tooltip=character(), underline=logical())
    current_color <- NA
    current_tooltip <- NA
    current_underline <- FALSE
    block_start <- NA

    for (pos in 1:sequence_length) {
      col <- color_map[pos]
      tip <- tooltip_map[pos]
      ul <- underline_map[pos]

      if (!identical(col, current_color) || !identical(tip, current_tooltip) || !identical(ul, current_underline)) {
        if (!is.na(current_color) || current_underline) {
          matches_df <- rbind(matches_df, data.frame(
            start = block_start - 1,
            end = pos - 1,
            color = current_color,
            tooltip = current_tooltip,
            underline = current_underline
          ))
        }
        current_color <- col
        current_tooltip <- tip
        current_underline <- ul
        block_start <- pos
      }
    }

    if (!is.na(current_color) || current_underline) {
      matches_df <- rbind(matches_df, data.frame(
        start = block_start - 1,
        end = sequence_length,
        color = current_color,
        tooltip = current_tooltip,
        underline = current_underline
      ))
    }

    if (nrow(matches_df) > 0) {
      cov_lines <- sapply(1:nrow(matches_df), function(i) {
        bg <- if (is.na(matches_df$color[i])) "transparent" else matches_df$color[i]
        paste0("{start: ", matches_df$start[i], ", end: ", matches_df$end[i],
               ", color: \"black\", bgcolor: \"", bg,
               "\", underscore: ", tolower(matches_df$underline[i]),
               ", tooltip: \"", matches_df$tooltip[i], "\"}")
      })
      coverage_str <- paste(cov_lines, collapse = ",\n")
    }

    if ("#ccd1d1" %in% color_map) {
      colors <- c(colors, "#ccd1d1")
      tooltips <- c(tooltips, "Multiple Annotations")
    }

    legend_df <- data.frame(name = as.character(tooltips), color = as.character(colors), underscore = rep(FALSE, length(tooltips)))
    if (any(underline_map)) {
      legend_df <- rbind(legend_df, data.frame(name = underline_label, color = "#ffffff", underscore = TRUE))
    }
    unique_legends <- unique(legend_df)
    if (nrow(unique_legends) > 0) {
      leg_lines <- sapply(1:nrow(unique_legends), function(i) {
        paste0("{name: \"", unique_legends$name[i], "\", color: \"", unique_legends$color[i],
               "\", underscore: ", tolower(unique_legends$underscore[i]), "}")
      })
      legend_str <- paste(leg_lines, collapse = ",\n")
    }
  }

  html_code <- paste0(
    "<div id=\"", div_id, "\"></div>\n",
    "<script type=\"text/javascript\">\n",
    "  function renderSeq_", div_id, "() {\n",
    "    if (typeof Sequence === 'undefined') return;\n",
    "    window.seq_", div_id, " = new Sequence('", sequence, "');\n",
    "    window.seq_", div_id, ".render('#", div_id, "', {\n",
    "      'title': '", seq_title, "',\n",
    "      'showLineNumbers': true,\n",
    "      'wrapAminoAcids': true,\n",
    "      'charsPerLine': 100,\n",
    "      'toolbar': false,\n",
    "      'search': false,\n",
    "      'sequenceMaxHeight': '400px',\n",
    "      'badge': false\n",
    "    });\n",
    "    window.cov_", div_id, " = [", coverage_str, "];\n",
    "    window.leg_", div_id, " = [", legend_str, "];\n",
    "    if (window.cov_", div_id, ".length > 0) window.seq_", div_id, ".coverage(window.cov_", div_id, ");\n",
    "    if (window.leg_", div_id, ".length > 0) window.seq_", div_id, ".addLegend(window.leg_", div_id, ");\n",
    "    \n",
    "    document.getElementById('", div_id, "').addEventListener('mouseup', function() {\n",
    "      setTimeout(function() {\n",
    "        if (window.cov_", div_id, ".length > 0) window.seq_", div_id, ".coverage(window.cov_", div_id, ");\n",
    "      }, 50);\n",
    "    });\n",
    "  }\n",
    "  if (typeof Sequence === 'undefined') {\n",
    "    var script = document.createElement('script');\n",
    "    script.src = 'https://cdn.jsdelivr.net/gh/calipho-sib/sequence-viewer@master/dist/sequence-viewer.min.js';\n",
    "    script.onload = function() { setTimeout(renderSeq_", div_id, ", 150); };\n",
    "    document.head.appendChild(script);\n",
    "  } else {\n",
    "    setTimeout(renderSeq_", div_id, ", 150);\n",
    "  }\n",
    "</script>"
  )

  return(shiny::HTML(html_code))
}

#' Find Codons That Are Not the Best Codon
#'
#' For the "High expression" (CAI1) method every codon is recoded to the best codon
#' for its amino acid, so the only codons that are not best are those that had to be
#' changed afterwards (to remove a restriction site or a piRNA match). This finds
#' them and converts their positions to coordinates in the final sequence, which also
#' contains any promoter, UTRs, tags and introns. A codon counts as best when its
#' `CAI1` weight is 1; stop codons all have weight 1 there, so they are never reported.
#'
#' @param cds Character string; the optimized coding sequence, before introns were inserted.
#' @param cds_start Integer; position in the final sequence of the first base of `cds`.
#' @param intron_insertions Numeric vector; bases of `cds` after which each intron was inserted.
#' @param intron_lengths Numeric vector; the length of each inserted intron.
#' @param best_column Character string; the codon-table column that is 1 for best codons.
#' @return A data frame with 1-based `start` and `end` columns (final-sequence
#'   coordinates), one row per codon that is not the best codon.
#' @noRd
non_optimal_codon_ranges <- function(cds, cds_start = 1, intron_insertions = numeric(0),
                                     intron_lengths = numeric(0), best_column = "CAI1") {
  empty <- data.frame(start = integer(0), end = integer(0))
  cds <- toupper(cds)
  usable <- nchar(cds) - nchar(cds) %% 3
  if (usable < 3) return(empty)

  codon_starts <- seq(1, usable - 2, by = 3)
  codons <- substring(cds, codon_starts, codon_starts + 2)
  cf <- WormTransgeneBuilder::codon_freqs
  idx <- match(codons, toupper(cf$Codon))
  non_best <- !is.na(idx) & cf[[best_column]][idx] != 1
  if (!any(non_best)) return(empty)

  s <- codon_starts[non_best]
  # Every intron inserted before a base pushes it further along the final sequence.
  shift <- function(pos) vapply(pos, function(b) sum(intron_lengths[intron_insertions < b]), numeric(1))
  data.frame(start = s + shift(s) + cds_start - 1,
             end = s + 2 + shift(s + 2) + cds_start - 1)
}
