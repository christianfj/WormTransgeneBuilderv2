#' Safely Load CSV Data
#'
#' Searches the official package extdata folder to load UI metadata on the cloud.
#'
#' @param filename Character; the name of the CSV file to load.
#' @return A data frame if found, otherwise NULL.
#' @noRd
safe_load_csv <- function(filename) {
  # 1. The official way to find files inside a compiled R package on Render
  cloud_path <- system.file("extdata", filename, package = "WormTransgeneBuilder")

  if (cloud_path != "" && file.exists(cloud_path)) {
    return(read.csv(cloud_path, stringsAsFactors = FALSE))
  }

  # 2. Fallback paths for when you are testing locally on your Mac
  local_paths <- c(
    paste0("inst/extdata/", filename),
    paste0("data-raw/", filename),
    paste0("data/", filename)
  )
  for (p in local_paths) {
    if (file.exists(p)) return(read.csv(p, stringsAsFactors = FALSE))
  }

  return(NULL)
}

#' Find Sequence Matches
#'
#' Locates exact or fuzzy matches of a pattern within a target sequence,
#' computing coordinates and formatting them for the sequence visualizer.
#'
#' @param pattern Character; the DNA string to search for.
#' @param seq_to_search Character; the target DNA sequence.
#' @param color Character; hex color code for annotation mapping.
#' @param tooltip Character; text to display on hover in the viewer.
#' @param feat_type Character; type of feature being annotated.
#' @param is_enzyme Logical; if TRUE, also searches for the reverse complement.
#' @param max_mismatch Integer; maximum allowed mismatches.
#' @return A list containing starts, ends, colors, tooltips, and feature types.
#' @noRd
add_matches <- function(pattern, seq_to_search, color, tooltip, feat_type, is_enzyme=FALSE, max_mismatch=0) {
  res <- list(starts=c(), ends=c(), colors=c(), tooltips=c(), feats=c())

  if (is_enzyme) {
    # Restriction sites may contain IUPAC ambiguity codes (e.g. BanI GGYRCC), so
    # they are matched with the same helper that the site-removal step uses.
    starts_fwd <- find_enzyme_site_starts(toupper(seq_to_search), toupper(pattern))
    if (length(starts_fwd) > 0) {
      res$starts <- as.numeric(starts_fwd)
      res$ends <- res$starts + nchar(pattern) - 1
    }
  } else if (max_mismatch == 0) {
    locs <- gregexpr(toupper(pattern), toupper(seq_to_search), fixed=TRUE)[[1]]
    if (locs[1] != -1) {
      res$starts <- as.numeric(locs)
      res$ends <- res$starts + nchar(pattern) - 1
    }
  } else {
    # Windows of the sequence that differ from the pattern at no more than max_mismatch places.
    # Only whole windows inside the sequence count.
    starts_approx <- approximate_match_starts(pattern, seq_to_search, max_mismatch)
    if (length(starts_approx) > 0) {
      res$starts <- starts_approx
      res$ends <- starts_approx + nchar(pattern) - 1L
    }
  }

  if (length(res$starts) > 0) {
    res$colors <- rep(color, length(res$starts))
    res$tooltips <- rep(tooltip, length(res$starts))
    res$feats <- rep(feat_type, length(res$starts))
  }

  if (is_enzyme) {
    rev_pat <- reverse_complement(toupper(pattern))
    if (rev_pat != toupper(pattern)) {
      starts_rev <- as.numeric(find_enzyme_site_starts(toupper(seq_to_search), rev_pat))
      if (length(starts_rev) > 0) {
        ends_rev <- starts_rev + nchar(pattern) - 1
        res$starts <- c(res$starts, starts_rev)
        res$ends <- c(res$ends, ends_rev)
        res$colors <- c(res$colors, rep(color, length(starts_rev)))
        res$tooltips <- c(res$tooltips, rep(paste(tooltip, "(rev)"), length(starts_rev)))
        res$feats <- c(res$feats, rep(feat_type, length(starts_rev)))
      }
    }
  }
  return(res)
}

# Pre-load datasets
ui_data_introns   <- safe_load_csv("Introns.csv")
ui_data_5putrs    <- safe_load_csv("5pUTRs.csv")
ui_data_3putrs    <- safe_load_csv("3pUTRs.csv")
ui_data_promoters <- safe_load_csv("Promoters.csv")
ui_data_pirnas    <- safe_load_csv("piRNA_Database.csv")
ui_data_enzymes   <- safe_load_csv("Enzymes.csv")
ui_data_tags      <- safe_load_csv("Tags.csv")

#' The application server-side
#'
#' Orchestrates data flow, reactive inputs, optimization calls, and dynamic UI rendering.
#'
#' @param input,output,session Internal parameters for `{shiny}`.
#' @noRd
app_server <- function(input, output, session) {

  # Shown in the About tab, so any result can be tied to the software that made it.
  output$software_versions <- shiny::renderText({
    v <- software_versions()
    sprintf("WormTransgeneBuilder %s, R %s, ViennaRNA %s", v[["app"]], v[["R"]], v[["ViennaRNA"]])
  })

  # Character counter logic
  output$char_count <- shiny::renderText({
    is_prot <- (input$intypeinput == 2)
    current_seq <- if (!is_prot) input$seqDNA else input$seqPROT
    paste("Characters:", nchar(current_seq), "/ 20,000")
  })

  shiny::observe({
    if (!is.null(ui_data_introns)) {
      clean_labels <- gsub("<[^>]+>", "", ui_data_introns$UI_HTML)
      choice_vals <- c(ui_data_introns$Name, "Synthetic")
      choice_names <- c(clean_labels, "Synthetic introns (random from 1K)")
      shiny::updateSelectInput(session, "intropt", choices = setNames(choice_vals, choice_names))
    }
    if (!is.null(ui_data_tags)) {
      clean_labels_tags <- gsub("<[^>]+>", "", ui_data_tags$Name)
      shiny::updateSelectizeInput(session, "n_tags", choices = setNames(ui_data_tags$Name, clean_labels_tags))
      shiny::updateSelectizeInput(session, "c_tags", choices = setNames(ui_data_tags$Name, clean_labels_tags))
    }
    if (!is.null(ui_data_enzymes)) {
      # Grouped choices; nothing is pre-selected.
      shiny::updateSelectizeInput(session, "Oenzymes", choices = enzyme_menu_choices(ui_data_enzymes))
    }
    if (!is.null(ui_data_5putrs))    shiny::updateSelectInput(session, "p5UTR", choices = setNames(ui_data_5putrs$Name, ui_data_5putrs$Name))
    if (!is.null(ui_data_3putrs))    shiny::updateSelectInput(session, "p3UTR", choices = setNames(ui_data_3putrs$Name, ui_data_3putrs$Name))
    if (!is.null(ui_data_promoters)) shiny::updateSelectInput(session, "oppromo", choices = promoter_menu_choices(ui_data_promoters))
  })

  shiny::observe({
    if (input$checkboxRibo == TRUE) {
      shiny::updateCheckboxInput(session, "checkfouras", value = TRUE)
      shinyjs::disable("checkfouras")
    } else {
      shinyjs::enable("checkfouras")
    }
  })

  shiny::observeEvent(input$actionRESET, {
    shinyjs::reset("panels")
    output$AllResults <- shiny::renderUI({ NULL })
  })

  shiny::observeEvent(input$actionSeq, {
    set.seed(NULL)

    is_prot <- (input$intypeinput == 2)
    raw_seq <- if (!is_prot) input$seqDNA else input$seqPROT
    raw_seq <- gsub("[^A-Z*]", "", toupper(raw_seq))

    # Safety checks: Empty or too long
    if (nchar(raw_seq) == 0) {
      shiny::showNotification("Please enter a valid sequence.", type = "error")
      return()
    }
    if (nchar(raw_seq) > 20000) {
      shiny::showNotification("Sequence is too long! Please limit input to 20,000 characters (20 kb).", type = "error")
      return()
    }

    shinyjs::disable("actionSeq")
    shiny::showNotification("Optimizing Sequence... Please wait.", id = "job_loader", duration = NULL, type = "message")

    opt_method_choice <- switch(input$selectCAI,
                                "100" = "None",
                                "1"  = "Ubiq",
                                "2"  = "Chance",
                                "7"  = "CAI1",
                                "6"  = "GLO")

    # The "High expression" (CAI1) method makes no random choices (unless the ribosome-binding-site step is on, see below), so it needs no seed.
    # "No codon optimization" only makes random choices when it removes restriction or
    # piRNA sites; those are kept repeatable with a seed that applies to the
    # optimization call alone (see below). The seed must never be set globally: it
    # would also fix later random choices, such as which synthetic introns are picked,
    # and would affect every other user of this R process.
    # With the ribosome-binding-site step "High expression" also draws random candidates for the start,
    # so that run gets a fixed seed too and stays repeatable.
    seed_for_method <- if (opt_method_choice == "None" || (opt_method_choice == "CAI1" && isTRUE(input$checkboxRibo))) 12345 else NULL

    enzymes_to_remove <- if (input$checkEnzySites) input$Oenzymes else NULL

    tag_n_seqs <- c()
    tag_c_seqs <- c()
    tag_n_names <- c()
    tag_c_names <- c()

    if (input$checkTags && !is.null(ui_data_tags)) {
      if (!is.null(input$n_tags)) {
        for (t in input$n_tags) {
          row_idx <- which(ui_data_tags$Name == t)
          if (length(row_idx) > 0) {
            tag_n_seqs <- c(tag_n_seqs, ui_data_tags[row_idx[1], "Sequence"])
            tag_n_names <- c(tag_n_names, gsub("<[^>]+>", "", t))
          }
        }
      }
      if (!is.null(input$c_tags)) {
        for (t in input$c_tags) {
          row_idx <- which(ui_data_tags$Name == t)
          if (length(row_idx) > 0) {
            tag_c_seqs <- c(tag_c_seqs, ui_data_tags[row_idx[1], "Sequence"])
            tag_c_names <- c(tag_c_names, gsub("<[^>]+>", "", t))
          }
        }
      }
    }

    introns_to_pass <- NULL
    if (input$checkIntron) {
      if (input$intropt == "Synthetic") {
        local_synth <- safe_load_csv("Synthetic_Introns.csv")
        if (!is.null(local_synth) && "Sequence" %in% colnames(local_synth)) {
          n_req <- as.numeric(input$num_synth_introns)
          n_req <- min(n_req, nrow(local_synth))
          introns_to_pass <- sample(local_synth$Sequence, n_req, replace = FALSE)
        } else {
          shiny::showNotification("Could not find 'Synthetic_Introns.csv' or 'Sequence' column.", type = "warning")
        }
      } else {
        local_introns <- safe_load_csv("Introns.csv")
        if (!is.null(local_introns) && !is.null(input$intropt)) {
          target_intron <- input$intropt
          if (target_intron == "1") target_intron <- "Fire_lab_introns"
          if (target_intron == "2") target_intron <- "rps-0"
          row_idx <- which(local_introns$Name == target_intron)
          if (length(row_idx) > 0) {
            introns_to_pass <- c(local_introns[row_idx, "Sequence_1"],
                                 local_introns[row_idx, "Sequence_2"],
                                 local_introns[row_idx, "Sequence_3"])
          }
        }
      }

      if (!is.null(introns_to_pass)) {
        introns_to_pass <- gsub("[^ATGCatgc]", "", introns_to_pass)
        introns_to_pass <- introns_to_pass[nchar(introns_to_pass) > 0]
      }

      if (!is.null(introns_to_pass) && length(introns_to_pass) > 0) {
        dna_len <- if (is_prot) nchar(raw_seq) * 3 else nchar(raw_seq)
        avg_exon <- dna_len / (length(introns_to_pass) + 1)
        if (avg_exon < 150) {
          shiny::showNotification(paste0("Warning: Average distance between introns will be ~", round(avg_exon), " bp, which is shorter than the recommended 150 bp limit."), type = "warning", duration = 15)
        }
      }
    }

    promoter_to_pass <- NULL
    utr5_to_pass <- NULL
    utr3_to_pass <- NULL

    if (input$checkPromoters && !is.null(input$oppromo)) {
      local_promoters <- safe_load_csv("Promoters.csv")
      if (!is.null(local_promoters)) {
        # Group headings have an empty Sequence and can never be chosen.
        row_idx <- which(local_promoters$Name == input$oppromo & nzchar(local_promoters$Sequence))
        if (length(row_idx) > 0) promoter_to_pass <- local_promoters[row_idx[1], "Sequence"]
      }
    }

    if (input$checkUTRs) {
      local_5putrs <- safe_load_csv("5pUTRs.csv")
      if (!is.null(input$p5UTR) && !is.null(local_5putrs)) {
        row_idx <- which(local_5putrs$Name == input$p5UTR)
        if (length(row_idx) > 0) utr5_to_pass <- local_5putrs[row_idx, "Sequence"]
      }
      local_3putrs <- safe_load_csv("3pUTRs.csv")
      if (!is.null(input$p3UTR) && !is.null(local_3putrs)) {
        row_idx <- which(local_3putrs$Name == input$p3UTR)
        if (length(row_idx) > 0) utr3_to_pass <- local_3putrs[row_idx, "Sequence"]
      }
    }

    spacing_choice <- "early"
    if (!is.null(input$intdistop)) spacing_choice <- if(input$intdistop == "2") "equidistant" else "early"
    force_frame <- if (!is.null(input$checkintframe)) input$checkintframe else TRUE
    target_hamming <- if (!is.null(input$selectPiMM)) as.numeric(input$selectPiMM) else 3

    tryCatch({
      opt_results <- with_user_notifications(with_local_seed(seed_for_method, optimize_transgene(
        dna_seq = raw_seq,
        opt_method = opt_method_choice,
        enzymes = enzymes_to_remove,
        remove_pirnas = input$checkPirna,
        target_min_hamming = target_hamming,
        introns = introns_to_pass,
        spacing_type = spacing_choice,
        force_reading_frame = force_frame,
        add_consensus_start = input$checkfouras,
        rbs_opt = input$checkboxRibo,
        promoter_seq = promoter_to_pass,
        utr5_seq = utr5_to_pass,
        utr3_seq = utr3_to_pass,
        is_protein = is_prot,
        tag_n_seqs = tag_n_seqs,
        tag_c_seqs = tag_c_seqs
      )))

      final_seq <- opt_results$final_seq
      optimized_cds <- opt_results$optimized_cds

      cai_scorer_col <- if(opt_method_choice == "CAI1") "BringMansWeights" else "Ubiq"
      cai_orig <- if (!is_prot) calc_cai(raw_seq, usage_column = cai_scorer_col) else NA
      gc_orig <- if (!is_prot) calc_gc_content(raw_seq) else NA
      glo_orig <- if (!is_prot) calc_glo_score(raw_seq) else NA

      cai_new  <- calc_cai(optimized_cds, usage_column = cai_scorer_col)
      gc_new  <- calc_gc_content(optimized_cds)
      glo_new  <- calc_glo_score(optimized_cds)

      twist_ui <- NULL
      if (isTRUE(input$checkTwisty)) {
        twist_report <- validate_for_synthesis(final_seq)
        if (twist_report$pass) {
          twist_ui <- shiny::div(shiny::p(shiny::strong("Status:"), " Sequence passes basic Twist synthesis guidelines."))
        } else {
          warnings_html <- lapply(twist_report$warnings, function(w) shiny::tags$li(w))
          twist_ui <- shiny::div(shiny::p(shiny::strong("Status:"), " Potential synthesis issues detected:"), shiny::tags$ul(warnings_html))
        }
      }

      opt_starts <- c(); opt_ends <- c(); opt_colors <- c(); opt_tooltips <- c(); opt_feats <- c()
      ori_starts <- c(); ori_ends <- c(); ori_colors <- c(); ori_tooltips <- c(); ori_feats <- c()

      if (!is_prot && nchar(raw_seq) >= 3) {
        ori_starts <- c(1, nchar(raw_seq)-2)
        ori_ends   <- c(3, nchar(raw_seq))
        ori_colors <- c("#a9dfbf", "#e6b0aa")
        ori_tooltips <- c("Start", "Stop")
        ori_feats  <- c("Start", "Stop")
      }

      if (!is.null(opt_results$cds_start) && !is.null(opt_results$cds_end)) {
        opt_starts <- c(opt_results$cds_start, opt_results$cds_end - 2)
        opt_ends   <- c(opt_results$cds_start + 2, opt_results$cds_end)
        opt_colors <- c("#a9dfbf", "#e6b0aa")
        opt_tooltips <- c("Start", "Stop")
        opt_feats  <- c("Start", "Stop")
      }

      if (!is.null(opt_results$tag_n_starts) && length(opt_results$tag_n_starts) > 0) {
        for (i in seq_along(opt_results$tag_n_starts)) {
          opt_starts <- c(opt_starts, opt_results$tag_n_starts[i])
          opt_ends <- c(opt_ends, opt_results$tag_n_ends[i])
          opt_colors <- c(opt_colors, "#d2b4de")
          opt_tooltips <- c(opt_tooltips, tag_n_names[i])
          opt_feats <- c(opt_feats, "misc_feature")
        }
      }

      if (!is.null(opt_results$tag_c_starts) && length(opt_results$tag_c_starts) > 0) {
        for (i in seq_along(opt_results$tag_c_starts)) {
          opt_starts <- c(opt_starts, opt_results$tag_c_starts[i])
          opt_ends <- c(opt_ends, opt_results$tag_c_ends[i])
          opt_colors <- c(opt_colors, "#d2b4de")
          opt_tooltips <- c(opt_tooltips, tag_c_names[i])
          opt_feats <- c(opt_feats, "misc_feature")
        }
      }

      if (!is.null(introns_to_pass) && length(introns_to_pass) > 0) {
        for (i in seq_along(introns_to_pass)) {
          res <- add_matches(introns_to_pass[i], final_seq, "#f5cba7", "Intron", "intron")
          opt_starts <- c(opt_starts, res$starts); opt_ends <- c(opt_ends, res$ends)
          opt_colors <- c(opt_colors, res$colors); opt_tooltips <- c(opt_tooltips, res$tooltips); opt_feats <- c(opt_feats, res$feats)
        }
      }
      if (!is.null(utr5_to_pass) && nchar(utr5_to_pass) > 0 && !tolower(input$p5UTR) %in% c("none", "")) {
        res <- add_matches(utr5_to_pass, final_seq, "#a3e4d7", input$p5UTR, "5'UTR")
        opt_starts <- c(opt_starts, res$starts); opt_ends <- c(opt_ends, res$ends)
        opt_colors <- c(opt_colors, res$colors); opt_tooltips <- c(opt_tooltips, res$tooltips); opt_feats <- c(opt_feats, res$feats)
      }
      if (!is.null(utr3_to_pass) && nchar(utr3_to_pass) > 0 && !tolower(input$p3UTR) %in% c("none", "")) {
        res <- add_matches(utr3_to_pass, final_seq, "#a9cce3", input$p3UTR, "3'UTR")
        opt_starts <- c(opt_starts, res$starts); opt_ends <- c(opt_ends, res$ends)
        opt_colors <- c(opt_colors, res$colors); opt_tooltips <- c(opt_tooltips, res$tooltips); opt_feats <- c(opt_feats, res$feats)
      }
      if (!is.null(promoter_to_pass) && nchar(promoter_to_pass) > 0 && !tolower(input$oppromo) %in% c("none", "")) {
        res <- add_matches(promoter_to_pass, final_seq, "#fad7a0", input$oppromo, "promoter")
        opt_starts <- c(opt_starts, res$starts); opt_ends <- c(opt_ends, res$ends)
        opt_colors <- c(opt_colors, res$colors); opt_tooltips <- c(opt_tooltips, res$tooltips); opt_feats <- c(opt_feats, res$feats)
      }

      if (input$checkEnzySites && length(input$Oenzymes) > 0 && !is.null(ui_data_enzymes)) {
        for (enz in input$Oenzymes) {
          site <- ui_data_enzymes$Site[ui_data_enzymes$Enzyme == enz]
          if (length(site) > 0 && !is_prot) {
            res_ori <- add_matches(site, raw_seq, "#d7bde2", paste(enz, "site"), "RE_site", is_enzyme = TRUE)
            ori_starts <- c(ori_starts, res_ori$starts); ori_ends <- c(ori_ends, res_ori$ends)
            ori_colors <- c(ori_colors, res_ori$colors); ori_tooltips <- c(ori_tooltips, res_ori$tooltips); ori_feats <- c(ori_feats, res_ori$feats)
          }
        }
      }

      if (input$checkPirna) {
        if (is_prot) {
          output$pirnaTable <- shiny::renderTable({ data.frame(Status = "Protein input detected: Final sequence was screened and cleared of all piRNA sites.") }, width = "100%")
        } else {
          init_dists <- opt_results$initial_pirna_dists
          fin_dists <- opt_results$final_pirna_dists
          problematic_idx <- which(init_dists < target_hamming)

          if (!is.null(ui_data_pirnas) && length(problematic_idx) > 0) {
            seqs_to_show <- WormTransgeneBuilder::pirna_sequences[problematic_idx]
            pirna_df <- ui_data_pirnas[ui_data_pirnas$Target_20mer %in% seqs_to_show, c("Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment")]
            match_idx <- match(pirna_df$Target_20mer, WormTransgeneBuilder::pirna_sequences)

            orig_ham_num <- as.integer(init_dists[match_idx])
            new_ham_num <- as.integer(fin_dists[match_idx])

            status_vec <- sapply(new_ham_num, function(val) {
              if (val >= target_hamming) "Safely Recoded" else "Unremovable"
            })
            pirna_df$Status <- status_vec

            survivors <- pirna_df$Target_20mer[new_ham_num < target_hamming]
            for (surv in survivors) {
              res <- add_matches(surv, final_seq, "#f9e79f", "piRNA (Unremovable)", "piRNA", max_mismatch = max(0, target_hamming - 1))
              opt_starts <- c(opt_starts, res$starts); opt_ends <- c(opt_ends, res$ends)
              opt_colors <- c(opt_colors, res$colors); opt_tooltips <- c(opt_tooltips, res$tooltips); opt_feats <- c(opt_feats, res$feats)
            }

            orig_targets <- pirna_df$Target_20mer[orig_ham_num < target_hamming]
            for (orig in orig_targets) {
              res <- add_matches(orig, raw_seq, "#f9e79f", "piRNA Target", "piRNA", max_mismatch = max(0, target_hamming - 1))
              ori_starts <- c(ori_starts, res$starts); ori_ends <- c(ori_ends, res$ends)
              ori_colors <- c(ori_colors, res$colors); ori_tooltips <- c(ori_tooltips, res$tooltips); ori_feats <- c(ori_feats, res$feats)
            }

            pirna_df$Original_Hamming <- as.character(orig_ham_num)
            pirna_df$New_Hamming <- as.character(new_ham_num)
            pirna_df$Target_20mer <- NULL
            pirna_df <- pirna_df[, c("Name", "Full_21mer", "Original_Hamming", "New_Hamming", "Abundance", "Enrichment", "Status")]
            colnames(pirna_df) <- c("Name", "piRNA", "Hamming (original)", "Hamming (new)", "Abundance", "Enrichment", "Status")

            output$pirnaTable <- shiny::renderTable({ pirna_df }, striped = TRUE, hover = TRUE, width = "100%", align = 'l')
          } else {
            output$pirnaTable <- shiny::renderTable({ data.frame(Status = paste("No piRNA sites found with a Hamming distance lower than", target_hamming)) }, width = "100%")
          }
        }
      } else {
        output$pirnaTable <- shiny::renderTable({ data.frame(Status = "piRNA minimization was disabled for this run.") }, width = "100%")
      }

      # For the High-expression method, underline the codons that are not the best codon
      # (those changed to remove restriction sites or piRNA matches).
      underline <- if (opt_method_choice == "CAI1") {
        non_optimal_codon_ranges(optimized_cds, opt_results$cds_start,
                                 opt_results$intron_insertions, opt_results$intron_lengths)
      } else {
        data.frame(start = integer(0), end = integer(0))
      }
      output$opt_sequence_viewer <- shiny::renderUI({ generate_sequence_viewer("opt_viewer_div", final_seq, opt_starts, opt_ends, opt_colors, opt_tooltips, seq_title = "DNA Sequence", underline_starts = underline$start, underline_ends = underline$end) })
      output$ori_sequence_viewer <- shiny::renderUI({ generate_sequence_viewer("ori_viewer_div", raw_seq, ori_starts, ori_ends, ori_colors, ori_tooltips, seq_title = ifelse(is_prot, "Protein Sequence", "DNA Sequence")) })

      output$AllResults <- shiny::renderUI({
        name_str <- ifelse(input$nameinput == "", "Sequence", input$nameinput)

        opt_stats <- shiny::div(
          style = "margin-bottom: 20px; font-size: 1.1em;",
          shiny::strong("Codon Adaptation Index: "), round(cai_new, 2), shiny::br(),
          shiny::strong("GC content: "), paste0(round(gc_new * 100), "%"), shiny::br(),
          shiny::strong("GLO Score: "), format(round(glo_new), big.mark=",")
        )

        ori_stats <- shiny::div(
          style = "margin-bottom: 20px; font-size: 1.1em;",
          shiny::strong("Codon Adaptation Index: "), if(!is.na(cai_orig)) round(cai_orig, 2) else "N/A", shiny::br(),
          shiny::strong("GC content: "), if(!is.na(gc_orig)) paste0(round(gc_orig * 100), "%") else "N/A", shiny::br(),
          shiny::strong("GLO Score: "), if(!is.na(glo_orig)) format(round(glo_orig), big.mark=",") else "N/A"
        )

        shiny::div(
          shiny::h2(paste("Optimized", name_str)),
          opt_stats,
          shiny::uiOutput("opt_sequence_viewer"),
          shiny::br(),
          shiny::downloadButton("downloadOptGenbank", "Genbank with annotations", class = "btn-secondary"),
          shiny::br(), shiny::hr(), shiny::br(),

          shiny::h2(paste(ifelse(is_prot, "Input Protein", "Input"), name_str)),
          ori_stats,
          shiny::uiOutput("ori_sequence_viewer"),
          shiny::br(),
          if(!is_prot) shiny::downloadButton("downloadOriGenbank", "Genbank with annotations", class = "btn-secondary") else NULL,
          shiny::br(), shiny::hr(), shiny::br(),

          if (!is.null(twist_ui)) shiny::div(shiny::h3("Twist Bioscience Synthesis Check"), twist_ui, shiny::hr(), shiny::br()) else NULL,

          shiny::h3("piRNA Target Analysis"),
          shiny::tableOutput("pirnaTable"),
          shiny::br(), shiny::br()
        )
      })

      param_log <- c()
      param_log <- c(param_log, paste("Optimization routine:", opt_method_choice))
      param_log <- c(param_log, paste("RBS optimization:", ifelse(input$checkboxRibo, "Yes", "No")))

      param_log <- c(param_log, paste("Add molecular tags:", ifelse(input$checkTags, "Yes", "No")))
      if (input$checkTags) {
        if (length(tag_n_names) > 0) param_log <- c(param_log, paste("N-terminal Tags:", paste(tag_n_names, collapse=", ")))
        if (length(tag_c_names) > 0) param_log <- c(param_log, paste("C-terminal Tags:", paste(tag_c_names, collapse=", ")))
      }

      param_log <- c(param_log, paste("Remove piRNA sites:", ifelse(input$checkPirna, "Yes", "No")))
      if (input$checkPirna) param_log <- c(param_log, paste("Maximum number of mismatches to 20-mer piRNA binding region:", target_hamming))
      param_log <- c(param_log, paste("Remove RE sites:", ifelse(input$checkEnzySites && length(input$Oenzymes) > 0, "Yes", "No")))
      if (input$checkEnzySites && length(input$Oenzymes) > 0) param_log <- c(param_log, paste("Restriction sites removed:", paste(input$Oenzymes, collapse=", ")))
      param_log <- c(param_log, paste("Add introns:", ifelse(input$checkIntron, "Yes", "No")))
      if (input$checkIntron) {
        if (input$intropt == "Synthetic") {
          param_log <- c(param_log, paste("Introns added:", length(introns_to_pass), "Random Synthetic Introns"))
        } else {
          param_log <- c(param_log, paste("Introns added: 3 x", input$intropt))
        }
      }
      param_log <- c(param_log, paste("Add UTRs:", ifelse(input$checkUTRs, "Yes", "No")))
      if (input$checkUTRs) {
        param_log <- c(param_log, paste("5' UTR:", input$p5UTR))
        param_log <- c(param_log, paste("3' UTR:", input$p3UTR))
      }
      param_log <- c(param_log, paste("Add promoter:", ifelse(input$checkPromoters, "Yes", "No")))
      if (input$checkPromoters) param_log <- c(param_log, paste("Promoter:", input$oppromo))
      param_log <- c(param_log, paste("Add consensus start site:", ifelse(input$checkfouras, "Yes", "No")))

      output$downloadOptGenbank <- shiny::downloadHandler(
        filename = function() { paste0("Optimized_", ifelse(input$nameinput == "", "Sequence", input$nameinput), ".gb") },
        content = function(file) {
          gb_lines <- export_genbank(
            locus_name = ifelse(input$nameinput == "", "Sequence", input$nameinput), sequence = final_seq,
            starts = opt_starts, ends = opt_ends, colors = opt_colors, tooltips = opt_tooltips, feat_types = opt_feats,
            cai_score = round(cai_new, 2), gc_score = paste0(round(gc_new * 100), "%"), glo_score = format(round(glo_new), big.mark=","),
            input_seq = raw_seq, extracomments = param_log
          )
          writeLines(gb_lines, file)
        }
      )
      output$downloadOriGenbank <- shiny::downloadHandler(
        filename = function() { paste0("Original_", ifelse(input$nameinput == "", "Sequence", input$nameinput), ".gb") },
        content = function(file) {
          gb_lines <- export_genbank(
            locus_name = ifelse(input$nameinput == "", "Sequence", input$nameinput), sequence = raw_seq,
            starts = ori_starts, ends = ori_ends, colors = ori_colors, tooltips = ori_tooltips, feat_types = ori_feats,
            cai_score = round(cai_orig, 2), gc_score = paste0(round(gc_orig * 100), "%"), glo_score = format(round(glo_orig), big.mark=","),
            extracomments = c("Original unmodified sequence")
          )
          writeLines(gb_lines, file)
        }
      )

      shiny::removeNotification("job_loader")
      shinyjs::enable("actionSeq")

    }, error = function(e) {
      shiny::removeNotification("job_loader")
      shinyjs::enable("actionSeq")
      shiny::showNotification(paste("Optimization Failed:", e$message), type = "error", duration = 10)
    })
  })
}
