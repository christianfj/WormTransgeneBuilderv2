#' The application User-Interface
#'
#' Defines the UI for the WormTransgeneBuilder Shiny application using the bslib framework.
#'
#' @param request Internal parameter for `{shiny}`.
#' @import shiny
#' @import bslib
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    shinyjs::useShinyjs(),

    page_navbar(
      theme = bs_theme(version = 5, primary = "#2c3e50", secondary = "#18bc9c"),
      title = "Wormbuilder",
      id = "panels",

      nav_panel(
        title = "Sequence Adaptation",

        layout_sidebar(
          sidebar = sidebar(
            width = "25%",
            title = "Transgene Parameters",

            shinyWidgets::radioGroupButtons("intypeinput", "Input Format", choices = list("DNA" = 1, "Protein" = 2), justified = TRUE),
            textInput("nameinput", "Sequence Name (Optional)", placeholder = "e.g., Fire lab gfp"),

            conditionalPanel("input.intypeinput == 1",
                             shiny::tagAppendAttributes(
                               textAreaInput("seqDNA", "DNA Sequence (Max 20kb)", height = "200px", placeholder = "ATG..."),
                               maxlength = "20000"
                             )),
            conditionalPanel("input.intypeinput == 2",
                             shiny::tagAppendAttributes(
                               textAreaInput("seqPROT", "Protein Sequence (Max 20kb)", height = "200px", placeholder = "M..."),
                               maxlength = "20000"
                             )),

            # Character counter display
            span(textOutput("char_count"), style = "color: #7f8c8d; font-size: 0.85em; margin-bottom: 10px; display: block;"),

            accordion(
              accordion_panel("1. Optimization Method",
                              selectizeInput("selectCAI", NULL,
                                             choices = list("No codon optimization" = 100,
                                                            "Ubiquitous expression (Serizay)" = 1,
                                                            "Random recoding" = 2,
                                                            "High expression (Redemann et al., 2011)" = 7,
                                                            "Germline optimized (GLO)" = 6))),

              accordion_panel("2. Molecular Tags",
                              checkboxInput("checkTags", "Add Molecular Tags", value = FALSE),
                              conditionalPanel("input.checkTags == 1",
                                               selectizeInput("n_tags", "N-terminal Tags (added after ATG)", choices = c("Loading..."), multiple = TRUE),
                                               selectizeInput("c_tags", "C-terminal Tags (added before STOP)", choices = c("Loading..."), multiple = TRUE))),

              accordion_panel("3. Biological Shields",
                              checkboxInput("checkPirna", "Minimize piRNA homology", value = FALSE),
                              conditionalPanel("input.checkPirna == 1",
                                               sliderInput("selectPiMM", "Minimum Hamming to known piRNAs", min=0, max=4, value=3, step=1)),
                              checkboxInput("checkEnzySites", "Remove Restriction Sites", value = FALSE),
                              conditionalPanel("input.checkEnzySites == 1",
                                               # The list of enzymes (in groups) is filled from Enzymes.csv in app_server.R.
                                               # Nothing is pre-selected: the user chooses.
                                               selectizeInput("Oenzymes", "Select Enzymes to Remove",
                                                              choices = NULL, selected = NULL, multiple = TRUE,
                                                              options = list(placeholder = "Choose enzymes...")))),

              accordion_panel("4. Gene Architecture",
                              checkboxInput("checkboxRibo", "Optimize Ribosomal Binding", value = FALSE),
                              checkboxInput("checkfouras", "Add Consensus Start (aaaaATG)", value = FALSE),

                              checkboxInput("checkIntron", "Add Introns", value = FALSE),
                              conditionalPanel("input.checkIntron == 1",
                                               selectInput("intropt", "Intron Type", choices = c("Loading...")),
                                               conditionalPanel("input.intropt == 'Synthetic'",
                                                                sliderInput("num_synth_introns", "Number of Synthetic Introns", min=1, max=10, value=3, step=1),
                                                                helpText("Average C. elegans exon length is 200 bp (Speith et al., Wormbook, 2014)")),
                                               radioButtons("intdistop", "Intron Placement", choices = list("Early start" = 1, "Equidistant" = 2), inline = TRUE),
                                               checkboxInput("checkintframe", "Force in reading frame", value = FALSE)),

                              checkboxInput("checkUTRs", "Append UTRs", value = FALSE),
                              conditionalPanel("input.checkUTRs == 1",
                                               selectInput("p5UTR", "5' UTR", choices = c("Loading...")),
                                               selectInput("p3UTR", "3' UTR", choices = c("Loading..."))),

                              checkboxInput("checkPromoters", "Add Promoter Element", value = FALSE),
                              conditionalPanel("input.checkPromoters == 1",
                                               selectInput("oppromo", "Select Promoter", choices = c("Loading...")))
              ),

              accordion_panel("5. Synthesis Checks",
                              checkboxInput("checkAnal", "Analytical Mode Only (No changes)", value = FALSE),
                              checkboxInput("checkTwisty", "Check Twist Bioscience Guidelines", value = FALSE))
            ),

            br(),
            actionButton("actionSeq", "Optimize Sequence", class = "btn-success btn-lg w-100"),
            br(),
            actionButton("actionRESET", "Reset Form", class = "btn-outline-danger w-100")
          ),

          div(
            id = "results_area",
            uiOutput("AllResults")
          )
        )
      ),

      nav_panel(
        title = "About",
        card(
          card_header("Transgene Builder (Beta)"),
          card_body(
            p("This application optimizes C. elegans transgenes using modern codon frequency tables and biological constraints."),
            p("Developed by Amhed Velazquez & Christian-Froekjaer Jensen, Laboratory of Synthetic Genome Biology, KAUST."),
            p(style = "color: #7f8c8d; font-size: 0.85em;", "Software versions: ", textOutput("software_versions", inline = TRUE))
          )
        )
      )
    )
  )
}

#' Add external resources to the application
#'
#' @importFrom golem add_resource_path favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path("www", app_sys("app/www"))
  tags$head(
    favicon(),
    bundle_resources(path = app_sys("app/www"), app_title = "WormTransgeneBuilder"),
    tags$script(src ="www/sequence-viewer.bundle.js")
  )
}
