#' The piRNA Databases of the Supported Species
#'
#' The app can screen a transgene against the piRNAs of one or several *Caenorhabditis* species. Each species has its
#' own database file in `inst/extdata/`, and the databases are independent: no entry is shared and nothing is combined
#' unless the user selects more than one species. See `docs/PIRNA_DATABASE.md` and `docs/PIRNA_OTHER_SPECIES.md`.
#'
#' * `elegans`: `piRNA_Database.csv`, type 1 (WormBase WS298), type 2 and PRG-1-IP entries; the type 2 entries can be left out;
#' * `briggsae`: PRG-1 IP loci (Pastore et al. 2022), high confidence;
#' * `remanei`, `brenneri`: called from small-RNA sequencing (Shi et al. 2013), low confidence.
#' @noRd
NULL

# id, label shown in the app, file, confidence
pirna_species_table <- function() {
  data.frame(
    id = c("elegans", "briggsae", "remanei", "brenneri"),
    label = c("C. elegans", "C. briggsae", "C. remanei", "C. brenneri"),
    file = c("piRNA_Database.csv", "piRNA_Database_briggsae.csv.gz", "piRNA_Database_remanei.csv.gz", "piRNA_Database_brenneri.csv.gz"),
    confidence = c("main", "high", "low", "low"),
    stringsAsFactors = FALSE)
}

# Reads one species' database (kept after the first read; a shared R process serves every session).
pirna_species_cache <- new.env(parent = emptyenv())
read_pirna_species <- function(id) {
  if (!is.null(pirna_species_cache[[id]])) return(pirna_species_cache[[id]])
  info <- pirna_species_table()
  info <- info[info$id == id, ]
  if (nrow(info) != 1) stop("Unknown species for the piRNA database: ", id, call. = FALSE)
  path <- system.file("extdata", info$file, package = "WormTransgeneBuilder")
  if (!nzchar(path)) path <- file.path("inst", "extdata", info$file)
  d <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
  d$Species <- info$label
  if (id != "elegans") {
    d$Type <- if (info$confidence == "high") "PRG-1 IP (high confidence)" else "small-RNA-seq call (low confidence)"
  }
  d <- d[, c("Species", "Type", "Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment")]
  pirna_species_cache[[id]] <- d
  d
}

#' Select the piRNAs to Screen Against
#'
#' @param species Character vector of species ids (`"elegans"`, `"briggsae"`, `"remanei"`, `"brenneri"`).
#' @param include_type2 Logical; whether the *C. elegans* type 2 piRNAs are included (ignored for other species).
#' @return A data frame with the columns `Species`, `Type`, `Name`, `Target_20mer`, `Full_21mer`, `Abundance`, `Enrichment`,
#'   one row per piRNA, *C. elegans* first, then the other species in the order briggsae, remanei, brenneri. The row order is
#'   the order of the distances that `get_pirna_distances()` returns for `Target_20mer`.
#' @noRd
select_pirna_database <- function(species = "elegans", include_type2 = TRUE) {
  info <- pirna_species_table()
  species <- info$id[info$id %in% species]
  if (length(species) == 0) stop("Select at least one species to screen piRNAs against.", call. = FALSE)
  parts <- lapply(species, function(id) {
    d <- read_pirna_species(id)
    if (id == "elegans" && !isTRUE(include_type2)) d <- d[d$Type != "type 2", ]
    d
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
