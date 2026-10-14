# Compares the weights returned by the Codon Adapter (codon_adapter_weights.csv, made by query_codon_adapter.R) with the app's codon table.
cf <- read.csv("inst/extdata/Codon_frequencies.csv", stringsAsFactors = FALSE); cf$Codon <- toupper(cf$Codon)
tool <- read.csv("data-raw/codon_adapter/codon_adapter_weights.csv", stringsAsFactors = FALSE)
m <- merge(tool, cf, by.x = "codon", by.y = "Codon")
cat("codons in both:", nrow(m), "\n")
sense <- m$Amino != "X"
d <- abs(m$weight - m$BringMansWeights)
cat(sprintf("weights (sense codons): max |tool - app| = %.4f; codons that differ by more than 0.0005: %d of %d\n", max(d[sense]), sum(d[sense] > 5e-4), sum(sense)))
print(m[sense & d > 5e-4, c("codon", "Amino", "weight", "BringMansWeights")], row.names = FALSE)
cat("stop codons (tool weight / app BringMansWeights):", paste(m$codon[!sense], m$weight[!sense], m$BringMansWeights[!sense], collapse = "; "), "\n")
# best codon per amino acid: the tool's weight-1 codon, and the codon its optimizer chose
best_tool <- tapply(m$codon[sense & m$weight == 1], m$Amino[sense & m$weight == 1], paste, collapse = "/")
best_app <- tapply(m$codon[sense & m$CAI1 == 1], m$Amino[sense & m$CAI1 == 1], paste, collapse = "/")
cat("amino acids with a different best codon (weight 1 in the tool vs CAI1 in the app):", sum(best_tool != best_app[names(best_tool)]), "\n")
chosen <- tapply(m$optimized[sense], m$Amino[sense], function(x) paste(sort(unique(x)), collapse = "/"))
cat("the tool's optimizer always chose the weight-1 codon:", all(chosen == best_tool[names(chosen)]), "\n")
