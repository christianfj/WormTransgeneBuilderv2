# NOTE: this script builds data/glo_full.rda, a named vector of GLO scores, from the raw text file
# data-raw/sequence_lib_scores.txt (not in the repository). That data file has since been removed:
# the app now reads inst/extdata/glo_scores.rds, the same scores without the 12-mer names (see
# dev/make_glo_scores.R and R/fct_glo_scores.R). Kept as a record of where the scores came from.

cat("1. Reading text file (this uses the most RAM)...\n")
lines <- readLines("data-raw/sequence_lib_scores.txt")
start_idx <- which(lines == "HEADER=END")[1] + 1
data_lines <- lines[start_idx:length(lines)]
# Clean up any trailing empty lines or 'DATA=END'
data_lines <- data_lines[data_lines != "" & data_lines != "DATA=END"]
rm(lines); gc()

cat("2. Extracting Keys and Values...\n")
# Extract alternating lines
key_idx <- seq(1, length(data_lines), by = 2)
val_idx <- seq(2, length(data_lines), by = 2)

# Ensure lengths match (The Off-by-One Fix)
common_len <- min(length(key_idx), length(val_idx))
key_idx <- key_idx[1:common_len]
val_idx <- val_idx[1:common_len]

raw_keys <- trimws(data_lines[key_idx])
raw_values <- as.numeric(data_lines[val_idx])
rm(data_lines); gc()

cat("3. Creating and Saving Full Dictionary (xz compression)...\n")
glo_full <- raw_values
names(glo_full) <- raw_keys
# This is the line that will take a few minutes!
save(glo_full, file = "data/glo_full.rda", compress = "xz")
rm(glo_full); gc()

cat("4. Creating and Saving Fast Dictionary (xz compression)...\n")
q_low <- quantile(raw_values, 0.05, na.rm = TRUE)
q_high <- quantile(raw_values, 0.95, na.rm = TRUE)
glo_median_score <- median(raw_values, na.rm = TRUE)

fast_idx <- which(raw_values <= q_low | raw_values >= q_high)
glo_fast <- raw_values[fast_idx]
names(glo_fast) <- raw_keys[fast_idx]

save(glo_fast, file = "data/glo_fast.rda", compress = "xz")
save(glo_median_score, file = "data/glo_median_score.rda")

cat("\n--- SUCCESS! 16,777,216 entries processed. ---\n")
