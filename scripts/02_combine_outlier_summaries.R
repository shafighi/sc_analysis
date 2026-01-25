#!/usr/bin/env Rscript
# Usage: Rscript scripts/02_combine_outlier_summaries.R <samples_csv> <out_base> <obj_base> [bin_size]

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
out_base    <- if (length(args) >= 2) args[[2]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024"
obj_base    <- if (length(args) >= 3) args[[3]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute"
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(dplyr)
library(knitr)

# Require samples CSV
if (is.null(samples_csv) || !file.exists(samples_csv)) {
  stop("samples_csv is required and must exist. Usage: Rscript 02_combine_outlier_summaries.R <samples_csv> <out_base> <obj_base> [bin_size]")
}

samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")

# Normalize column names to: sample, category, feature1, feature2
normalize_colname <- function(df, candidates, target) {
  for (c in candidates) if (c %in% names(df)) { names(df)[names(df) == c] <- target; break }
  df
}
samples_tbl <- normalize_colname(samples_tbl, c("Cell line", "Cell_line", "category"), "category")
samples_tbl <- normalize_colname(samples_tbl, c("HRD_HRP"), "feature1")
samples_tbl <- normalize_colname(samples_tbl, c("Resistance"), "feature2")

# Ensure required columns exist
for (col in c("category", "feature1", "feature2")) {
  if (!col %in% names(samples_tbl)) samples_tbl[[col]] <- NA
}

message("Processing ", nrow(samples_tbl), " samples from ", samples_csv)
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Combine outlier summaries
combined <- tibble()
for (i in seq_len(nrow(samples_tbl))) {
  sample_id  <- samples_tbl$sample[i]
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  rds_path   <- file.path(obj_base, sample_obj, "outlier_summary.rds")

  if (!file.exists(rds_path)) {
    warning("Missing: ", rds_path)
    next
  }

  summary <- tryCatch(readRDS(rds_path), error = function(e) { warning(e$message); NULL })
  if (is.null(summary)) next

  summary <- as.data.frame(summary)
  summary$Sample      <- as.character(sample_id)
  summary$`Cell line` <- as.character(samples_tbl$category[i])
  summary$feature1    <- as.character(samples_tbl$feature1[i])
  summary$feature2    <- as.character(samples_tbl$feature2[i])
  combined <- bind_rows(combined, summary)
}

if (nrow(combined) == 0) stop("No outlier summaries found")

# Add totals row
totals <- combined %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(Sample = "Total", `Cell line` = "", feature1 = "", feature2 = "")
combined <- bind_rows(combined, totals)

# Write output
outfile <- file.path(out_base, "outlier_summary_table_meta_combined.csv")
write.csv(combined, outfile, row.names = FALSE)
message("Wrote: ", outfile)

# Optional PNG rendering
if (requireNamespace("kableExtra", quietly = TRUE) && requireNamespace("webshot2", quietly = TRUE)) {
  outfile_png <- file.path(out_base, "outlier_summary_table_meta.png")
  kable(combined, format = "html") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    kableExtra::row_spec(nrow(combined), bold = TRUE, background = "#EFEFEF") %>%
    kableExtra::save_kable(file = outfile_png, zoom = 2)
  message("Wrote: ", outfile_png)
}
