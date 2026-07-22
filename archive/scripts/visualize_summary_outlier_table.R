#!/usr/bin/env Rscript

# Adapted visualization script — reads a samples CSV (or uses fallback),
# loads per-sample `outlier_summary.rds` files produced by the pipeline,
# and writes combined tables (CSV, PNG, PDF).

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
out_base <- if (length(args) >= 2) args[[2]] else "/Volumes/Fl/sc_analysis/all_samples_23july2024"
obj_base <- if (length(args) >= 3) args[[3]] else "/Volumes/Fl/sc_analysis/post-scAbsolute"
bin_size <- if (length(args) >= 4) args[[4]] else "100"

library(dplyr)
library(knitr)
# Optional rendering packages are used only when available; script will fall back to CSV when missing.

# source local helpers (if needed for future extensions)
if (file.exists("R/core.R")) source("R/core.R")
if (file.exists("R/summary_helpers.R")) source("R/summary_helpers.R")

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Read sample manifest if provided, otherwise use a small fallback
if (!is.null(samples_csv) && file.exists(samples_csv)) {
  samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE)
  if (!"sample" %in% colnames(samples_tbl)) stop("samples CSV must contain a 'sample' column")
  if (!"category" %in% colnames(samples_tbl)) samples_tbl$category <- NA
} else {
  message("No samples CSV provided or not found — using a small fallback list. To run for many samples, pass a CSV as first arg.")
  samples_tbl <- data.frame(sample = c("23355", "25146", "25148"), category = NA, stringsAsFactors = FALSE)
}

# initialize
combined_outlier_summary <- tibble()

for (i in seq_len(nrow(samples_tbl))) {
  SAMPLENAME <- samples_tbl$sample[i]
  CATEGORY <- samples_tbl$category[i]
  SAMPLEONJ <- paste0("SLX-", SAMPLENAME, "_", bin_size)
  OUTPUT <- file.path(obj_base, SAMPLEONJ)
  outlier_rds <- file.path(OUTPUT, "outlier_summary.rds")

  if (!file.exists(outlier_rds)){
    warning("Missing outlier summary: ", outlier_rds, " — skipping sample ", SAMPLENAME)
    next
  }

  outlier_summary <- tryCatch(readRDS(outlier_rds), error = function(e){ warning("Failed to read: ", outlier_rds); return(NULL) })
  if (is.null(outlier_summary)) next

  # Ensure outlier_summary is a data.frame / tibble
  outlier_summary <- as.data.frame(outlier_summary)
  outlier_summary$Sample <- SAMPLENAME
  outlier_summary$Category <- ifelse(is.na(CATEGORY) || CATEGORY == "", NA, CATEGORY)

  combined_outlier_summary <- bind_rows(combined_outlier_summary, outlier_summary)
}

if (nrow(combined_outlier_summary) == 0){
  stop("No outlier summaries found for any samples — nothing to visualize.")
}

# Optionally order by category
sorted_data <- combined_outlier_summary %>% arrange(Category)

# Ensure Sample/Category are character to avoid type conflicts when adding totals row
if ("Sample" %in% colnames(sorted_data)) sorted_data$Sample <- as.character(sorted_data$Sample)
if ("Category" %in% colnames(sorted_data)) sorted_data$Category <- as.character(sorted_data$Category)

# Add totals row summarizing numeric columns
sum_row <- sorted_data %>% summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(Sample = "Total", Category = "")
sorted_data <- bind_rows(sorted_data, sum_row)

# write results
outfile_csv <- file.path(out_base, "combined_outlier_summary.csv")
write.csv(sorted_data, outfile_csv, row.names = FALSE)
message("Wrote: ", outfile_csv)

# render table to PDF (LaTeX) if kableExtra is available, otherwise skip
if (requireNamespace("kableExtra", quietly = TRUE)) {
  outfile_pdf <- file.path(out_base, "combined_outlier_summary.pdf")
  pdf(outfile_pdf, height = 6, width = 8)
  print(kable(sorted_data, format = "latex", booktabs = TRUE) %>% kableExtra::kable_styling(latex_options = c("striped", "hold_position")))
  dev.off()
  message("Wrote: ", outfile_pdf)
} else {
  message("Package 'kableExtra' not available — skipping PDF rendering. CSV was written to: ", outfile_csv)
}

# try to render PNG via flextable if available
if (requireNamespace("flextable", quietly = TRUE) && (requireNamespace("webshot", quietly = TRUE) || requireNamespace("webshot2", quietly = TRUE))) {
  outfile_png_html <- file.path(out_base, "combined_outlier_summary.png")
  ft <- flextable::flextable(sorted_data)
  tryCatch({
    flextable::save_as_image(ft, path = outfile_png_html)
    message("Wrote: ", outfile_png_html)
  }, error = function(e){
    message("Could not write PNG via flextable/webshot: ", conditionMessage(e))
  })
} else {
  message("Packages 'flextable' and/or 'webshot' not available — skipping PNG rendering. CSV was written to: ", outfile_csv)
}

invisible(sorted_data)
