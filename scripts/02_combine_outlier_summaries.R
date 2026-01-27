#!/usr/bin/env Rscript
# ==============================================================================
# Combine Outlier Summaries from Multiple Samples
# ==============================================================================
#
# Usage: Rscript scripts/02_combine_outlier_summaries.R <samples_csv> <out_base> <obj_base> [bin_size]
#
# Output: Single CSV file with metadata header:
#   {output_folder}_{config}_{YYYYMMDD}.csv
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
out_base    <- if (length(args) >= 2) args[[2]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024"
obj_base    <- if (length(args) >= 3) args[[3]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute"
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(dplyr)

run_date <- format(Sys.time(), "%Y%m%d")

if (is.null(samples_csv) || !file.exists(samples_csv)) {
  stop("samples_csv is required and must exist.")
}

samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")

# Normalize column names
normalize_colname <- function(df, candidates, target) {
  for (c in candidates) if (c %in% names(df)) { names(df)[names(df) == c] <- target; break }
  df
}
samples_tbl <- normalize_colname(samples_tbl, c("Cell line", "Cell_line", "category"), "category")
samples_tbl <- normalize_colname(samples_tbl, c("HRD_HRP"), "feature1")
samples_tbl <- normalize_colname(samples_tbl, c("Resistance"), "feature2")

for (col in c("category", "feature1", "feature2")) {
  if (!col %in% names(samples_tbl)) samples_tbl[[col]] <- NA
}

message("Processing ", nrow(samples_tbl), " samples from ", samples_csv)
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Collect data and QC params
combined <- tibble()
qc_params_list <- list()
config_name <- "unknown"

for (i in seq_len(nrow(samples_tbl))) {
  sample_id  <- samples_tbl$sample[i]
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  sample_dir <- file.path(obj_base, sample_obj)

  # Read summary (RDS or CSV)
  csv_path <- file.path(sample_dir, "outlier_summary.csv")
  rds_path <- file.path(sample_dir, "outlier_summary.rds")
  qc_csv_path <- file.path(sample_dir, "qc_criteria_used.csv")

  summary <- NULL
  if (file.exists(csv_path)) {
    summary <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  } else if (file.exists(rds_path)) {
    summary <- tryCatch(as.data.frame(readRDS(rds_path)), error = function(e) NULL)
  }

  if (is.null(summary)) {
    warning("Missing summary for: ", sample_id)
    next
  }

  # Remove any QC columns from data (shouldn't be there now, but just in case)
  qc_cols <- grep("^QC_", names(summary), value = TRUE)
  if (length(qc_cols) > 0) {
    summary <- summary %>% select(-any_of(qc_cols))
  }

  # Add metadata columns
  summary$Sample      <- as.character(sample_id)
  summary$`Cell line` <- as.character(samples_tbl$category[i])
  summary$feature1    <- as.character(samples_tbl$feature1[i])
  summary$feature2    <- as.character(samples_tbl$feature2[i])

  # Read QC criteria from first sample
  if (length(qc_params_list) == 0 && file.exists(qc_csv_path)) {
    qc_params <- tryCatch(read.csv(qc_csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(qc_params)) {
      for (col in names(qc_params)) {
        qc_params_list[[col]] <- qc_params[[col]][1]
      }
      if ("config_file" %in% names(qc_params)) {
        config_name <- tools::file_path_sans_ext(qc_params$config_file[1])
      }
    }
  }

  combined <- bind_rows(combined, summary)
}

if (nrow(combined) == 0) stop("No outlier summaries found")

# Add totals row
numeric_cols <- names(combined)[sapply(combined, is.numeric)]
totals <- combined %>%
  summarise(across(all_of(numeric_cols), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(Sample = "Total", `Cell line` = "", feature1 = "", feature2 = "")
combined <- bind_rows(combined, totals)

# Reorder columns
metadata_cols <- c("Sample", "Cell line", "feature1", "feature2")
metric_cols <- setdiff(names(combined), metadata_cols)
final_order <- c(intersect(metadata_cols, names(combined)), metric_cols)
combined <- combined %>% select(all_of(final_order))

# Generate output filename
output_name <- basename(normalizePath(out_base, mustWork = FALSE))
outfile <- file.path(out_base, paste0(output_name, "_", config_name, "_", run_date, ".csv"))

# Write CSV with metadata header
con <- file(outfile, "w")
writeLines("# ==============================================================================", con)
writeLines("# QC Summary Report", con)
writeLines("# ==============================================================================", con)
writeLines(paste0("# Run Date: ", run_date), con)
writeLines(paste0("# Samples CSV: ", basename(samples_csv)), con)
writeLines(paste0("# Samples Processed: ", nrow(combined) - 1), con)
writeLines(paste0("# Bin Size: ", bin_size), con)
writeLines("#", con)
writeLines("# QC Parameters:", con)

if (length(qc_params_list) > 0) {
  for (param in names(qc_params_list)) {
    writeLines(paste0("#   ", param, ": ", qc_params_list[[param]]), con)
  }
}

writeLines("# ==============================================================================", con)
close(con)

write.table(combined, outfile, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)

message("Wrote: ", outfile)
message("\nSamples: ", nrow(combined) - 1)
message("Config: ", config_name)
