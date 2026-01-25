#!/usr/bin/env Rscript
# ==============================================================================
# Combine Outlier Summaries from Multiple Samples
# ==============================================================================
#
# Usage: Rscript scripts/02_combine_outlier_summaries.R <samples_csv> <out_base> <obj_base> [bin_size]
#
# This script:
#   - Reads individual outlier_summary.csv files (which include QC criteria)
#   - Combines them into a single table with metadata
#   - Creates timestamped archive to prevent overwriting
#   - Saves run metadata for reproducibility
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
out_base    <- if (length(args) >= 2) args[[2]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024"
obj_base    <- if (length(args) >= 3) args[[3]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute"
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(dplyr)
library(knitr)

# Generate run identifier based on date only (not time)
# This means same-day runs with same config will overwrite each other
run_date <- format(Sys.time(), "%Y%m%d")

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

# We'll determine the config name from the first sample's data later
config_name <- "unknown"

# Combine outlier summaries - prefer CSV (has QC columns) over RDS
combined <- tibble()
qc_params_combined <- tibble()

for (i in seq_len(nrow(samples_tbl))) {
  sample_id  <- samples_tbl$sample[i]
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  sample_dir <- file.path(obj_base, sample_obj)

  # Try CSV first (has QC criteria), fall back to RDS
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

  # Add metadata
  summary$Sample      <- as.character(sample_id)
  summary$`Cell line` <- as.character(samples_tbl$category[i])
  summary$feature1    <- as.character(samples_tbl$feature1[i])
  summary$feature2    <- as.character(samples_tbl$feature2[i])

  # Read QC criteria if available
  if (file.exists(qc_csv_path)) {
    qc_params <- tryCatch(read.csv(qc_csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(qc_params)) {
      qc_params$Sample <- as.character(sample_id)
      qc_params_combined <- bind_rows(qc_params_combined, qc_params)
      # Extract config name from first sample
      if (config_name == "unknown" && "config_file" %in% names(qc_params)) {
        config_name <- tools::file_path_sans_ext(qc_params$config_file[1])
      }
    }
  }

  # Also try to get config name from summary CSV (QC_config_file column)
  if (config_name == "unknown" && "QC_config_file" %in% names(summary)) {
    cfg <- summary$QC_config_file[1]
    if (!is.na(cfg) && nzchar(cfg)) {
      config_name <- tools::file_path_sans_ext(cfg)
    }
  }

  combined <- bind_rows(combined, summary)
}

if (nrow(combined) == 0) stop("No outlier summaries found")

# Create archive directory based on: output folder + config + date
# Same output + same config + same date = same folder (overwrites)
output_name <- basename(normalizePath(out_base, mustWork = FALSE))
archive_dir <- file.path(out_base, "archive", paste0(output_name, "_", config_name, "_", run_date))
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

# Identify QC columns (start with QC_)
qc_cols <- grep("^QC_", names(combined), value = TRUE)

# Add totals row (only for numeric non-QC columns)
numeric_cols <- setdiff(names(combined)[sapply(combined, is.numeric)], qc_cols)
totals <- combined %>%
  summarise(across(all_of(numeric_cols), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(Sample = "Total", `Cell line` = "", feature1 = "", feature2 = "")

# For QC columns, just take the first value (they should all be the same if using same config)
for (col in qc_cols) {
  if (col %in% names(combined)) {
    totals[[col]] <- combined[[col]][1]
  }
}

combined <- bind_rows(combined, totals)

# Reorder columns: metadata first, then metrics, then QC params
metadata_cols <- c("Sample", "Cell line", "feature1", "feature2")
metric_cols <- setdiff(setdiff(names(combined), metadata_cols), qc_cols)
final_order <- c(intersect(metadata_cols, names(combined)),
                 intersect(metric_cols, names(combined)),
                 intersect(qc_cols, names(combined)))
combined <- combined %>% select(all_of(final_order))

# Write outputs
# 1. Main output (will be overwritten each run)
outfile <- file.path(out_base, "outlier_summary_table_meta_combined.csv")
write.csv(combined, outfile, row.names = FALSE)
message("Wrote: ", outfile)

# 2. Timestamped archive (never overwritten)
archive_file <- file.path(archive_dir, "outlier_summary_table_meta_combined.csv")
write.csv(combined, archive_file, row.names = FALSE)
message("Archived: ", archive_file)

# 3. Save run metadata
run_info <- data.frame(
  run_timestamp = run_timestamp,
  samples_csv = normalizePath(samples_csv, mustWork = FALSE),
  obj_base = obj_base,
  out_base = out_base,
  bin_size = bin_size,
  n_samples = nrow(samples_tbl),
  n_processed = nrow(combined) - 1,  # exclude totals row
  stringsAsFactors = FALSE
)
run_info_file <- file.path(archive_dir, "run_info.csv")
write.csv(run_info, run_info_file, row.names = FALSE)
message("Run info: ", run_info_file)

# 4. Save combined QC parameters (if we collected any)
if (nrow(qc_params_combined) > 0) {
  qc_file <- file.path(archive_dir, "qc_params_all_samples.csv")
  write.csv(qc_params_combined, qc_file, row.names = FALSE)
  message("QC params: ", qc_file)
}

# Optional PNG rendering
if (requireNamespace("kableExtra", quietly = TRUE) && requireNamespace("webshot2", quietly = TRUE)) {
  outfile_png <- file.path(out_base, "outlier_summary_table_meta.png")
  kable(combined, format = "html") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    kableExtra::row_spec(nrow(combined), bold = TRUE, background = "#EFEFEF") %>%
    kableExtra::save_kable(file = outfile_png, zoom = 2)
  message("Wrote: ", outfile_png)
}

message("\n=== Run Complete ===")
message("Date: ", run_date)
message("Output: ", output_name)
message("Config: ", config_name)
message("Samples processed: ", nrow(combined) - 1)
message("Archive: ", basename(archive_dir))
if (length(qc_cols) > 0) {
  message("QC columns: ", paste(qc_cols, collapse = ", "))
}
message("====================\n")
