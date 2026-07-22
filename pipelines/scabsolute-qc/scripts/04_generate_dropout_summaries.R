#!/usr/bin/env Rscript
# ==============================================================================
# Generate Dropout Summaries for scAbsolute Objects
# ==============================================================================
#
# Analyzes copy number dropouts (bins with segVal == 0) per cell and chromosome.
# Requires step 01 output (qc_cell_labels.rds) for QC status assignment.
#
# Usage:
#   Rscript scripts/04_generate_dropout_summaries.R <samples_csv> <obj_base> <out_base> [bin_size]
#
# Arguments:
#   samples_csv - CSV with 'sample' column (and optional metadata)
#   obj_base    - Directory containing scAbsolute RDS objects
#   out_base    - Directory containing step 01 output (also used for dropout output)
#   bin_size    - Bin size (default: 100)
#
# Per-sample output (in out_base/SLX-{sample}_{bin_size}/):
#   - cn_binned.rds              : Long-format copy number data (bin-level, all cells)
#   - dropout_per_cell_chr.rds   : Dropout counts per cell x chromosome (with QC status)
#   - dropout_per_cell.csv       : Per-cell total dropout count + QC status
#
# Combined output (in out_base/):
#   - dropout_summary_combined.csv : Dropout freq per cell by sample x chromosome
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript scripts/04_generate_dropout_summaries.R <samples_csv> <obj_base> <out_base> [bin_size]")
}
samples_csv <- args[[1]]
obj_base    <- args[[2]]
out_base    <- args[[3]]
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(Biobase)
library(dplyr)
library(tidyr)
source("R/dropout_helpers.R")

# ==============================================================================
# Load samples
# ==============================================================================
if (!file.exists(samples_csv)) {
  stop("Samples CSV does not exist: ", samples_csv)
}

samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")

# Normalize metadata column names
normalize_colname <- function(df, candidates, target) {
  for (c in candidates) if (c %in% names(df)) { names(df)[names(df) == c] <- target; break }
  df
}
samples_tbl <- normalize_colname(samples_tbl, c("Cell line", "Cell_line", "category"), "category")
for (col in c("category")) {
  if (!col %in% names(samples_tbl)) samples_tbl[[col]] <- NA
}

message("Processing ", nrow(samples_tbl), " samples from ", samples_csv)

# ==============================================================================
# Process each sample
# ==============================================================================
combined_dropout <- data.frame()

for (i in seq_len(nrow(samples_tbl))) {
  sample_id  <- samples_tbl$sample[i]
  category   <- if ("category" %in% names(samples_tbl)) samples_tbl$category[i] else NA
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  input_file <- file.path(obj_base, paste0(sample_obj, ".rds"))
  output_dir <- file.path(out_base, sample_obj)

  if (!file.exists(input_file)) {
    warning("Object not found: ", input_file, " - skipping")
    next
  }

  # qc_cell_labels from step 01
  outlier_path <- file.path(output_dir, "qc_cell_labels.rds")
  if (!file.exists(outlier_path)) {
    warning("qc_cell_labels.rds not found for ", sample_id,
            " - run step 01 first. Skipping.")
    next
  }

  message("Processing: ", sample_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  object <- tryCatch(readRDS(input_file), error = function(e) {
    warning(e$message); NULL
  })
  if (is.null(object)) next

  cellbased_outliers <- readRDS(outlier_path)

  # Compute dropouts
  res <- tryCatch(
    get_dropout_by_chromosome(object),
    error = function(e) { warning("Dropout analysis failed for ", sample_id, ": ", e$message); NULL }
  )
  if (is.null(res)) next

  # Assign QC status to each cell
  all_cells <- unique(res$cell_chr_dropout$cell)
  cell_status <- data.frame(
    cell   = all_cells,
    status = assign_cell_status(all_cells, cellbased_outliers),
    stringsAsFactors = FALSE
  )
  res$cell_chr_dropout <- res$cell_chr_dropout %>%
    left_join(cell_status, by = "cell")

  # Save per-sample outputs
  saveRDS(res$cn_binned, file.path(output_dir, "cn_binned.rds"))
  saveRDS(res$cell_chr_dropout, file.path(output_dir, "dropout_per_cell_chr.rds"))

  # Per-cell summary: total dropouts + status
  cell_dropout_status <- res$cell_chr_dropout %>%
    group_by(cell, status) %>%
    summarise(total_dropout = sum(dropout_count, na.rm = TRUE), .groups = "drop")
  write.csv(cell_dropout_status, file.path(output_dir, "dropout_per_cell.csv"),
            row.names = FALSE)

  # Aggregate for combined summary: dropout freq per chromosome (normalized by cell count)
  n_cells <- length(unique(res$cn_binned$cell))
  chr_summary <- res$cn_binned %>%
    filter(segVal == 0) %>%
    count(chromosome, name = "total_dropouts") %>%
    mutate(
      freq_per_cell = total_dropouts / n_cells,
      sample        = sample_id,
      category      = category,
      n_cells       = n_cells
    )

  combined_dropout <- bind_rows(combined_dropout, chr_summary)
  message("Done: ", sample_id, " (", n_cells, " cells)")
}

# ==============================================================================
# Write combined summary
# ==============================================================================
if (nrow(combined_dropout) > 0) {
  # Ensure all sample x chromosome combinations exist
  all_chroms <- as.character(1:22)
  all_samples <- unique(combined_dropout$sample)
  complete_grid <- expand.grid(
    sample     = all_samples,
    chromosome = all_chroms,
    stringsAsFactors = FALSE
  )

  # Get sample metadata for filling
  sample_meta <- combined_dropout %>%
    select(sample, category, n_cells) %>%
    distinct()

  combined_dropout <- complete_grid %>%
    left_join(combined_dropout, by = c("sample", "chromosome")) %>%
    left_join(sample_meta, by = "sample", suffix = c("", ".meta")) %>%
    mutate(
      total_dropouts = replace_na(total_dropouts, 0),
      freq_per_cell  = replace_na(freq_per_cell, 0),
      category       = coalesce(category, category.meta),
      n_cells        = coalesce(n_cells, n_cells.meta)
    ) %>%
    select(sample, category, chromosome, total_dropouts, freq_per_cell, n_cells) %>%
    arrange(sample, as.integer(chromosome))

  outfile <- file.path(out_base, "dropout_summary_combined.csv")
  write.csv(combined_dropout, outfile, row.names = FALSE)
  message("\nWrote combined summary: ", outfile)
} else {
  message("\nNo dropout data generated.")
}

message("\nAll samples processed.")
