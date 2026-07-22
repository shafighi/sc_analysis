#!/usr/bin/env Rscript
# ==============================================================================
# Generate Read Distribution Summaries and Plots
# ==============================================================================
#
# Usage:
#   Rscript scripts/07_generate_read_distribution.R <samples_csv> <obj_base> <out_base> [bin_size]
#
# Arguments:
#   samples_csv - CSV file with 'sample' column (and optional 'category')
#   obj_base    - Directory containing scAbsolute RDS objects
#   out_base    - Output directory for results
#   bin_size    - Bin size (default: 100)
#
# Output per sample:
#   - figures/read_distribution_total_reads.pdf : Histogram of total.reads
#   - figures/read_distribution_rpc.pdf         : Histogram of rpc
#   - figures/read_distribution_scatter.pdf     : Scatter plot of rpc vs total.reads
#
# Cross-sample output:
#   - read_distribution_summary.csv : Summary statistics per sample
#   - combined_reads.rds            : Combined read counts across all samples
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript scripts/07_generate_read_distribution.R <samples_csv> <obj_base> <out_base> [bin_size]")
}
samples_csv <- args[[1]]
obj_base    <- args[[2]]
out_base    <- args[[3]]
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(Biobase)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Load samples
# ==============================================================================
if (file.exists(samples_csv)) {
  samples_tbl <- read.csv(samples_csv)
  if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")
} else {
  stop("Samples CSV file is required")
}

# ==============================================================================
# Initialize data structures
# ==============================================================================
sample_summary <- data.frame(
  SLXNumber = character(),
  SampleName = character(),
  MeanReads = numeric(),
  MedianReads = numeric(),
  MeanRPC = numeric(),
  MedianRPC = numeric(),
  NumCells = numeric(),
  stringsAsFactors = FALSE
)

combined_reads <- c()
combined_rpc <- c()

# ==============================================================================
# Process each sample
# ==============================================================================
for (i in seq_len(nrow(samples_tbl))) {
  sample_id <- samples_tbl$sample[i]
  category  <- if ("category" %in% names(samples_tbl)) samples_tbl$category[i] else NA
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  input_file <- file.path(obj_base, paste0(sample_obj, ".rds"))
  output_dir <- file.path(out_base, sample_obj)

  if (!file.exists(input_file)) {
    warning("File not found: ", input_file, " - skipping")
    next
  }

  message("Processing: ", sample_id)

  object <- tryCatch(readRDS(input_file), error = function(e) { warning(e$message); NULL })
  if (is.null(object)) next

  object_df <- Biobase::pData(object)

  # Calculate statistics
  mean_reads <- mean(object_df$total.reads, na.rm = TRUE)
  median_reads <- median(object_df$total.reads, na.rm = TRUE)
  mean_rpc <- mean(object_df$rpc, na.rm = TRUE)
  median_rpc <- median(object_df$rpc, na.rm = TRUE)
  num_cells <- nrow(object_df)

  # Add to summary
  sample_summary <- rbind(
    sample_summary,
    data.frame(
      SLXNumber = sample_id,
      SampleName = if (!is.na(category)) category else sample_id,
      MeanReads = mean_reads,
      MedianReads = median_reads,
      MeanRPC = mean_rpc,
      MedianRPC = median_rpc,
      NumCells = num_cells
    )
  )

  # Add to combined
  combined_reads <- c(combined_reads, object_df$total.reads)
  combined_rpc <- c(combined_rpc, object_df$rpc)

  # Create figures directory
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  # Generate plots
  # Histogram of total.reads
  p1 <- ggplot(object_df, aes(x = total.reads)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
    labs(title = paste("Total Reads Distribution -", sample_id),
         x = "Total Reads", y = "Frequency") +
    theme_minimal()
  ggsave(file.path(figures_dir, "read_distribution_total_reads.pdf"), p1, width = 6, height = 4)

  # Histogram of rpc
  p2 <- ggplot(object_df, aes(x = rpc)) +
    geom_histogram(bins = 50, fill = "green", alpha = 0.7) +
    labs(title = paste("RPC Distribution -", sample_id),
         x = "RPC (Reads Per Cell)", y = "Frequency") +
    theme_minimal()
  ggsave(file.path(figures_dir, "read_distribution_rpc.pdf"), p2, width = 6, height = 4)

  # Scatter plot of rpc vs total.reads
  p3 <- ggplot(object_df, aes(x = total.reads, y = rpc)) +
    geom_point(alpha = 0.6) +
    labs(title = paste("RPC vs Total Reads -", sample_id),
         x = "Total Reads", y = "RPC (Reads Per Cell)") +
    theme_minimal()
  ggsave(file.path(figures_dir, "read_distribution_scatter.pdf"), p3, width = 6, height = 4)
}

# ==============================================================================
# Save cross-sample outputs
# ==============================================================================
output_dir <- file.path(out_base, "read_distribution")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(sample_summary, file.path(output_dir, "read_distribution_summary.csv"), row.names = FALSE)
saveRDS(list(reads = combined_reads, rpc = combined_rpc), file.path(output_dir, "combined_reads.rds"))

message("Read distribution summaries and plots saved.")
