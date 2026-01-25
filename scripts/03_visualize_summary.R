#!/usr/bin/env Rscript
# ==============================================================================
# Visualize QC Summary and Generate Final Output
# ==============================================================================
#
# Usage: Rscript scripts/03_visualize_summary.R <input_csv> <out_base> [group_col] [label_col]
#
# This script:
#   - Reads combined outlier summary CSV (with QC parameters)
#   - Generates visualization plots
#   - Creates QC_summary.csv with all columns including QC parameters
#   - Archives outputs with timestamps
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else NULL
out_base  <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "."
group_col <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else NULL
label_col <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else NULL

if (is.null(input_csv) || !file.exists(input_csv)) {
  stop("Usage: Rscript 03_visualize_summary.R <input_csv> <out_base> [group_col] [label_col]")
}

# Generate run date (not timestamp) - same day runs with same config overwrite
run_date <- format(Sys.time(), "%Y%m%d")

library(ggplot2)
library(tidyr)
library(dplyr)

# Helper: find first matching column
pick_col <- function(df, candidates) {
  for (c in candidates) if (c %in% names(df)) return(c)
  NULL
}

# Load and prepare data
message("Reading: ", input_csv)
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
df <- df %>% filter(is.na(Sample) | Sample != "Total")

# Normalize key columns (internal names for processing)
col_map <- list(
  Sample      = c("Sample", "sample", "SampleID", "ID"),
  `post-scAbsolute` = c("Processed Cells", "Processed.Cells", "post-scAbsolute", "post-scAbsolute.Cells", "Total"),
  Good_Quality = c("Good Quality Cells", "Good.Quality.Cells", "Good.Quality", "GoodQuality"),
  Replicating = c("Replicating", "Replicating...RPC", "Replicating..RPC", "RPC"),
  `Outliers(RPC)` = c("RPC Outliers", "RPC.Outliers"),
  Normal      = c("Normal Cells", "Normal.Cells", "Normal", "Normal_Cells"),
  `Outliers(Alpha/Mapd/Gini)` = c("Alpha/Mapd/Gini", "Alpha.Mapd.Gini", "Outliers")
)
for (target in names(col_map)) {
  src <- pick_col(df, col_map[[target]])
  if (!is.null(src)) df[[target]] <- suppressWarnings(as.numeric(df[[src]]))
  else df[[target]] <- NA_real_
}

# Setup label and group columns (for plotting only, not in output)
if (is.null(label_col)) label_col <- pick_col(df, c("Sample", "Cell line", "category")) %||% "Sample"
label_values <- if (label_col %in% names(df)) as.character(df[[label_col]]) else seq_len(nrow(df))

if (is.null(group_col)) group_col <- pick_col(df, c("feature1", "feature2", "category", "Cell line"))
group_values <- if (!is.null(group_col) && group_col %in% names(df)) as.character(df[[group_col]]) else "All"

# Compute derived metrics with descriptive names
df <- df %>% mutate(
  Filtered = coalesce(Good_Quality, 0) - coalesce(Normal, 0),
  `(PassedQC+Replicating)/post-scAbsolute%` = round((coalesce(Good_Quality, 0) + coalesce(Replicating, 0)) * 100 / pmax(`post-scAbsolute`, 1)),
  `PassedQC/post-scAbsolute%` = round(coalesce(Good_Quality, 0) * 100 / pmax(`post-scAbsolute`, 1)),
  `(PassedQC-Normal)/post-scAbsolute%` = round((coalesce(Good_Quality, 0) - coalesce(Normal, 0)) * 100 / pmax(`post-scAbsolute`, 1))
)

# Rename columns to be more descriptive
df <- df %>% rename(
  `Outliers(Alpha/Mapd/Gini,post-RPC)` = `Outliers(Alpha/Mapd/Gini)`,
  `PassedQC(incl.Normal)` = Good_Quality,
  `PassedQC-Normal` = Filtered
)

# Identify QC columns from input (start with QC_)
qc_cols <- grep("^QC_", names(df), value = TRUE)

# Select and reorder columns for output: metadata first, then metrics, then QC params
metadata_cols <- c("Sample", "Cell line", "feature1", "feature2")
metric_cols <- c("post-scAbsolute", "Replicating", "Outliers(RPC)", "Outliers(Alpha/Mapd/Gini,post-RPC)",
                 "PassedQC(incl.Normal)", "Normal", "PassedQC-Normal",
                 "(PassedQC+Replicating)/post-scAbsolute%", "PassedQC/post-scAbsolute%", "(PassedQC-Normal)/post-scAbsolute%")
output_cols <- c(intersect(metadata_cols, names(df)),
                 intersect(metric_cols, names(df)),
                 intersect(qc_cols, names(df)))
df_output <- df %>% select(all_of(output_cols))

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
base_name <- tools::file_path_sans_ext(basename(input_csv))

# Determine config name from QC columns (if present)
config_name <- "unknown"
if (length(qc_cols) > 0 && "QC_config_file" %in% names(df)) {
  cfg <- unique(df$QC_config_file)
  cfg <- cfg[!is.na(cfg) & nzchar(cfg)]
  if (length(cfg) > 0) config_name <- tools::file_path_sans_ext(cfg[1])
}

# Use input CSV name as identifier for the sample set
input_name <- tools::file_path_sans_ext(basename(input_csv))

# Create archive directory based on: input + config + date
# Same input + same config + same date = same folder (overwrites)
archive_dir <- file.path(out_base, "archive", paste0(input_name, "_", config_name, "_", run_date))
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

# Output 1: Main QC_summary.csv (includes QC parameters)
qc_summary_file <- file.path(out_base, "QC_summary.csv")
write.csv(df_output, qc_summary_file, row.names = FALSE)
message("Wrote: ", qc_summary_file)

# Output 2: Timestamped archive
archive_qc_file <- file.path(archive_dir, "QC_summary.csv")
write.csv(df_output, archive_qc_file, row.names = FALSE)
message("Archived: ", archive_qc_file)

# Output 3: Also keep normalized version for backward compatibility
write.csv(df_output, file.path(out_base, paste0(base_name, "_normalized.csv")), row.names = FALSE)
message("Wrote: ", base_name, "_normalized.csv")

# Log QC parameters if present
if (length(qc_cols) > 0) {
  message("\nQC parameters in output:")
  for (col in qc_cols) {
    val <- unique(df_output[[col]])
    val <- val[!is.na(val)]
    if (length(val) > 0) message("  ", col, ": ", val[1])
  }
}

# Add plotting columns to df (not saved to CSV)
df$Label <- label_values
df$Group <- group_values
df$High_Quality_Pct <- df$`(PassedQC+Replicating)/post-scAbsolute%`

# Plot 1: Quality distribution
if (any(!is.na(df$High_Quality_Pct))) {
  p <- ggplot(df, aes(x = High_Quality_Pct)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    labs(title = "Distribution of (PassedQC+Replicating)/post-scAbsolute %", x = "(PassedQC+Replicating)/post-scAbsolute %", y = "Density") +
    theme_minimal()
  ggsave(file.path(out_base, paste0(base_name, "_quality_distribution.pdf")), p, width = 6, height = 4)
  message("Wrote: ", base_name, "_quality_distribution.pdf")
}

# Plot 2-3: Group comparisons (if multiple groups)
if (length(unique(df$Group)) > 1 && any(!is.na(df$High_Quality_Pct))) {
  # Density by group
  p <- ggplot(df, aes(x = High_Quality_Pct, fill = Group)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("(PassedQC+Replicating)/post-scAbsolute % by", group_col), x = "(PassedQC+Replicating)/post-scAbsolute %", y = "Density") +
    theme_minimal()
  ggsave(file.path(out_base, paste0(base_name, "_quality_by_group.pdf")), p, width = 7, height = 4)
  message("Wrote: ", base_name, "_quality_by_group.pdf")

  # Boxplot by group
  p <- ggplot(df, aes(x = Group, y = High_Quality_Pct, fill = Group)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    labs(title = paste("(PassedQC+Replicating)/post-scAbsolute % by", group_col), x = group_col, y = "(PassedQC+Replicating)/post-scAbsolute %") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(out_base, paste0(base_name, "_quality_boxplot.pdf")), p, width = 7, height = 5)
  message("Wrote: ", base_name, "_quality_boxplot.pdf")
}

# Plot 4: Composition bar chart
comp_cols <- c("Replicating", "Normal", "PassedQC-Normal")
if (any(comp_cols %in% names(df))) {
  df_plot <- head(df, 30)
  df_long <- df_plot %>%
    select(Label, any_of(comp_cols)) %>%
    pivot_longer(-Label, names_to = "Category", values_to = "Count")

  if (nrow(df_long) > 0) {
    p <- ggplot(df_long, aes(x = factor(Label, levels = df_plot$Label), y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c(Replicating = "#81C784", Normal = "#F0E68C", `PassedQC-Normal` = "#D1C4E9")) +
      labs(title = "Cell Composition by Sample", x = label_col, y = "Cell Count") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
    ggsave(file.path(out_base, paste0(base_name, "_composition.pdf")), p, width = max(8, nrow(df_plot) * 0.3), height = 6)
    message("Wrote: ", base_name, "_composition.pdf")
  }
}

# Plot 5: Summary by group
if (length(unique(df$Group)) > 1) {
  stats <- df %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean_pct = mean(High_Quality_Pct, na.rm = TRUE),
      sd_pct = sd(High_Quality_Pct, na.rm = TRUE),
      total_seq = sum(`post-scAbsolute`, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(stats, aes(x = Group, y = mean_pct, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct), width = 0.2) +
    geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3) +
    labs(title = paste("Mean (PassedQC+Replicating)/post-scAbsolute % by", group_col), x = group_col, y = "Mean %") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    ylim(0, 100)
  ggsave(file.path(out_base, paste0(base_name, "_summary_by_group.pdf")), p, width = 7, height = 5)
  write.csv(stats, file.path(out_base, paste0(base_name, "_summary_stats.csv")), row.names = FALSE)
  message("Wrote: ", base_name, "_summary_by_group.pdf, _summary_stats.csv")
}

# Plot 6: Metrics heatmap
if (nrow(df) >= 5) {
  heat_cols <- intersect(c("post-scAbsolute", "PassedQC(incl.Normal)", "Replicating", "Normal", "Outliers(Alpha/Mapd/Gini,post-RPC)", "High_Quality_Pct"), names(df))
  if (length(heat_cols) >= 2) {
    df_heat <- df %>%
      select(Label, all_of(heat_cols)) %>%
      head(30) %>%
      pivot_longer(-Label, names_to = "Metric", values_to = "Value") %>%
      group_by(Metric) %>%
      mutate(Scaled = (Value - min(Value, na.rm = TRUE)) / (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE) + 0.001)) %>%
      ungroup()

    p <- ggplot(df_heat, aes(x = Metric, y = factor(Label, levels = rev(unique(Label))), fill = Scaled)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue", midpoint = 0.5) +
      geom_text(aes(label = round(Value, 0)), size = 2) +
      labs(title = "Sample Metrics Overview", x = "Metric", y = label_col) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
    ggsave(file.path(out_base, paste0(base_name, "_metrics_heatmap.pdf")), p, width = 8, height = max(5, length(unique(df_heat$Label)) * 0.25))
    message("Wrote: ", base_name, "_metrics_heatmap.pdf")
  }
}

message("\n=== Complete ===")
message("Samples processed: ", nrow(df))
message("Date: ", run_date)
message("Input: ", input_name)
message("Config: ", config_name)
message("Main output: ", qc_summary_file)
message("Archive: ", basename(archive_dir))
message("================")
