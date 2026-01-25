#!/usr/bin/env Rscript
# Usage: Rscript scripts/03_visualize_summary.R <input_csv> <out_base> [group_col] [label_col]

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else NULL
out_base  <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "."
group_col <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else NULL
label_col <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else NULL

if (is.null(input_csv) || !file.exists(input_csv)) {
  stop("Usage: Rscript 03_visualize_summary.R <input_csv> <out_base> [group_col] [label_col]")
}

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

# Normalize key columns
col_map <- list(
  Sample      = c("Sample", "sample", "SampleID", "ID"),
  Sequenced   = c("Processed Cells", "Processed.Cells", "Sequenced", "Sequenced.Cells", "Total"),
  Good_Quality = c("Good Quality Cells", "Good.Quality.Cells", "Good.Quality", "GoodQuality"),
  Replicating = c("Replicating", "Replicating...RPC", "Replicating..RPC", "RPC"),
  Normal      = c("Normal Cells", "Normal.Cells", "Normal", "Normal_Cells"),
  Outliers    = c("Alpha/Mapd/Gini", "Alpha.Mapd.Gini", "Outliers")
)
for (target in names(col_map)) {
  src <- pick_col(df, col_map[[target]])
  if (!is.null(src)) df[[target]] <- suppressWarnings(as.numeric(df[[src]]))
  else df[[target]] <- NA_real_
}

# Setup label and group columns
if (is.null(label_col)) label_col <- pick_col(df, c("Sample", "Cell line", "category")) %||% "Sample"
df$Label <- if (label_col %in% names(df)) as.character(df[[label_col]]) else seq_len(nrow(df))

if (is.null(group_col)) group_col <- pick_col(df, c("feature1", "feature2", "category", "Cell line"))
df$Group <- if (!is.null(group_col) && group_col %in% names(df)) as.character(df[[group_col]]) else "All"

# Compute derived metrics
df <- df %>% mutate(
  Filtered = coalesce(Good_Quality, 0) - coalesce(Normal, 0),
  High_Quality_Pct = round((coalesce(Good_Quality, 0) + coalesce(Replicating, 0)) * 100 / pmax(Sequenced, 1)),
  Pass_Rate = round((coalesce(Sequenced, 0) - coalesce(Outliers, 0)) * 100 / pmax(Sequenced, 1))
)

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
base_name <- tools::file_path_sans_ext(basename(input_csv))

# Output normalized CSV
write.csv(df, file.path(out_base, paste0(base_name, "_normalized.csv")), row.names = FALSE)
message("Wrote: ", base_name, "_normalized.csv")

# Plot 1: Quality distribution
if (any(!is.na(df$High_Quality_Pct))) {
  p <- ggplot(df, aes(x = High_Quality_Pct)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    labs(title = "Distribution of High-Quality Cell Percentage", x = "High-Quality (%)", y = "Density") +
    theme_minimal()
  ggsave(file.path(out_base, paste0(base_name, "_quality_distribution.pdf")), p, width = 6, height = 4)
  message("Wrote: ", base_name, "_quality_distribution.pdf")
}

# Plot 2-3: Group comparisons (if multiple groups)
if (length(unique(df$Group)) > 1 && any(!is.na(df$High_Quality_Pct))) {
  # Density by group
  p <- ggplot(df, aes(x = High_Quality_Pct, fill = Group)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("High-Quality % by", group_col), x = "High-Quality (%)", y = "Density") +
    theme_minimal()
  ggsave(file.path(out_base, paste0(base_name, "_quality_by_group.pdf")), p, width = 7, height = 4)
  message("Wrote: ", base_name, "_quality_by_group.pdf")

  # Boxplot by group
  p <- ggplot(df, aes(x = Group, y = High_Quality_Pct, fill = Group)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    labs(title = paste("High-Quality % by", group_col), x = group_col, y = "High-Quality (%)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(out_base, paste0(base_name, "_quality_boxplot.pdf")), p, width = 7, height = 5)
  message("Wrote: ", base_name, "_quality_boxplot.pdf")
}

# Plot 4: Composition bar chart
comp_cols <- c("Replicating", "Normal", "Filtered")
if (any(comp_cols %in% names(df))) {
  df_plot <- head(df, 30)
  df_long <- df_plot %>%
    select(Label, any_of(comp_cols)) %>%
    pivot_longer(-Label, names_to = "Category", values_to = "Count")

  if (nrow(df_long) > 0) {
    p <- ggplot(df_long, aes(x = factor(Label, levels = df_plot$Label), y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c(Replicating = "#81C784", Normal = "#F0E68C", Filtered = "#D1C4E9")) +
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
      total_seq = sum(Sequenced, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(stats, aes(x = Group, y = mean_pct, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct), width = 0.2) +
    geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3) +
    labs(title = paste("Mean High-Quality % by", group_col), x = group_col, y = "Mean %") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    ylim(0, 100)
  ggsave(file.path(out_base, paste0(base_name, "_summary_by_group.pdf")), p, width = 7, height = 5)
  write.csv(stats, file.path(out_base, paste0(base_name, "_summary_stats.csv")), row.names = FALSE)
  message("Wrote: ", base_name, "_summary_by_group.pdf, _summary_stats.csv")
}

# Plot 6: Metrics heatmap
if (nrow(df) >= 5) {
  heat_cols <- intersect(c("Sequenced", "Good_Quality", "Replicating", "Normal", "Outliers", "High_Quality_Pct"), names(df))
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

message("\n=== Complete: ", nrow(df), " samples processed ===")
