#!/usr/bin/env Rscript
# ==============================================================================
# Visualize QC Summary
# ==============================================================================
#
# Usage: Rscript scripts/03_visualize_summary.R <input_csv> <out_base> [group_col] [label_col]
#
# Reads CSV with metadata header (# comments) and generates visualization plots.
# Does NOT create additional CSV files - the input CSV from step 2 is the final output.
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

library(ggplot2)
library(tidyr)
library(dplyr)

# Helper: find first matching column
pick_col <- function(df, candidates) {
  for (c in candidates) if (c %in% names(df)) return(c)
  NULL
}

# Read CSV, skipping comment lines
message("Reading: ", input_csv)
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
df <- df %>% filter(is.na(Sample) | Sample != "Total")

# Normalize key columns for plotting
col_map <- list(
  Sample      = c("Sample", "sample", "SampleID", "ID"),
  `post-scAbsolute` = c("Processed Cells", "Processed.Cells", "post-scAbsolute"),
  Good_Quality = c("PassedQC(incl.Normal)", "PassedQC.incl.Normal.", "Good Quality Cells", "Good.Quality.Cells", "Good.Quality"),
  Replicating = c("Replicating", "Replicating...RPC"),
  Normal      = c("Normal Cells", "Normal.Cells", "Normal"),
  `Outliers(RPC)` = c("Outliers(RPC)", "Outliers.RPC."),
  `Outliers(Alpha/Mapd/Gini)` = c("Outliers(Alpha/Mapd/Gini,post-RPC)", "Outliers.Alpha.Mapd.Gini.post.RPC.")
)
for (target in names(col_map)) {
  src <- pick_col(df, col_map[[target]])
  if (!is.null(src)) df[[target]] <- suppressWarnings(as.numeric(df[[src]]))
  else df[[target]] <- NA_real_
}
# Derive PassedQC (excluding Normal) for composition chart
df$PassedQC <- coalesce(df$Good_Quality, 0) - coalesce(df$Normal, 0)

# Setup label and group columns
if (is.null(label_col)) label_col <- pick_col(df, c("Sample", "Cell line", "category")) %||% "Sample"
label_values <- if (label_col %in% names(df)) as.character(df[[label_col]]) else seq_len(nrow(df))

if (is.null(group_col)) group_col <- pick_col(df, c("feature1", "feature2", "category", "Cell line"))
group_values <- if (!is.null(group_col) && group_col %in% names(df)) as.character(df[[group_col]]) else "All"

# Compute high quality percentage for plotting
df <- df %>% mutate(
  High_Quality_Pct = round((coalesce(Good_Quality, 0) + coalesce(Replicating, 0)) * 100 / pmax(`post-scAbsolute`, 1))
)

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
base_name <- tools::file_path_sans_ext(basename(input_csv))

# Add plotting columns
df$Label <- label_values
df$Group <- group_values

# Generate plots
# Plot 1: Quality distribution
if (any(!is.na(df$High_Quality_Pct))) {
  p <- ggplot(df, aes(x = High_Quality_Pct)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    labs(title = "Distribution of (PassedQC+Replicating)/post-scAbsolute %", x = "%", y = "Density") +
    theme_minimal()
  ggsave(file.path(out_base, paste0(base_name, "_quality_distribution.pdf")), p, width = 6, height = 4)
  message("Wrote: ", base_name, "_quality_distribution.pdf")
}

# Plot 2-3: Group comparisons
if (length(unique(df$Group)) > 1 && any(!is.na(df$High_Quality_Pct))) {
  p <- ggplot(df, aes(x = High_Quality_Pct, fill = Group)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Quality by", group_col), x = "%", y = "Density") +
    theme_minimal()
  ggsave(file.path(out_base, paste0(base_name, "_quality_by_group.pdf")), p, width = 7, height = 4)

  p <- ggplot(df, aes(x = Group, y = High_Quality_Pct, fill = Group)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    labs(title = paste("Quality by", group_col), x = group_col, y = "%") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(out_base, paste0(base_name, "_quality_boxplot.pdf")), p, width = 7, height = 5)
  message("Wrote: boxplot and group plots")
}

# Plot 4: Composition bar chart (categories matching outlier_summary.csv)
comp_cols <- c("Replicating", "Outliers(RPC)", "Outliers(Alpha/Mapd/Gini)", "PassedQC", "Normal")
# Okabe-Ito colorblind-safe palette (Nature recommended)
comp_colors <- c(
  "Replicating"                 = "#E69F00",
  "Outliers(RPC)"               = "#D55E00",
  "Outliers(Alpha/Mapd/Gini)"   = "#0072B2",
  "PassedQC"                    = "#009E73",
  "Normal"                      = "#CC79A7"
)

if (any(comp_cols %in% names(df))) {
  df_plot <- head(df, 30)
  df_long <- df_plot %>%
    select(Label, any_of(comp_cols)) %>%
    pivot_longer(-Label, names_to = "Category", values_to = "Count") %>%
    mutate(Category = factor(Category, levels = rev(names(comp_colors))))

  if (nrow(df_long) > 0) {
    p <- ggplot(df_long, aes(x = factor(Label, levels = df_plot$Label), y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = comp_colors, drop = FALSE) +
      labs(title = "Cell Composition by Sample", x = label_col, y = "Cell Count") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
    ggsave(file.path(out_base, paste0(base_name, "_composition.pdf")), p, width = max(8, nrow(df_plot) * 0.3), height = 6)
    message("Wrote: ", base_name, "_composition.pdf")
  }
}

message("\n=== Complete ===")
message("Samples: ", nrow(df))
message("Plots saved to: ", out_base)
message("===============")
