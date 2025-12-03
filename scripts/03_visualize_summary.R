#!/usr/bin/env Rscript

# General visualizer for outlier summary CSVs
# Not specific to any sample type (FFPE, cell lines, etc.)
#
# Usage:
#   Rscript scripts/04_visualize_summary_general.R <input_csv> <out_base> [group_col] [label_col]
#
# Arguments:
#   input_csv  - Path to combined outlier summary CSV
#   out_base   - Output directory for plots and normalized CSV
#   group_col  - (Optional) Column name to use for grouping/coloring (e.g., "feature1", "category")
#   label_col  - (Optional) Column name to use for x-axis labels (e.g., "Sample", "Cell line")
#
# Example:
#   Rscript scripts/04_visualize_summary_general.R results/summary.csv results/ feature1 Sample

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else NULL
out_base  <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "."
group_col <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else NULL
label_col <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else NULL

if (is.null(input_csv) || !file.exists(input_csv)) {
  stop("Usage: Rscript 04_visualize_summary_general.R <input_csv> <out_base> [group_col] [label_col]\n",
       "  input_csv must be provided and exist.")
}

library(ggplot2)
library(tidyr)
library(dplyr)

# ==============================================================================
# Helper functions
# ==============================================================================

#' Pick the first matching column name from candidates
pick_col <- function(df, candidates) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  NULL
}

#' Safe pivot_longer that returns empty tibble if no cols present
safe_pivot <- function(dat, cols_vec) {
  if (!any(cols_vec %in% colnames(dat))) return(tibble::tibble())
  pivot_longer(dat, cols = any_of(cols_vec), names_to = "Metric_Type", values_to = "Count")
}

#' Detect and normalize common column names
normalize_columns <- function(df) {
  # Sample/Label column
  sample_col <- pick_col(df, c("Sample", "sample", "SampleID", "sample_id", "ID", "id"))
  if (!is.null(sample_col)) df$Sample <- as.character(df[[sample_col]])
  
  # Sequenced cells
  seq_col <- pick_col(df, c("Sequenced", "Sequenced.Cells", "sequenced", "Total", "total_cells", "Processed Cells", "Processed.Cells"))
  if (!is.null(seq_col)) {
    df$Sequenced <- suppressWarnings(as.numeric(df[[seq_col]]))
  } else {
    df$Sequenced <- NA_real_
  }
  
  # Good quality cells
  good_col <- pick_col(df, c("Good.Quality.Cells", "Good.Quality", "good_quality", "GoodQuality", "High.Quality", "Good Quality Cells"))
  if (!is.null(good_col)) {
    df$Good_Quality <- suppressWarnings(as.numeric(df[[good_col]]))
  } else {
    df$Good_Quality <- NA_real_
  }
  
  # Replicating cells
  rep_col <- pick_col(df, c("Replicating", "Replicating...RPC", "Replicating..RPC", "replicating", "RPC"))
  if (!is.null(rep_col)) {
    df$Replicating <- suppressWarnings(as.numeric(df[[rep_col]]))
  } else {
    df$Replicating <- NA_real_
  }
  
  # Normal cells
  normal_col <- pick_col(df, c("Normal", "Normal.Cells", "normal", "Normal_Cells", "Normal Cells"))
  if (!is.null(normal_col)) {
    df$Normal <- suppressWarnings(as.numeric(df[[normal_col]]))
  } else {
    df$Normal <- NA_real_
  }
  
  # Outlier cells
  outlier_col <- pick_col(df, c("Outliers", "outliers", "Outlier", "outlier_count", "Alpha/Mapd/Gini", "Alpha.Mapd.Gini"))
  if (!is.null(outlier_col)) {
    df$Outliers <- suppressWarnings(as.numeric(df[[outlier_col]]))
  } else {
    df$Outliers <- NA_real_
  }
  
  df
}

# ==============================================================================
# Load and process data
# ==============================================================================

message("Reading: ", input_csv)
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

# Remove "Total" row if present (summary row)
if ("Sample" %in% colnames(df)) {
  df <- df %>% filter(Sample != "Total" | is.na(Sample))
}

# Normalize column names
df <- normalize_columns(df)

# Auto-detect label column if not specified
if (is.null(label_col)) {
  label_col <- pick_col(df, c("Sample", "Cell line", "Cell_line", "category", "Category", "ID"))
  if (is.null(label_col)) label_col <- "Sample"
}
if (!label_col %in% colnames(df)) {

  df$Label <- if ("Sample" %in% colnames(df)) df$Sample else seq_len(nrow(df))
} else {
  df$Label <- as.character(df[[label_col]])
}

# Auto-detect group column if not specified
if (is.null(group_col)) {
  group_col <- pick_col(df, c("feature1", "feature2", "category", "Category", "Cell line", "Cell_line", "Group"))
}
if (!is.null(group_col) && group_col %in% colnames(df)) {
  df$Group <- as.character(df[[group_col]])
} else {
  df$Group <- "All"
  group_col <- NULL
}

# Compute derived metrics
df <- df %>% mutate(
  # Filtered = Good quality minus normal
  Filtered = coalesce(Good_Quality, 0) - coalesce(Normal, 0),
  # High quality percentage
  High_Quality_Pct = round((coalesce(Good_Quality, 0) + coalesce(Replicating, 0)) * 100 / pmax(Sequenced, 1, na.rm = TRUE)),
  # Pass rate (non-outlier percentage)
  Pass_Rate = round((coalesce(Sequenced, 0) - coalesce(Outliers, 0)) * 100 / pmax(Sequenced, 1, na.rm = TRUE))
)

# Ensure output directory exists
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# Get base name for output files
input_basename <- tools::file_path_sans_ext(basename(input_csv))

# ==============================================================================
# Output: Normalized CSV
# ==============================================================================

out_csv <- file.path(out_base, paste0(input_basename, "_normalized.csv"))
write.csv(df, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)

# ==============================================================================
# Plot 1: Distribution of High Quality Percentage (Density)
# ==============================================================================

if (any(!is.na(df$High_Quality_Pct))) {
  p_density <- ggplot(df, aes(x = High_Quality_Pct)) +
    geom_density(fill = "steelblue", color = "darkblue", alpha = 0.7, linewidth = 0.6) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Distribution of High-Quality Cell Percentage",
      x = "High-Quality (%)",
      y = "Density"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  out_plot <- file.path(out_base, paste0(input_basename, "_quality_distribution.pdf"))
  ggsave(out_plot, p_density, width = 6, height = 4)
  message("Wrote: ", out_plot)
}

# ==============================================================================
# Plot 2: Distribution by Group (if group column exists)
# ==============================================================================

if (!is.null(group_col) && length(unique(df$Group)) > 1 && any(!is.na(df$High_Quality_Pct))) {
  p_group_density <- ggplot(df, aes(x = High_Quality_Pct, fill = Group)) +
    geom_density(alpha = 0.5, linewidth = 0.5) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("High-Quality % Distribution by", group_col),
      x = "High-Quality (%)",
      y = "Density",
      fill = group_col
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  out_plot <- file.path(out_base, paste0(input_basename, "_quality_by_group.pdf"))
  ggsave(out_plot, p_group_density, width = 7, height = 4)
  message("Wrote: ", out_plot)
  
  # Box plot by group
  p_group_box <- ggplot(df, aes(x = Group, y = High_Quality_Pct, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("High-Quality % by", group_col),
      x = group_col,
      y = "High-Quality (%)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  out_plot <- file.path(out_base, paste0(input_basename, "_quality_boxplot.pdf"))
  ggsave(out_plot, p_group_box, width = 7, height = 5)
  message("Wrote: ", out_plot)
}

# ==============================================================================
# Plot 3: Cell Composition Stacked Bar Chart
# ==============================================================================

composition_cols <- c("Replicating", "Normal", "Filtered")
has_composition <- any(composition_cols %in% colnames(df)) && 
                   any(!is.na(df[, intersect(composition_cols, colnames(df)), drop = FALSE]))

if (has_composition) {
  # Limit to reasonable number of samples for readability
  n_samples <- min(30, nrow(df))
  df_plot <- df[seq_len(n_samples), ]
  
  df_long <- safe_pivot(df_plot, composition_cols)
  
  if (nrow(df_long) > 0) {
    # Add Label column to long format
    df_long <- df_long %>%
      mutate(Label = rep(df_plot$Label, length(composition_cols)))
    
    # Color palette
    comp_colors <- c(
      "Replicating" = "#81C784",
      "Normal" = "#F0E68C", 
      "Filtered" = "#D1C4E9"
    )
    
    p_composition <- ggplot(df_long, aes(x = factor(Label, levels = unique(df_plot$Label)), y = Count, fill = Metric_Type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = comp_colors, na.value = "gray70") +
      theme_minimal(base_size = 12) +
      labs(
        title = "Cell Composition by Sample",
        x = label_col,
        y = "Cell Count",
        fill = "Category"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
      )
    
    # Add percentage labels if High_Quality_Pct available
    if (any(!is.na(df_plot$High_Quality_Pct)) && any(!is.na(df_plot$Sequenced))) {
      p_composition <- p_composition +
        geom_text(
          data = df_plot,
          aes(x = factor(Label, levels = unique(Label)), 
              y = Sequenced + max(Sequenced, na.rm = TRUE) * 0.02,
              label = paste0(High_Quality_Pct, "%")),
          inherit.aes = FALSE,
          size = 2.5,
          angle = 90,
          vjust = 0.5,
          hjust = 0
        )
    }
    
    out_plot <- file.path(out_base, paste0(input_basename, "_composition.pdf"))
    ggsave(out_plot, p_composition, width = max(8, n_samples * 0.3), height = 6)
    message("Wrote: ", out_plot)
  }
}

# ==============================================================================
# Plot 4: Summary Statistics Bar Chart
# ==============================================================================

# Create summary by group (or overall if no group)
if (!is.null(group_col) && length(unique(df$Group)) > 1) {
  summary_stats <- df %>%
    group_by(Group) %>%
    summarise(
      n_samples = n(),
      mean_quality_pct = mean(High_Quality_Pct, na.rm = TRUE),
      sd_quality_pct = sd(High_Quality_Pct, na.rm = TRUE),
      median_quality_pct = median(High_Quality_Pct, na.rm = TRUE),
      total_sequenced = sum(Sequenced, na.rm = TRUE),
      total_good_quality = sum(Good_Quality, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_summary <- ggplot(summary_stats, aes(x = Group, y = mean_quality_pct, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_quality_pct - sd_quality_pct, 
                      ymax = mean_quality_pct + sd_quality_pct),
                  width = 0.2) +
    geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.5, size = 3) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Mean High-Quality % by", group_col),
      x = group_col,
      y = "Mean High-Quality (%)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ylim(0, 100)
  
  out_plot <- file.path(out_base, paste0(input_basename, "_summary_by_group.pdf"))
  ggsave(out_plot, p_summary, width = 7, height = 5)
  message("Wrote: ", out_plot)
  
  # Write summary stats CSV
  out_summary_csv <- file.path(out_base, paste0(input_basename, "_summary_stats.csv"))
  write.csv(summary_stats, out_summary_csv, row.names = FALSE)
  message("Wrote: ", out_summary_csv)
}

# ==============================================================================
# Plot 5: Heatmap-style Overview (if many samples)
# ==============================================================================

if (nrow(df) >= 5) {
  # Select numeric columns for heatmap
  heatmap_cols <- c("Sequenced", "Good_Quality", "Replicating", "Normal", "Outliers", "High_Quality_Pct")
  heatmap_cols <- intersect(heatmap_cols, colnames(df))
  
  if (length(heatmap_cols) >= 2) {
    df_heatmap <- df %>%
      select(Label, all_of(heatmap_cols)) %>%
      pivot_longer(cols = -Label, names_to = "Metric", values_to = "Value")
    
    # Scale values within each metric for better visualization
    df_heatmap <- df_heatmap %>%
      group_by(Metric) %>%
      mutate(Value_Scaled = (Value - min(Value, na.rm = TRUE)) / 
               (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE) + 0.001)) %>%
      ungroup()
    
    n_show <- min(30, length(unique(df_heatmap$Label)))
    df_heatmap_sub <- df_heatmap %>% filter(Label %in% unique(df_heatmap$Label)[1:n_show])
    
    p_heatmap <- ggplot(df_heatmap_sub, aes(x = Metric, y = factor(Label, levels = rev(unique(Label))), fill = Value_Scaled)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue", midpoint = 0.5,
                           na.value = "gray90", name = "Scaled\nValue") +
      geom_text(aes(label = round(Value, 0)), size = 2, color = "black") +
      theme_minimal(base_size = 10) +
      labs(
        title = "Sample Metrics Overview",
        x = "Metric",
        y = label_col
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7)
      )
    
    out_plot <- file.path(out_base, paste0(input_basename, "_metrics_heatmap.pdf"))
    ggsave(out_plot, p_heatmap, width = 8, height = max(5, n_show * 0.25))
    message("Wrote: ", out_plot)
  }
}

# ==============================================================================
# Summary message
# ==============================================================================

message("\n=== Visualization Complete ===")
message("Input: ", input_csv)
message("Output directory: ", out_base)
message("Samples processed: ", nrow(df))
if (!is.null(group_col)) message("Grouped by: ", group_col)
message("Label column: ", label_col)
