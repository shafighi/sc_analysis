#!/usr/bin/env Rscript
# ==============================================================================
# Visualize Dropout Heatmap (Publication-Ready)
# ==============================================================================
#
# Generates a compact heatmap of dropout frequency per cell across samples and
# chromosomes. Can be run standalone or as pipeline step 05.
#
# Usage:
#   Rscript scripts/05_visualize_dropout_heatmap.R <samples_csv> <obj_base> [out_base] [bin_size] [group_col]
#
# Arguments:
#   samples_csv - CSV with 'sample' column and metadata (category, etc.)
#   obj_base    - Directory containing step 04 output (cn_binned.rds per sample)
#   out_base    - Output directory for heatmap PDF (default: obj_base)
#   bin_size    - Bin size (default: 100)
#   group_col   - Column from samples_csv to group/facet by (default: "Cell line")
#
# Output:
#   - dropout_heatmap_{date}.pdf : Publication-ready heatmap
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
obj_base    <- if (length(args) >= 2) args[[2]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/analysis_per_sample"
out_base    <- if (length(args) >= 3) args[[3]] else obj_base
bin_size    <- if (length(args) >= 4) args[[4]] else "100"
group_col   <- if (length(args) >= 5) args[[5]] else "Cell line"

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# ==============================================================================
# Load samples
# ==============================================================================
if (is.null(samples_csv) || !file.exists(samples_csv)) {
  stop("Usage: Rscript scripts/05_visualize_dropout_heatmap.R <samples_csv> <obj_base> [out_base] [bin_size] [group_col]")
}

samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")

# Resolve group column
normalize_colname <- function(df, candidates, target) {
  for (c in candidates) if (c %in% names(df)) { names(df)[names(df) == c] <- target; break }
  df
}
samples_tbl <- normalize_colname(samples_tbl, c("Cell line", "Cell_line", "category"), "Category")

if (!group_col %in% names(samples_tbl)) {
  # Try normalized name
  if ("Category" %in% names(samples_tbl)) {
    group_col <- "Category"
  } else {
    message("Group column '", group_col, "' not found. Using ungrouped layout.")
    group_col <- NULL
  }
}

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Collect dropout data from per-sample cn_binned.rds
# ==============================================================================
dropout_summary <- data.frame()

for (i in seq_len(nrow(samples_tbl))) {
  sample_id  <- samples_tbl$sample[i]
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  cn_path    <- file.path(obj_base, sample_obj, "cn_binned.rds")

  if (!file.exists(cn_path)) {
    warning("cn_binned.rds not found for ", sample_id, " - run step 04 first. Skipping.")
    next
  }

  cn_binned <- readRDS(cn_path)
  n_cells   <- length(unique(cn_binned$cell))

  # Count dropout bins per chromosome
  dropout <- cn_binned %>%
    filter(segVal == 0) %>%
    count(chromosome, name = "Freq") %>%
    mutate(
      freq_per_cell = Freq / n_cells,
      sample        = sample_id
    )

  dropout_summary <- bind_rows(dropout_summary, dropout)
}

if (nrow(dropout_summary) == 0) stop("No dropout data found. Run step 04 first.")

# ==============================================================================
# Complete missing chromosomes and attach metadata
# ==============================================================================
all_chroms  <- as.character(1:22)
all_samples <- unique(dropout_summary$sample)

complete_grid <- expand.grid(
  sample     = all_samples,
  chromosome = all_chroms,
  stringsAsFactors = FALSE
)

heatmap_data <- complete_grid %>%
  left_join(dropout_summary, by = c("sample", "chromosome")) %>%
  mutate(
    Freq          = replace_na(Freq, 0),
    freq_per_cell = replace_na(freq_per_cell, 0)
  )

# Attach group metadata from samples_tbl
if (!is.null(group_col)) {
  sample_meta <- samples_tbl %>%
    select(sample, all_of(group_col)) %>%
    distinct()
  heatmap_data <- heatmap_data %>%
    left_join(sample_meta, by = "sample")
}

# Order chromosomes numerically
heatmap_data$chromosome <- factor(heatmap_data$chromosome,
                                  levels = all_chroms)
heatmap_data$sample <- factor(heatmap_data$sample,
                              levels = rev(all_samples))

# Round display values
heatmap_data$display_val <- round(heatmap_data$freq_per_cell, 1)

# ==============================================================================
# Generate publication-ready heatmap
# ==============================================================================

# Use sqrt transform for color scale to handle outlier samples
# This preserves detail in the low range while still showing high values
heatmap_data$fill_val <- sqrt(heatmap_data$freq_per_cell)

# Adaptive text color: white on dark tiles, black on light
# Use a low threshold so text is white on most non-zero tiles
heatmap_data$text_color <- ifelse(
  heatmap_data$freq_per_cell >= 0.3, "white", "grey40"
)

# Smart label: integer for >= 10, one decimal for > 0, blank for 0
heatmap_data$display_val <- ifelse(
  heatmap_data$freq_per_cell == 0, "",
  ifelse(heatmap_data$freq_per_cell >= 10,
         as.character(round(heatmap_data$freq_per_cell)),
         format(round(heatmap_data$freq_per_cell, 1), nsmall = 1))
)

# Breaks for the color bar legend (in original scale)
legend_breaks <- c(0, 5, 20, 50, 100, 200)
legend_breaks <- legend_breaks[legend_breaks <= max(heatmap_data$freq_per_cell, na.rm = TRUE) * 1.1]

# Build plot
p <- ggplot(heatmap_data, aes(x = chromosome, y = sample, fill = fill_val)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = display_val, color = text_color),
            size = 2.0, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_viridis(
    option    = "D",
    name      = "Dropouts\nper cell",
    begin     = 0.02,
    end       = 0.95,
    na.value  = "grey90",
    breaks    = sqrt(legend_breaks),
    labels    = legend_breaks,
    guide     = guide_colorbar(
      barwidth  = 0.5,
      barheight = 3.5,
      title.theme = element_text(size = 6.5, face = "bold"),
      label.theme = element_text(size = 5.5)
    )
  ) +
  labs(x = "Chromosome", y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    panel.grid          = element_blank(),
    panel.background    = element_rect(fill = "white", color = NA),
    plot.background     = element_rect(fill = "white", color = NA),
    panel.spacing       = unit(1, "pt"),
    plot.margin         = margin(2, 2, 2, 2, "mm"),
    axis.text.x         = element_text(size = 6, color = "grey20"),
    axis.text.y         = element_text(size = 6, color = "grey20"),
    axis.title.x        = element_text(size = 7, margin = margin(t = 2)),
    axis.ticks          = element_blank(),
    legend.position     = "right",
    legend.margin       = margin(0, 0, 0, 1, "mm"),
    strip.background    = element_rect(fill = "grey92", color = NA),
    strip.text.y.right  = element_text(size = 6, face = "bold", angle = 0,
                                       margin = margin(0, 1, 0, 1, "mm"))
  )

# Add faceting if group column exists
if (!is.null(group_col) && group_col %in% names(heatmap_data)) {
  p <- p + facet_grid(
    reformulate(".", group_col),
    scales = "free_y",
    space  = "free_y"
  )
}

# Calculate compact dimensions
n_samples <- length(all_samples)
n_groups  <- if (!is.null(group_col) && group_col %in% names(heatmap_data)) {
  length(unique(heatmap_data[[group_col]]))
} else { 0 }

# Find widest strip label to size the right margin
max_label_len <- if (!is.null(group_col) && group_col %in% names(heatmap_data)) {
  max(nchar(as.character(unique(heatmap_data[[group_col]]))), na.rm = TRUE)
} else { 0 }

plot_height <- max(1.5, n_samples * 0.25 + n_groups * 0.08 + 0.5)
plot_width  <- 5.5 + max(0, (max_label_len - 6) * 0.06)

# Save
run_date <- format(Sys.time(), "%Y%m%d")
outfile  <- file.path(out_base, paste0("dropout_heatmap_", run_date, ".pdf"))

ggsave(outfile, p, width = plot_width, height = plot_height,
       units = "in", dpi = 300, device = cairo_pdf)

message("Wrote: ", outfile)
message("Dimensions: ", plot_width, " x ", round(plot_height, 1), " inches")
message("Samples: ", n_samples)
