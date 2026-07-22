#!/usr/bin/env Rscript
# ==============================================================================
# Post-scUnique visualizations
# ==============================================================================
#
# Generates, from one completed scUnique result directory:
#   1. Copy-number heatmap ordered by the scUnique evolutionary tree.
#   2. Distribution and per-cell counts of validated unique events (freq == 1).
#   3. CSV containing the per-cell frequency-1 event counts.
#
# Usage:
#   Rscript scripts/generate_visualizations.R \
#     <scunique_result_dir> [output_dir] [prefix]
#
# The result directory must contain:
#   <prefix>.finalCN.RDS
#   <prefix>.tree.RDS
#   <prefix>.df_pass_post.RDS
#
# If prefix is omitted, it is inferred when exactly one *.finalCN.RDS exists.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L) {
  stop(
    "Usage: Rscript scripts/generate_visualizations.R ",
    "<scunique_result_dir> [output_dir] [prefix]"
  )
}

result_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(result_dir, "post_scunique")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

infer_prefix <- function(path) {
  candidates <- list.files(path, pattern = "\\.finalCN\\.RDS$", full.names = FALSE)
  if (length(candidates) != 1L) {
    stop(
      "Expected exactly one *.finalCN.RDS in ", path,
      "; found ", length(candidates),
      ". Supply prefix as the third argument."
    )
  }
  sub("\\.finalCN\\.RDS$", "", candidates[[1]])
}

prefix <- if (length(args) >= 3L) args[[3]] else infer_prefix(result_dir)
base <- file.path(result_dir, prefix)

suppressPackageStartupMessages({
  library(Biobase)
  library(QDNAseq)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

read_output <- function(suffix) {
  path <- paste0(base, ".", suffix, ".RDS")
  if (!file.exists(path)) stop("Missing scUnique output: ", path)
  readRDS(path)
}

common_prefix <- function(values) {
  values <- as.character(values)
  if (!length(values) || length(values) == 1L) return("")

  limit <- min(nchar(values))
  if (limit == 0L) return("")

  same <- vapply(
    seq_len(limit),
    function(i) length(unique(substr(values, i, i))) == 1L,
    logical(1)
  )
  first_difference <- which(!same)[1]
  raw_prefix <- if (is.na(first_difference)) {
    substr(values[[1]], 1L, limit)
  } else if (first_difference == 1L) {
    ""
  } else {
    substr(values[[1]], 1L, first_difference - 1L)
  }

  # Avoid removing a partial token such as the first digit of a plate number.
  sub("[^_:/-]*$", "", raw_prefix)
}

shorten_cell_names <- function(values) {
  prefix_removed <- common_prefix(values)
  labels <- if (nchar(prefix_removed)) {
    substring(values, nchar(prefix_removed) + 1L)
  } else {
    values
  }

  if (any(!nzchar(labels)) || anyDuplicated(labels)) {
    warning("Shared-prefix removal produced empty or duplicate labels; using full names")
    prefix_removed <- ""
    labels <- values
  }

  list(labels = labels, prefix = prefix_removed)
}

cn <- read_output("finalCN")
tree <- read_output("tree")
events <- read_output("df_pass_post")
dend <- if (inherits(tree, "dendrogram")) tree else as.dendrogram(tree)

tree_cells <- labels(dend)
cn_cells <- colnames(cn)
missing_from_cn <- setdiff(tree_cells, cn_cells)
missing_from_tree <- setdiff(cn_cells, tree_cells)
if (length(missing_from_cn) || length(missing_from_tree)) {
  stop(
    "Tree/CN cell mismatch. Missing from CN: ",
    paste(missing_from_cn, collapse = ", "),
    "; missing from tree: ",
    paste(missing_from_tree, collapse = ", ")
  )
}

label_info <- shorten_cell_names(tree_cells)
cell_labels <- label_info$labels
names(cell_labels) <- tree_cells

cn_ordered <- cn[, tree_cells]
cn_matrix <- t(round(Biobase::assayDataElement(cn_ordered, "copynumber")))
rownames(cn_matrix) <- tree_cells

chromosome <- factor(
  vapply(strsplit(colnames(cn_matrix), ":", fixed = TRUE), `[[`, character(1), 1L),
  levels = c(as.character(seq_len(22)), "X", "Y")
)
keep <- !is.na(chromosome)
cn_matrix <- cn_matrix[, keep, drop = FALSE]
chromosome <- droplevels(chromosome[keep])

cutoff <- 10L
cn_plot <- cn_matrix
cn_plot[cn_plot > cutoff] <- cutoff + 1L
observed_states <- sort(unique(as.integer(cn_plot[!is.na(cn_plot)])))
legend_at <- as.character(observed_states)
legend_labels <- ifelse(observed_states == cutoff + 1L,
                        paste0(">", cutoff), as.character(observed_states))

cn_colors <- c(
  "0" = "#3182BD", "1" = "#9ECAE1", "2" = "#CCCCCC",
  "3" = "#FDCC8A", "4" = "#FC8D59", "5" = "#E34A33",
  "6" = "#B30000", "7" = "#980043", "8" = "#DD1C77",
  "9" = "#DF65B0", "10" = "#C994C7", "11" = "#D4B9DA"
)

# Force categorical copy-number states. Older ComplexHeatmap releases may
# otherwise interpret a numeric matrix with a named palette as continuous and
# fail specifically on state 0.
cn_plot <- matrix(
  as.character(cn_plot),
  nrow = nrow(cn_plot),
  ncol = ncol(cn_plot),
  dimnames = dimnames(cn_plot)
)


chromosome_levels <- levels(chromosome)[levels(chromosome) %in% as.character(chromosome)]
chromosome_annotation <- ComplexHeatmap::HeatmapAnnotation(
  chromosome = ComplexHeatmap::anno_block(
    labels = chromosome_levels,
    gp = grid::gpar(fill = NA, col = NA),
    labels_gp = grid::gpar(fontsize = 8)
  ),
  show_annotation_name = FALSE
)

heatmap <- ComplexHeatmap::Heatmap(
  cn_plot,
  name = "Copy number",
  col = cn_colors,
  na_col = "#000000",
  cluster_rows = dend,
  cluster_columns = FALSE,
  column_split = chromosome,
  bottom_annotation = chromosome_annotation,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_labels = unname(cell_labels[rownames(cn_plot)]),
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 6.2),
  row_names_max_width = grid::unit(34, "mm"),
  row_dend_side = "left",
  row_dend_width = grid::unit(65, "mm"),
  row_dend_gp = grid::gpar(col = "#303030", lwd = 1.1),
  column_title = "Chromosome",
  column_title_side = "bottom",
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 10,
  heatmap_legend_param = list(
    at = legend_at,
    labels = legend_labels,
    title_gp = grid::gpar(fontsize = 10),
    labels_gp = grid::gpar(fontsize = 8)
  )
)

heatmap_pdf <- file.path(output_dir, paste0(prefix, "_tree_cn_heatmap_labeled.pdf"))
pdf(heatmap_pdf, width = 20, height = 24, useDingbats = FALSE)
ComplexHeatmap::draw(heatmap)
invisible(dev.off())

heatmap_png <- file.path(output_dir, paste0(prefix, "_tree_cn_heatmap_labeled.png"))
png(heatmap_png, width = 6000, height = 7200, res = 300, type = "cairo")
ComplexHeatmap::draw(heatmap)
invisible(dev.off())

required_columns <- c("cellname", "freq")
if (!all(required_columns %in% names(events))) {
  stop(
    "df_pass_post is missing required columns: ",
    paste(setdiff(required_columns, names(events)), collapse = ", ")
  )
}

freq1 <- events[!is.na(events$freq) & as.numeric(events$freq) == 1, , drop = FALSE]
counts <- table(factor(freq1$cellname, levels = cn_cells))
event_counts <- data.frame(
  cellname = cn_cells,
  cell_label = unname(cell_labels[cn_cells]),
  n_unique_events = as.integer(counts),
  stringsAsFactors = FALSE
)
event_counts <- event_counts[order(-event_counts$n_unique_events, event_counts$cellname), ]

counts_csv <- file.path(output_dir, paste0(prefix, "_freq1_unique_events_per_cell.csv"))
write.csv(event_counts, counts_csv, row.names = FALSE, quote = TRUE)

summary_text <- sprintf(
  "n = %d cells; total events = %d; median = %.1f; mean = %.2f; zero-event cells = %d",
  nrow(event_counts),
  sum(event_counts$n_unique_events),
  median(event_counts$n_unique_events),
  mean(event_counts$n_unique_events),
  sum(event_counts$n_unique_events == 0L)
)

max_count <- max(event_counts$n_unique_events)
axis_step <- max(1L, ceiling(max_count / 12))

p_hist <- ggplot(event_counts, aes(x = n_unique_events)) +
  geom_histogram(
    breaks = seq(-0.5, max_count + 0.5, by = 1),
    fill = "#377EB8",
    color = "white",
    linewidth = 0.25
  ) +
  scale_x_continuous(breaks = seq(0, max_count, by = axis_step)) +
  labs(
    title = "Private (freq = 1) rCNA burden per cell",
    subtitle = summary_text,
    x = "Private events per cell",
    y = "Cells with this many private events"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

plot_counts <- event_counts[order(event_counts$n_unique_events, event_counts$cellname), ]
plot_counts$cell_label <- factor(plot_counts$cell_label, levels = plot_counts$cell_label)

p_cell <- ggplot(plot_counts, aes(x = n_unique_events, y = cell_label)) +
  geom_col(fill = "#377EB8", width = 0.72) +
  geom_text(
    aes(label = n_unique_events),
    hjust = -0.25,
    size = 1.9,
    color = "#222222"
  ) +
  scale_x_continuous(
    breaks = seq(0, max_count, by = axis_step),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = "Private events per cell (freq = 1)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 5.5)
  )

event_plot <- cowplot::plot_grid(
  p_hist,
  p_cell,
  ncol = 1,
  rel_heights = c(1, 5),
  align = "v",
  axis = "lr"
)

events_pdf <- file.path(output_dir, paste0(prefix, "_freq1_unique_events_distribution.pdf"))
ggsave(events_pdf, event_plot, width = 11, height = 26, units = "in", device = cairo_pdf)

events_png <- file.path(output_dir, paste0(prefix, "_freq1_unique_events_distribution.png"))
ggsave(events_png, event_plot, width = 11, height = 26, units = "in", dpi = 300, limitsize = FALSE)

run_summary <- data.frame(
  metric = c(
    "prefix_removed_from_cell_labels", "cells", "bins", "tree_labels",
    "validated_freq1_events", "cells_with_freq1_events", "zero_event_cells",
    "median_events_per_cell", "mean_events_per_cell", "max_events_per_cell"
  ),
  value = c(
    label_info$prefix, ncol(cn), nrow(cn), length(tree_cells),
    sum(event_counts$n_unique_events), sum(event_counts$n_unique_events > 0L),
    sum(event_counts$n_unique_events == 0L), median(event_counts$n_unique_events),
    sprintf("%.4f", mean(event_counts$n_unique_events)), max_count
  ),
  stringsAsFactors = FALSE
)
write.csv(
  run_summary,
  file.path(output_dir, paste0(prefix, "_post_scunique_summary.csv")),
  row.names = FALSE,
  quote = TRUE
)

cat("prefix=", prefix, "\n", sep = "")
cat("prefix_removed_from_labels=", label_info$prefix, "\n", sep = "")
cat("cells=", ncol(cn), "\n", sep = "")
cat("bins=", nrow(cn), "\n", sep = "")
cat("tree_labels=", length(tree_cells), "\n", sep = "")
cat("validated_freq1_events=", sum(event_counts$n_unique_events), "\n", sep = "")
cat("cells_with_freq1_events=", sum(event_counts$n_unique_events > 0L), "\n", sep = "")
cat("zero_event_cells=", sum(event_counts$n_unique_events == 0L), "\n", sep = "")
cat("median_events_per_cell=", median(event_counts$n_unique_events), "\n", sep = "")
cat("mean_events_per_cell=", sprintf("%.4f", mean(event_counts$n_unique_events)), "\n", sep = "")
cat("max_events_per_cell=", max_count, "\n", sep = "")
cat("output_dir=", normalizePath(output_dir), "\n", sep = "")
