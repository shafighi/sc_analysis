#!/usr/bin/env Rscript
# ==============================================================================
# Visualize Per-Sample Dropout Heatmaps (Publication-Ready)
# ==============================================================================
#
# For each sample, generates compact, publication-ready plots:
#   - Dropout frequency barplot by chromosome (PDF)
#   - Cell x Chromosome dropout heatmap, faceted by QC status (PDF)
#   - Well-plate layout heatmap of total dropouts per cell (PNG)
#
# Requires step 04 output (cn_binned.rds, dropout_per_cell_chr.rds).
# Optionally reads cell_pos_dropout.rds for plate layout; if absent,
# derives well positions from cell names (e.g. "..._A12" -> row A, col 12).
#
# Usage:
#   Rscript scripts/06_visualize_sample_dropouts.R <samples_csv> <obj_base> [out_base] [bin_size]
#
# Arguments:
#   samples_csv - CSV with 'sample' column (and optional metadata)
#   obj_base    - Directory containing step 04 per-sample output
#   out_base    - Output directory for plots (default: obj_base)
#   bin_size    - Bin size (default: 100)
#
# Per-sample output (in out_base/SLX-{sample}_{bin_size}/figures/):
#   - dropout_barplot.pdf       : Chromosome-level dropout frequency barplot
#   - dropout_heatmap.pdf       : Cell x Chromosome dropout heatmap by QC status
#   - dropout_plate_heatmap.png : Well-plate layout heatmap of total dropouts
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript scripts/06_visualize_sample_dropouts.R <samples_csv> <obj_base> [out_base] [bin_size]")
}
samples_csv <- args[[1]]
obj_base    <- args[[2]]
out_base    <- if (length(args) >= 3) args[[3]] else obj_base
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
source("R/visualization_helpers.R")
source("R/summary_helpers.R")

# ==============================================================================
# Load samples
# ==============================================================================
if (!file.exists(samples_csv)) {
  stop("Samples CSV does not exist: ", samples_csv)
}

samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")

message("Processing ", nrow(samples_tbl), " samples from ", samples_csv)

# ==============================================================================
# Shared theme for publication-quality plots
# ==============================================================================
theme_pub <- function(base_size = 7) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm"),
      axis.ticks        = element_line(color = "grey70", linewidth = 0.3),
      axis.ticks.length = unit(1, "mm"),
      axis.text         = element_text(size = 6, color = "grey20"),
      axis.title        = element_text(size = 7),
      plot.title        = element_text(size = 8, face = "bold", hjust = 0,
                                       margin = margin(b = 2)),
      plot.subtitle     = element_text(size = 6.5, color = "grey40",
                                       margin = margin(b = 2)),
      legend.title      = element_text(size = 6.5, face = "bold"),
      legend.text       = element_text(size = 5.5),
      strip.background  = element_rect(fill = "grey92", color = NA),
      strip.text        = element_text(size = 6.5, face = "bold",
                                       margin = margin(1.5, 2, 1.5, 2, "mm"))
    )
}

# ==============================================================================
# Plot functions
# ==============================================================================

#' Dropout frequency barplot per chromosome
plot_dropout_barplot <- function(cn_binned, sample_id, outfile) {
  n_cells <- length(unique(cn_binned$cell))

  dropout <- cn_binned %>%
    filter(segVal == 0) %>%
    count(chromosome, name = "freq") %>%
    mutate(freq_per_cell = freq / n_cells)

  # Complete missing chromosomes
  all_chroms <- as.character(1:22)
  dropout <- data.frame(chromosome = all_chroms, stringsAsFactors = FALSE) %>%
    left_join(dropout, by = "chromosome") %>%
    mutate(
      freq         = replace_na(freq, 0),
      freq_per_cell = replace_na(freq_per_cell, 0),
      chromosome   = factor(chromosome, levels = all_chroms)
    )

  p <- ggplot(dropout, aes(x = chromosome, y = freq_per_cell)) +
    geom_col(fill = "#4A90D9", width = 0.7) +
    geom_text(aes(label = ifelse(freq_per_cell >= 0.5,
                                  round(freq_per_cell, 1), "")),
              vjust = -0.3, size = 1.8, color = "grey30") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = paste0("SLX-", sample_id),
      subtitle = paste0("Dropout bins per cell (", n_cells, " cells)"),
      x = "Chromosome", y = "Dropouts / cell"
    ) +
    theme_pub()

  ggsave(outfile, p, width = 4.5, height = 2, units = "in",
         dpi = 300, device = cairo_pdf)
}

#' Cell x Chromosome dropout heatmap faceted by QC status
plot_dropout_heatmap <- function(cell_chr_dropout, sample_id, outfile) {
  df <- cell_chr_dropout

  # Order chromosomes numerically
  all_chroms <- as.character(1:22)
  df$chromosome <- factor(df$chromosome, levels = all_chroms)
  df <- df %>% filter(!is.na(chromosome))

  # Abbreviate cell names: extract well position (e.g. A10) from cell name
  m <- regexpr("[A-P][0-9]+", df$cell)
  df$cell_label <- ifelse(m > 0, regmatches(df$cell, m), sub("^.*_", "", df$cell))
  # If extraction failed or labels aren't unique, fall back to longer suffix
  if (length(unique(df$cell_label)) < length(unique(df$cell))) {
    df$cell_label <- sub("^.*Plate_", "", df$cell)
  }
  # Final fallback: make unique
  cell_lookup <- data.frame(
    cell = unique(df$cell),
    cell_label = make.unique(
      df$cell_label[match(unique(df$cell), df$cell)]
    ),
    stringsAsFactors = FALSE
  )
  df$cell_label <- cell_lookup$cell_label[match(df$cell, cell_lookup$cell)]

  # Order cells: within each status group, sort by total dropout (descending)
  cell_order <- df %>%
    group_by(cell, cell_label, status) %>%
    summarise(total = sum(dropout_count, na.rm = TRUE), .groups = "drop") %>%
    arrange(status, desc(total))
  df$cell_label <- factor(df$cell_label, levels = unique(cell_order$cell_label))

  n_cells  <- length(unique(df$cell))
  n_groups <- length(unique(df$status))
  max_val  <- max(df$dropout_count, na.rm = TRUE)

  # Adaptive text: show numbers only when cells <= 60 and numbers > 0
  show_text <- n_cells <= 60
  if (show_text) {
    df$display_val <- ifelse(df$dropout_count == 0, "",
                             as.character(df$dropout_count))
  }

  # Adaptive font size based on number of cells
  text_size <- if (n_cells <= 20) 2.2 else if (n_cells <= 40) 1.6 else 1.2

  # Sqrt-transform for better visual spread
  df$fill_val <- sqrt(df$dropout_count)

  # Adaptive text color: white on dark tiles, dark on light tiles
  df$text_color <- ifelse(df$dropout_count >= max(1, max_val * 0.2),
                          "white", "grey30")

  # Build color scale breaks in original scale
  breaks_orig <- c(0, 2, 5, 10, 25, 50, 100)
  breaks_orig <- breaks_orig[breaks_orig <= max_val * 1.1]
  if (length(breaks_orig) < 2) breaks_orig <- c(0, max(1, max_val))

  p <- ggplot(df, aes(x = chromosome, y = cell_label, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis(
      option   = "D",
      name     = "Dropout\nbins",
      begin    = 0.02,
      end      = 0.92,
      na.value = "grey90",
      breaks   = sqrt(breaks_orig),
      labels   = breaks_orig,
      guide    = guide_colorbar(
        barwidth  = 0.5,
        barheight = 3,
        title.theme = element_text(size = 6, face = "bold"),
        label.theme = element_text(size = 5)
      )
    ) +
    scale_color_identity() +
    facet_grid(status ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = paste0("SLX-", sample_id),
      subtitle = paste0("Dropout bins per cell per chromosome (", n_cells, " cells)"),
      x = "Chromosome", y = NULL
    ) +
    theme_pub() +
    theme(
      axis.text.y        = element_text(size = max(3, min(6, 120 / n_cells)),
                                        color = "grey20"),
      panel.spacing      = unit(2, "mm"),
      strip.text.y.right = element_text(angle = 0)
    )

  if (show_text) {
    p <- p + geom_text(aes(label = display_val, color = text_color),
                        size = text_size, fontface = "bold",
                        show.legend = FALSE)
  }

  # Adaptive dimensions: compact but readable
  row_height <- if (n_cells <= 20) 0.22 else if (n_cells <= 60) 0.16 else 0.10
  plot_height <- max(2, n_cells * row_height + n_groups * 0.3 + 0.8)
  plot_height <- min(plot_height, 20)  # cap at 20 inches
  plot_width  <- 5.5

  ggsave(outfile, p, width = plot_width, height = plot_height,
         units = "in", dpi = 300, device = cairo_pdf)
}

#' Render a single plate heatmap to PNG
render_plate_heatmap <- function(plate, sample_id, plate_label, n_cells,
                                  outfile) {
  all_rows <- sort(unique(plate$AZ))
  all_cols <- sort(unique(plate$Number))

  # Complete the grid so empty wells show as grey
  full_grid <- expand.grid(
    AZ = all_rows, Number = all_cols, stringsAsFactors = FALSE
  )
  plate <- full_grid %>%
    left_join(plate, by = c("AZ", "Number"))

  # Order: rows A-P top to bottom, columns numerically left to right
  row_letters <- LETTERS[1:16]
  plate$AZ     <- factor(plate$AZ, levels = rev(intersect(row_letters, all_rows)))
  plate$Number <- factor(plate$Number, levels = all_cols)

  max_val <- max(plate$total_dropout, na.rm = TRUE)

  plate$fill_val    <- sqrt(plate$total_dropout)
  plate$display_val <- ifelse(is.na(plate$total_dropout), "",
                              as.character(plate$total_dropout))
  # Inferno (option "C") is dark at low values and bright yellow at high values,
  # so LOW dropout tiles are dark and need white text; HIGH values are bright and
  # need dark text.
  mid_thresh <- max(1, max_val * 0.35)
  plate$text_color <- ifelse(
    is.na(plate$total_dropout), "grey70",
    ifelse(plate$total_dropout < mid_thresh, "white", "grey20")
  )

  n_cols <- length(all_cols)
  n_rows <- length(all_rows)
  text_size <- if (n_cols <= 12) 2.8 else if (n_cols <= 18) 2.2 else 1.9

  breaks_orig <- c(0, 10, 25, 50, 100, 200, 500)
  breaks_orig <- breaks_orig[breaks_orig <= max_val * 1.1]
  if (length(breaks_orig) < 2) breaks_orig <- c(0, max(1, max_val))

  p <- ggplot(plate, aes(x = Number, y = AZ, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = display_val, color = text_color),
              size = text_size, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_viridis(
      option = "C", name = "Total\ndropout\nbins",
      begin = 0.05, end = 0.95, na.value = "grey93",
      breaks = sqrt(breaks_orig), labels = breaks_orig,
      guide = guide_colorbar(
        barwidth = 0.6, barheight = 4,
        title.theme = element_text(size = 7, face = "bold"),
        label.theme = element_text(size = 6)
      )
    ) +
    labs(
      title    = paste0("SLX-", sample_id, "  ", plate_label),
      subtitle = paste0("Total dropout bins per cell (", n_cells, " cells)"),
      x = "Column", y = "Row"
    ) +
    theme_pub(base_size = 8) +
    theme(
      axis.text.x  = element_text(size = 7),
      axis.text.y  = element_text(size = 7),
      aspect.ratio = n_rows / n_cols
    )

  plot_width  <- max(4, n_cols * 0.32 + 1.5)
  plot_height <- max(2.5, n_rows * 0.32 + 1.2)
  ggsave(outfile, p, width = plot_width, height = plot_height,
         units = "in", dpi = 300)
}

#' Well-plate layout heatmap of total dropouts per cell
#'
#' Detects plates from cell names (e.g. "...Plate_537_A10_1") and generates
#' one PNG per plate. Falls back to cell_pos_dropout.rds for single-plate
#' samples when available.
plot_dropout_plate <- function(cell_pos_dropout_path, cell_chr_dropout,
                               sample_id, fig_dir) {
  # Compute total dropout per cell
  cell_totals <- cell_chr_dropout %>%
    group_by(cell) %>%
    summarise(total_dropout = sum(dropout_count, na.rm = TRUE), .groups = "drop")

  # Extract plate ID and well position from cell names
  # Pattern: "Plate_<id>_<well>" where well is like A10
  plate_match <- regexpr("Plate_[0-9]+", cell_totals$cell)
  cell_totals$plate_id <- ifelse(plate_match > 0,
                                  regmatches(cell_totals$cell, plate_match),
                                  NA)
  well_match <- regexpr("[A-P][0-9]+", cell_totals$cell)
  cell_totals$well <- ifelse(well_match > 0,
                              regmatches(cell_totals$cell, well_match), NA)
  cell_totals$AZ     <- sub("[0-9]+$", "", cell_totals$well)
  cell_totals$Number <- as.integer(sub("^[A-Z]+", "", cell_totals$well))

  valid <- cell_totals %>%
    filter(grepl("^[A-P]$", AZ), !is.na(Number))

  # If no valid well positions from cell names, try cell_pos_dropout.rds
  if (nrow(valid) == 0 && file.exists(cell_pos_dropout_path)) {
    rds <- readRDS(cell_pos_dropout_path)
    rds$AZ     <- as.character(rds$AZ)
    rds$Number <- as.integer(as.character(rds$Number))
    rds <- rds %>% filter(!is.na(AZ), !is.na(Number), nchar(AZ) == 1)
    if (nrow(rds) > 1) {
      rds$plate_id <- "Plate"
      valid <- rds %>% select(plate_id, AZ, Number, total_dropout)
    }
  }

  if (nrow(valid) == 0) {
    warning("Could not get well positions for ", sample_id,
            " - skipping plate heatmap.")
    return(invisible(NULL))
  }

  # Split by plate and render each
  plates <- sort(unique(valid$plate_id))

  for (pid in plates) {
    plate_data <- valid %>%
      filter(plate_id == pid) %>%
      select(AZ, Number, total_dropout)
    n_cells <- nrow(plate_data)

    if (length(plates) == 1) {
      outfile <- file.path(fig_dir, "dropout_plate_heatmap.png")
      label <- ""
    } else {
      safe_pid <- gsub("[^A-Za-z0-9_]", "", pid)
      outfile <- file.path(fig_dir, paste0("dropout_plate_heatmap_",
                                            safe_pid, ".png"))
      label <- paste0("(", pid, ")")
    }

    render_plate_heatmap(plate_data, sample_id, label, n_cells, outfile)
    message("  Plate heatmap: ", basename(outfile))
  }
}

# ==============================================================================
# Process each sample
# ==============================================================================
for (i in seq_len(nrow(samples_tbl))) {
  sample_id  <- samples_tbl$sample[i]
  sample_obj <- paste0("SLX-", sample_id, "_", bin_size)
  sample_dir <- file.path(obj_base, sample_obj)
  out_dir    <- file.path(out_base, sample_obj)

  cn_binned_path     <- file.path(sample_dir, "cn_binned.rds")
  cell_chr_drop_path <- file.path(sample_dir, "dropout_per_cell_chr.rds")

  if (!file.exists(cn_binned_path) || !file.exists(cell_chr_drop_path)) {
    warning("Step 04 output not found for ", sample_id, " - skipping.")
    next
  }

  message("Plotting: ", sample_id)

  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  cn_binned        <- readRDS(cn_binned_path)
  cell_chr_dropout <- readRDS(cell_chr_drop_path)

  # Barplot
  plot_dropout_barplot(
    cn_binned, sample_id,
    file.path(fig_dir, "dropout_barplot.pdf")
  )

  # Heatmap
  plot_dropout_heatmap(
    cell_chr_dropout, sample_id,
    file.path(fig_dir, "dropout_heatmap.pdf")
  )

  # Well-plate layout heatmap (one per plate)
  cell_pos_path <- file.path(sample_dir, "cell_pos_dropout.rds")
  plot_dropout_plate(cell_pos_path, cell_chr_dropout, sample_id, fig_dir)

  # Regenerate QC summary panel with CN heatmap + dropout plate heatmap
  qc_metric_path    <- file.path(sample_dir, "qc_metric_data.rds")
  plate_heatmap_png <- file.path(fig_dir, "dropout_plate_heatmap.png")
  if (file.exists(qc_metric_path)) {
    qc_data <- tryCatch(readRDS(qc_metric_path), error = function(e) NULL)
    if (!is.null(qc_data)) {
      cn_png <- if (!is.null(qc_data$cn_heatmap_png) && file.exists(qc_data$cn_heatmap_png))
                  qc_data$cn_heatmap_png else NULL
      plate_png <- if (file.exists(plate_heatmap_png)) plate_heatmap_png else NULL
      tryCatch(
        plot_qc_summary_panel(
          df                 = qc_data$df,
          iq                 = qc_data$iq,
          summary_df         = qc_data$summary_df,
          sample_id          = sample_id,
          outfile            = file.path(fig_dir, "qc_summary_panel.pdf"),
          cn_heatmap_path    = cn_png,
          plate_heatmap_path = plate_png
        ),
        error = function(e) warning("QC panel update failed for ", sample_id, ": ", e$message)
      )
      message("  Updated QC summary panel with CN heatmap + dropout plate.")
    }
  }

  # Per-condition panels (if condition_patterns column is present in samples CSV)
  condition_patterns_str <- if ("condition_patterns" %in% names(samples_tbl))
    samples_tbl$condition_patterns[i] else NA
  cond_patterns <- parse_condition_patterns(condition_patterns_str)

  if (length(cond_patterns) > 0) {
    # Extract conditions from cell names present in dropout data
    all_cells_in_dropout <- unique(cn_binned$cell)
    cond_df   <- extract_cell_conditions(all_cells_in_dropout, cond_patterns)
    conditions <- sort(setdiff(unique(cond_df$condition), "all"))

    message("  Condition panels: ", paste(conditions, collapse=", "))

    for (cond in conditions) {
      safe_cond  <- gsub("[^A-Za-z0-9_-]", "_", cond)
      cond_cells <- cond_df$cell[cond_df$condition == cond]

      if (length(cond_cells) < 3) {
        message("  Skipping condition '", cond, "' — too few cells")
        next
      }

      # Subset dropout data to this condition
      cn_binned_cond    <- cn_binned        %>% filter(cell %in% cond_cells)
      cell_chr_drop_cond <- cell_chr_dropout %>% filter(cell %in% cond_cells)

      # Per-condition plate heatmap
      cond_plate_png <- file.path(fig_dir, paste0("dropout_plate_heatmap_cond_", safe_cond, ".png"))
      tryCatch(
        plot_dropout_plate(
          cell_pos_dropout_path = file.path(sample_dir, "cell_pos_dropout.rds"),
          cell_chr_dropout      = cell_chr_drop_cond,
          sample_id             = paste0(sample_id, " [", cond, "]"),
          fig_dir               = fig_dir
        ),
        error = function(e) warning("Plate heatmap failed for condition '", cond, "': ", e$message)
      )
      # Rename the generic output to condition-specific name if it was written
      generic_plate <- file.path(fig_dir, paste0("dropout_plate_heatmap_", gsub("[^A-Za-z0-9_]", "", sample_id), "_.png"))
      # Simpler: just pass a dedicated fig_dir temp and rename, or call render_plate_heatmap directly
      # Instead, re-generate using render_plate_heatmap on the subsetted cells
      cell_totals_cond <- cell_chr_drop_cond %>%
        group_by(cell) %>%
        summarise(total_dropout = sum(dropout_count, na.rm = TRUE), .groups = "drop")
      plate_match <- regexpr("Plate_[0-9]+", cell_totals_cond$cell)
      cell_totals_cond$plate_id <- ifelse(plate_match > 0,
                                           regmatches(cell_totals_cond$cell, plate_match), NA)
      well_match <- regexpr("[A-P][0-9]+", cell_totals_cond$cell)
      cell_totals_cond$well <- ifelse(well_match > 0,
                                       regmatches(cell_totals_cond$cell, well_match), NA)
      cell_totals_cond$AZ     <- sub("[0-9]+$", "", cell_totals_cond$well)
      cell_totals_cond$Number <- as.integer(sub("^[A-Z]+", "", cell_totals_cond$well))
      valid_cond <- cell_totals_cond %>% filter(grepl("^[A-P]$", AZ), !is.na(Number))

      if (nrow(valid_cond) > 0) {
        for (pid in sort(unique(valid_cond$plate_id))) {
          plate_data <- valid_cond %>%
            filter(plate_id == pid) %>%
            select(AZ, Number, total_dropout)
          render_plate_heatmap(
            plate      = plate_data,
            sample_id  = sample_id,
            plate_label = paste0("[", cond, "]"),
            n_cells    = nrow(plate_data),
            outfile    = file.path(fig_dir, paste0("dropout_plate_heatmap_cond_", safe_cond, ".png"))
          )
        }
        cond_plate_png <- file.path(fig_dir, paste0("dropout_plate_heatmap_cond_", safe_cond, ".png"))
      } else {
        cond_plate_png <- NULL
      }

      # Load per-condition qc_metric_data saved by step 01
      cond_metric_path <- file.path(sample_dir, paste0("qc_metric_data_", safe_cond, ".rds"))
      if (file.exists(cond_metric_path)) {
        qc_cond <- tryCatch(readRDS(cond_metric_path), error = function(e) NULL)
        if (!is.null(qc_cond)) {
          cn_png_cond <- if (!is.null(qc_cond$cn_heatmap_png) &&
                             file.exists(qc_cond$cn_heatmap_png))
                           qc_cond$cn_heatmap_png else NULL
          plate_png_cond <- if (!is.null(cond_plate_png) && file.exists(cond_plate_png))
                              cond_plate_png else NULL
          tryCatch(
            plot_qc_summary_panel(
              df                 = qc_cond$df,
              iq                 = qc_cond$iq,
              summary_df         = qc_cond$summary_df,
              sample_id          = sample_id,
              condition          = cond,
              outfile            = file.path(fig_dir, paste0("qc_summary_panel_", safe_cond, ".pdf")),
              cn_heatmap_path    = cn_png_cond,
              plate_heatmap_path = plate_png_cond
            ),
            error = function(e) warning("QC panel failed for condition '", cond, "': ", e$message)
          )
          message("  Updated QC summary panel for condition: ", cond)
        }
      } else {
        message("  No qc_metric_data_", safe_cond, ".rds — run step 01 with condition_patterns first.")
      }
    }
  }

  message("Done: ", sample_id)
}

message("\nAll per-sample dropout plots generated.")
