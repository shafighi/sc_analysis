#!/usr/bin/env Rscript
# ==============================================================================
# Generate Outlier Summaries for scAbsolute Objects
# ==============================================================================
#
# Usage:
#   Rscript scripts/01_generate_outlier_summaries.R <samples_csv> <obj_base> <out_base> [bin_size] [qc_config]
#
# Arguments:
#   samples_csv - CSV file with 'sample' column (and optional 'category')
#   obj_base    - Directory containing scAbsolute RDS objects
#   out_base    - Output directory for results
#   bin_size    - Bin size (default: 100)
#   qc_config   - QC parameters config file (default: config/qc_params_default.csv)
#
# Output per sample (in out_base/SLX-{sample}_{bin}/):
#   - qc_summary.csv                        : Per-sample QC counts (replicating, outliers, PassedQC)
#   - qc_params.csv                         : QC thresholds used
#   - qc_cell_labels.rds                    : Per-cell QC flags (replicating, mapd/gini/alpha outlier, etc.)
#   - all_cells_qc.csv                      : Per-cell QC metrics, final status, and explicit failure reason(s)
#   - cells_passedqc.rds                    : QDNAseq object containing PassedQC cells only
#   - cells_borderline.rds                  : QDNAseq object containing cells queued for manual review
#   - borderline_cells.csv                  : Review table with QC metrics and blank decision/notes columns
#   - cells_normal.rds                      : Subset of PassedQC cells classified as diploid/normal
#   - figures/heatmap_clustered.pdf         : CN heatmap of all PassedQC cells
#   - figures/cn_profiles_passedqc.pdf      : Threshold cover + one annotated page per PassedQC cell
#   - figures/cn_profiles_outliers.pdf      : Threshold cover + one annotated page per outlier cell
#   - figures/cn_profiles_borderline.pdf    : Threshold cover + annotated profiles for manual review
#   - figures/normal_dist_alpha.pdf         : Alpha distribution with filtering cutoffs
#   - figures/cellcycle.jpg                 : Cell cycle activity (cycling_activity vs replicating)
#   - figures/cell_qc.jpg                   : Alpha coloured by outlier status
#   - figures/dmapd_cutoff.jpg              : MAPD coloured by outlier status
#   - figures/dgini_cutoff.jpg              : Gini coloured by outlier status
#   - figures/rpc_usedreads_ploidy_replicatings.jpg : RPC vs reads coloured by ploidy
#   - figures/rpc_usedreads_outliers_replicatings.jpg : RPC vs reads coloured by outlier
#   - figures/alpha_rpc_ploidy_replicatings.jpg : Alpha vs RPC coloured by ploidy
#   - figures/alpha_rpc_replicating_all.jpg : Alpha vs RPC coloured by replicating
#   - figures/alpha_used_reads.jpg          : Alpha vs used reads coloured by replicating
#   - figures/total.reads_all.pdf           : Total reads distribution (all cells)
#   - figures/total.reads_non_replicating.pdf : Total reads distribution (non-replicating)
#   - figures/total.reads_non_outliers.pdf  : Total reads distribution (PassedQC only)
#   - figures/qc_summary_panel.pdf          : Compact publication-quality QC panel (all metrics, one page)
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript scripts/01_generate_outlier_summaries.R <samples_csv> <obj_base> <out_base> [bin_size] [qc_config]")
}
samples_csv <- args[[1]]
obj_base    <- args[[2]]
out_base    <- args[[3]]
bin_size    <- if (length(args) >= 4) args[[4]] else "100"
qc_config   <- if (length(args) >= 5) args[[5]] else "config/qc_params_default.csv"

# Some upstream plotting helpers initialize R's default device even when their
# return value is used programmatically. Send any implicit device to a null PDF;
# all intended outputs below open explicit file devices.
options(device = function(...) grDevices::pdf(file = NULL))

suppressPackageStartupMessages({
  library(Biobase)
  library(QDNAseq)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(ggplot2)
  library(ggbeeswarm)
  library(ggpubr)
  library(patchwork)
  library(robustbase)
  library(readxl)
  library(future)
})
source("R/core.R")
source("R/visualization_helpers.R")
source("R/summary_helpers.R")

# ==============================================================================
# Load QC Parameters from config file
# ==============================================================================
load_qc_params <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("QC config file not found: ", config_path,
         "\nUse config/qc_params_default.csv or create your own.")
  }

  config <- read.csv(config_path, stringsAsFactors = FALSE)

  # Convert to named list with proper types
  params <- list()
  for (i in seq_len(nrow(config))) {
    param_name <- config$parameter[i]
    param_value <- config$value[i]

    # Convert to appropriate type
    if (param_value %in% c("TRUE", "FALSE")) {
      params[[param_name]] <- as.logical(param_value)
    } else if (!is.na(suppressWarnings(as.numeric(param_value)))) {
      params[[param_name]] <- as.numeric(param_value)
    } else {
      params[[param_name]] <- param_value
    }
  }

  return(params)
}

# Load QC parameters
message("Loading QC parameters from: ", qc_config)
QC_PARAMS <- load_qc_params(qc_config)
if (is.null(QC_PARAMS$borderline_extreme_cutoff)) {
  QC_PARAMS$borderline_extreme_cutoff <- 3
}
if (is.null(QC_PARAMS$borderline_max_metric_flags)) {
  QC_PARAMS$borderline_max_metric_flags <- 2
}

message("\n=== QC Parameters ===")
message("RPC cutoff: >= ", QC_PARAMS$rpc_cutoff)
message("MAPD cutoff: <= ", QC_PARAMS$mapd_cutoff, " (density_control: ", QC_PARAMS$mapd_density_control, ")")
message("Gini cutoff: <= ", QC_PARAMS$gini_norm_cutoff, " (density_control: ", QC_PARAMS$gini_density_control, ")")
message("Alpha cutoff: <= ", QC_PARAMS$alpha_cutoff, " SD (hard: ", QC_PARAMS$alpha_hard_cutoff, ")")
message("Borderline extreme cutoff: > ", QC_PARAMS$borderline_extreme_cutoff, " standardized residual")
message("Borderline maximum metric flags: ", QC_PARAMS$borderline_max_metric_flags)
message("Replicating cutoff: ", QC_PARAMS$replicating_cutoff_value, " SD")
message("Normal threshold: > ", QC_PARAMS$normal_threshold, "% bins with CN=2")
message("=====================\n")

# Plot the continuous bin signal as raw read counts / RPC, with the scAbsolute
# segmented absolute copy number overlaid in orange. In scAbsolute output the
# QDNAseq `calls` assay contains raw bin counts, while `copynumber` contains
# discrete integer calls; plotting `copynumber` as points hides them under the
# segmentation and makes a clean cell look artificially flat.
strip_cell_suffix <- function(x) sub("_[0-9]+$", "", x)

extract_segment_runs <- function(segmented_cn, chromosome) {
  if (length(segmented_cn) != length(chromosome)) {
    stop("Segment values and chromosome labels must have equal lengths")
  }
  empty <- data.frame(
    value = numeric(0), start_idx = integer(0), end_idx = integer(0)
  )
  if (!length(segmented_cn)) return(empty)

  same_as_previous <- c(
    FALSE,
    chromosome[-1] == chromosome[-length(chromosome)] &
      is.finite(segmented_cn[-1]) &
      is.finite(segmented_cn[-length(segmented_cn)]) &
      abs(segmented_cn[-1] - segmented_cn[-length(segmented_cn)]) < 1e-8
  )
  segment_id <- cumsum(!same_as_previous)
  finite_segment <- is.finite(segmented_cn)
  segment_indices <- split(which(finite_segment), segment_id[finite_segment])
  if (!length(segment_indices)) return(empty)

  do.call(rbind, lapply(segment_indices, function(idx) data.frame(
    value = segmented_cn[idx[1]],
    start_idx = min(idx),
    end_idx = max(idx)
  )))
}

format_qc_number <- function(x, digits = 2, missing = "not evaluated") {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), missing)
}

draw_qc_cover_page <- function(document_title, n_cells, qc_criteria) {
  density_label <- function(x) if (isTRUE(x)) "density-based" else "standardized residual"
  lines <- c(
    paste0("Cells in this PDF: ", n_cells, " | Bin size: ", bin_size, " kb"),
    "Filtering order: replication -> RPC -> MAPD/Gini/alpha -> Borderline review -> PassedQC",
    "",
    paste0(
      "Replication: cycling activity above the batch cutoff (",
      qc_criteria$replicating_cutoff_value, " SD; IQR multiplier ",
      qc_criteria$replicating_iqr_value, " where available)"
    ),
    paste0("RPC: fail when RPC < ", qc_criteria$rpc_cutoff),
    paste0(
      "MAPD: flag cutoff ", qc_criteria$mapd_cutoff, " (",
      density_label(qc_criteria$mapd_density_control), ")"
    ),
    paste0(
      "Gini: flag cutoff ", qc_criteria$gini_norm_cutoff, " (",
      density_label(qc_criteria$gini_density_control), ")"
    ),
    paste0(
      "Alpha: distribution cutoff ", qc_criteria$alpha_cutoff,
      " SD; hard cutoff > ", qc_criteria$alpha_hard_cutoff
    ),
    paste0(
      "Borderline: <= ", qc_criteria$borderline_max_metric_flags,
      " metric flags, MAPD/Gini residuals <= ", qc_criteria$borderline_extreme_cutoff,
      " SD, and alpha <= ", qc_criteria$alpha_hard_cutoff
    ),
    paste0("Normal classification: > ", qc_criteria$normal_threshold, "% of bins have integer CN=2"),
    "",
    "Profile legend: gray/black points = raw bin counts / RPC; orange = segmented absolute CN; dashed line = CN 2",
    "Cell pages report the final QC status, failure reason(s), and the metrics available at that filtering stage."
  )

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::text(0.05, 0.90, document_title, adj = c(0, 0.5), cex = 1.35, font = 2)
  graphics::text(
    0.05, 0.82, paste(lines, collapse = "\n"),
    adj = c(0, 1), cex = 0.72, family = "mono"
  )
}

plot_cn_profiles_pdf <- function(object, outfile, title_lookup = NULL,
                                 document_title = "Copy-number profiles",
                                 qc_criteria = NULL) {
  if (is.null(object) || ncol(object) == 0) return(invisible(0L))

  grDevices::pdf(outfile, width = 11, height = 3.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (!is.null(qc_criteria)) {
    draw_qc_cover_page(document_title, ncol(object), qc_criteria)
  }

  plotted <- 0L
  for (ci in seq_len(ncol(object))) {
    cell_obj <- object[, ci]
    full_name <- as.character(Biobase::pData(cell_obj)$name[1])
    lookup_name <- strip_cell_suffix(full_name)
    plot_title <- full_name
    if (!is.null(title_lookup) && lookup_name %in% names(title_lookup)) {
      plot_title <- unname(title_lookup[[lookup_name]])
    }

    ok <- tryCatch({
      valid <- QDNAseq:::binsToUse(cell_obj)
      fd <- Biobase::fData(cell_obj)[valid, , drop = FALSE]
      chromosome <- as.character(fd$chromosome)
      chromosome_levels <- unique(chromosome)

      chromosome_lengths <- vapply(chromosome_levels, function(chr) {
        max(fd$end[chromosome == chr], na.rm = TRUE)
      }, numeric(1))
      chromosome_offsets <- c(0, head(cumsum(chromosome_lengths), -1))
      names(chromosome_offsets) <- chromosome_levels
      chromosome_ends <- cumsum(chromosome_lengths)
      chromosome_centres <- chromosome_offsets + chromosome_lengths / 2

      x_start <- fd$start + chromosome_offsets[chromosome]
      x_end <- fd$end + chromosome_offsets[chromosome]
      x_mid <- (x_start + x_end) / 2

      raw_counts <- as.numeric(Biobase::assayDataElement(cell_obj, "calls"))[valid]
      rpc <- as.numeric(Biobase::pData(cell_obj)$rpc[1])
      if (!is.finite(rpc) || rpc <= 0) stop("RPC must be positive to scale raw bin counts")
      raw_absolute_cn <- raw_counts / rpc
      segmented_cn <- as.numeric(Biobase::assayDataElement(cell_obj, "segmented"))[valid]

      # Build runs directly. CGHbase::.makeSegments() is a plotting helper,
      # not a pure segment extractor, and can unexpectedly open Rplots.pdf.
      segment_table <- extract_segment_runs(segmented_cn, chromosome)
      finite_signal <- c(raw_absolute_cn, segmented_cn)
      finite_signal <- finite_signal[is.finite(finite_signal)]
      upper <- if (length(finite_signal)) {
        min(12, max(6, as.numeric(stats::quantile(finite_signal, 0.995, na.rm = TRUE)) + 1))
      } else {
        6
      }

      graphics::par(mar = c(4.2, 4.4, 5.4, 1.0))
      graphics::plot(
        NA_real_, NA_real_,
        xlim = c(0, max(chromosome_ends)), ylim = c(0, upper),
        xaxs = "i", yaxs = "i", xaxt = "n",
        xlab = "Chromosome", ylab = "Raw bin counts / RPC",
        main = plot_title, cex.main = 0.72
      )
      graphics::axis(1, at = chromosome_centres, labels = chromosome_levels)
      if (length(chromosome_ends) > 1) {
        graphics::abline(v = head(chromosome_ends, -1), lty = 3, col = "grey55")
      }
      graphics::abline(h = 2, lty = 2, col = "grey55")
      for (si in seq_len(nrow(segment_table))) {
        graphics::segments(
          x_start[segment_table$start_idx[si]], segment_table$value[si],
          x_end[segment_table$end_idx[si]], segment_table$value[si],
          col = "#D95F02", lwd = 2.2
        )
      }
      graphics::points(
        x_mid, raw_absolute_cn,
        pch = 16, cex = 0.28,
        col = grDevices::adjustcolor("black", alpha.f = 0.32)
      )
      graphics::mtext(
        paste0(sum(valid), " x ", bin_size, " kb bins, ", nrow(segment_table), " segments"),
        side = 3, adj = 0, line = 0.15, cex = 0.72, col = "grey30"
      )
      TRUE
    }, error = function(e) {
      message("  Skipping cell ", full_name, ": ", e$message)
      FALSE
    })
    plotted <- plotted + as.integer(ok)
  }
  invisible(plotted)
}

# plotCopynumberHeatmap(file = NULL) still initializes/draws on the current
# graphics device. Give it a disposable null device so it cannot create or
# overwrite a repository-level Rplots.pdf while returning the heatmap object.
build_cn_heatmap <- function(object, row_title) {
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotCopynumberHeatmap(
    object,
    file = NULL,
    cluster_rows = "ploidy", row_split = NULL, cutoff = 10,
    show_unobserved_states = TRUE, har = NULL, useCopynumber = TRUE,
    show_cell_names = FALSE, abbreviate_cell_names = TRUE,
    show_chromosome_names = FALSE,
    fontsize_row = 9, fontsize_col = 9, fontsize_chr = 12,
    fontsize_leg_title = 18, fontsize_leg_label = 14,
    row_title = row_title
  )
}

build_cell_qc_annotations <- function(object, res, qc_criteria) {
  labels <- res$summary_df_outliers
  phenotype <- res$df %>%
    dplyr::select(dplyr::any_of(c(
      "name", "total.reads", "used.reads", "rpc", "ploidy",
      "cycling_activity", "cycling_cutoff_sd"
    )))
  metric_qc <- res$iq %>%
    dplyr::transmute(
      name = name,
      mapd_z = as.numeric(dmapd),
      gini_z = as.numeric(dgini),
      alpha = as.numeric(alpha.select)
    )

  annotations <- labels %>%
    dplyr::left_join(phenotype, by = "name") %>%
    dplyr::left_join(metric_qc, by = "name")

  full_names <- as.character(Biobase::pData(object)$name)
  short_names <- strip_cell_suffix(full_names)
  full_name_lookup <- setNames(full_names, short_names)
  annotations$full_name <- unname(full_name_lookup[annotations$name])

  cn <- Biobase::assayDataElement(object, "copynumber")
  percent_cn2 <- colMeans(cn == 2, na.rm = TRUE) * 100
  names(percent_cn2) <- short_names
  annotations$percent_cn2 <- unname(percent_cn2[annotations$name])

  annotations$metric_flags <- vapply(seq_len(nrow(annotations)), function(i) {
    flags <- c(
      MAPD = annotations$dmap_outlier[i] == 1,
      Gini = annotations$gini_outlier[i] == 1,
      Alpha = annotations$alpha_outlier[i] == 1
    )
    paste(names(flags)[flags], collapse = "+")
  }, character(1))

  annotations$failure_reason <- vapply(seq_len(nrow(annotations)), function(i) {
    reasons <- character(0)
    if (annotations$replicating[i] == 1) {
      reasons <- c(reasons, "Replicating/S-phase")
    }
    if (annotations$rpc_outlier_non_replicating[i] == 1 ||
        annotations$rpc_outlier_replicating[i] == 1) {
      reasons <- c(reasons, paste0("Low RPC (<", qc_criteria$rpc_cutoff, ")"))
    }
    if (annotations$na_alpha_outliers[i] == 1) {
      reasons <- c(reasons, "Missing alpha")
    }
    if (nzchar(annotations$metric_flags[i])) {
      prefix <- if (annotations$borderline[i] == 1) "Borderline metric" else "Metric failure"
      reasons <- c(reasons, paste0(prefix, ": ", annotations$metric_flags[i]))
    }
    if (!length(reasons)) "None - all evaluated QC checks passed" else paste(reasons, collapse = "; ")
  }, character(1))

  annotations$status_display <- ifelse(
    annotations$qc_status == "PassedQC",
    ifelse(
      annotations$percent_cn2 > qc_criteria$normal_threshold,
      "PassedQC / Normal",
      "PassedQC / CNV"
    ),
    annotations$qc_status
  )

  annotations$title <- paste0(
    annotations$full_name,
    "\nQC: ", annotations$status_display,
    " | Reason: ", annotations$failure_reason,
    "\nRPC=", format_qc_number(annotations$rpc, 1),
    " | MAPD z=", format_qc_number(annotations$mapd_z, 2),
    " | Gini z=", format_qc_number(annotations$gini_z, 2),
    " | alpha=", format_qc_number(annotations$alpha, 4),
    ifelse(
      annotations$replicating == 1,
      paste0(
        " | cycle=", format_qc_number(annotations$cycling_activity, 4),
        " (cutoff ", format_qc_number(annotations$cycling_cutoff_sd, 4), ")"
      ),
      ""
    )
  )

  annotations
}

# ==============================================================================
# Load samples
# ==============================================================================
if (!is.null(samples_csv) && file.exists(samples_csv)) {
  samples_tbl <- read.csv(samples_csv)
  if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")
} else {
  samples_tbl <- data.frame(sample = "25393", category = "PEO1 parent")
  message("No CSV provided - using fallback sample. Provide CSV with 'sample' and optional 'category' columns.")
}

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
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  object <- tryCatch(readRDS(input_file), error = function(e) { warning(e$message); NULL })
  if (is.null(object)) next

  # Run QC pipeline with all parameters from config
  res <- tryCatch({
    get_summary_of_outliers(
      object = object,
      SLX = paste0("SLX-", sample_id),
      UID = if (!is.na(category) && nzchar(category)) category else sample_id,
      # Pass all QC parameters from config
      rpc_cutoff = QC_PARAMS$rpc_cutoff,
      mapd_cutoff = QC_PARAMS$mapd_cutoff,
      mapd_density_control = QC_PARAMS$mapd_density_control,
      gini_norm_cutoff = QC_PARAMS$gini_norm_cutoff,
      alpha_cutoff = QC_PARAMS$alpha_cutoff,
      alpha_hard_cutoff = QC_PARAMS$alpha_hard_cutoff,
      borderline_extreme_cutoff = QC_PARAMS$borderline_extreme_cutoff,
      borderline_max_metric_flags = QC_PARAMS$borderline_max_metric_flags,
      gini_density_control = QC_PARAMS$gini_density_control,
      densityCutoff = QC_PARAMS$densityCutoff,
      cutoff_percentile = QC_PARAMS$cutoff_percentile,
      replicating_cutoff_value = QC_PARAMS$replicating_cutoff_value,
      replicating_iqr_value = QC_PARAMS$replicating_iqr_value,
      normal_threshold = QC_PARAMS$normal_threshold,
      # Output paths
      save_results = TRUE,
      outlier_cellbase_path = file.path(output_dir, "qc_cell_labels.rds"),
      non_outlier_obj_path  = file.path(output_dir, "cells_passedqc.rds"),
      outlier_summary_path  = NULL,
      normals_path          = file.path(output_dir, "cells_normal.rds")
    )
  }, error = function(e) { warning("QC failed for ", sample_id, ": ", e$message); NULL })

  if (is.null(res)) next

  # Save QC criteria in separate file (for reference)
  qc_criteria_df <- as.data.frame(res$qc_criteria)
  qc_criteria_df$config_file <- basename(qc_config)
  write.csv(qc_criteria_df, file.path(output_dir, "qc_params.csv"), row.names = FALSE)

  write.csv(res$summary_df, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)

  # === Generate figures ===
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  # Build one auditable annotation table used by every profile PDF. Metrics that
  # were not evaluated because a cell failed an earlier stage are reported as NA.
  qc_annotations <- build_cell_qc_annotations(object, res, res$qc_criteria)
  write.csv(qc_annotations, file.path(output_dir, "all_cells_qc.csv"), row.names = FALSE)
  qc_title_lookup <- setNames(qc_annotations$title, qc_annotations$name)

  # 1. CN heatmap (PassedQC cells) — render once, save as PDF + PNG
  cn_heatmap_png <- NULL
  if (!is.null(res$non_outlier_object)) {
    cn_mat <- tryCatch(res$non_outlier_object@assayData$copynumber, error = function(e) NULL)
    if (!is.null(cn_mat) && ncol(cn_mat) > 0 && any(!is.na(cn_mat))) {
      cn_ht <- tryCatch(
        build_cn_heatmap(res$non_outlier_object, row_title = sample_id),
        error = function(e) { warning("Heatmap failed for ", sample_id, ": ", e$message); NULL }
      )
      if (!is.null(cn_ht)) {
        tryCatch({
          pdf(file.path(figures_dir, "heatmap_clustered.pdf"), width = 11, height = 4)
          ComplexHeatmap::draw(cn_ht)
          dev.off()
        }, error = function(e) { try(dev.off(), silent = TRUE); warning("Heatmap PDF failed: ", e$message) })
        png_path <- file.path(figures_dir, "heatmap_clustered.png")
        tryCatch({
          # 4000 × 1275 px @ 300 dpi ≈ 13.3" × 4.25" (ratio 3.14:1 matches panel)
          png(png_path, width = 4000, height = 1275, res = 300)
          ComplexHeatmap::draw(cn_ht)
          dev.off()
          cn_heatmap_png <- png_path
        }, error = function(e) { try(dev.off(), silent = TRUE); warning("Heatmap PNG failed: ", e$message) })
      }
    }
  }

  # 2. Per-cell CN profiles — PassedQC cells
  tryCatch(
    plot_cn_profiles_pdf(
      res$non_outlier_object,
      file.path(figures_dir, "cn_profiles_passedqc.pdf"),
      title_lookup = qc_title_lookup,
      document_title = paste0("SLX-", sample_id, " - PassedQC copy-number profiles"),
      qc_criteria = res$qc_criteria
    ),
    error = function(e) warning("CN profiles (PassedQC) failed for ", sample_id, ": ", e$message)
  )

  # 3. Per-cell CN profiles — outlier cells
  passedqc_names <- Biobase::pData(res$non_outlier_object)$name
  outlier_mask   <- !Biobase::pData(object)$name %in% passedqc_names
  if (sum(outlier_mask) > 0) {
    outlier_object <- object[, outlier_mask]
    tryCatch(
      plot_cn_profiles_pdf(
        outlier_object,
        file.path(figures_dir, "cn_profiles_outliers.pdf"),
        title_lookup = qc_title_lookup,
        document_title = paste0("SLX-", sample_id, " - Outlier copy-number profiles"),
        qc_criteria = res$qc_criteria
      ),
      error = function(e) warning("CN profiles (outliers) failed for ", sample_id, ": ", e$message)
    )
  }

  # 4. Borderline review bundle. Borderline cells still fail strict QC; the
  # separate object/table/PDF make the manual decision explicit and auditable.
  borderline_iq <- res$iq[res$iq$borderline, , drop = FALSE]
  borderline_names <- borderline_iq$name
  borderline_mask <- strip_cell_suffix(Biobase::pData(object)$name) %in% borderline_names
  borderline_object <- object[, borderline_mask]
  saveRDS(borderline_object, file.path(output_dir, "cells_borderline.rds"))

  flag_matrix <- cbind(
    MAPD = as.logical(borderline_iq$dmapd.outlier),
    Gini = as.logical(borderline_iq$dgini.outlier),
    Alpha = as.logical(borderline_iq$alpha.outlier)
  )
  borderline_reason <- if (nrow(borderline_iq)) {
    apply(flag_matrix, 1, function(x) paste(names(x)[x], collapse = "+"))
  } else {
    character(0)
  }
  full_name_lookup <- setNames(
    Biobase::pData(object)$name,
    strip_cell_suffix(Biobase::pData(object)$name)
  )
  borderline_object_idx <- match(
    unname(full_name_lookup[borderline_iq$name]),
    Biobase::pData(object)$name
  )
  borderline_cn <- Biobase::assayDataElement(object, "copynumber")[, borderline_object_idx, drop = FALSE]
  borderline_segments <- vapply(borderline_object_idx, function(ci) {
    valid <- QDNAseq:::binsToUse(object)
    seg <- Biobase::assayDataElement(object, "segmented")[valid, ci]
    chr <- QDNAseq:::chromosomes(object)[valid]
    nrow(extract_segment_runs(seg, chr))
  }, integer(1))
  borderline_review <- data.frame(
    name = borderline_iq$name,
    full_name = unname(full_name_lookup[borderline_iq$name]),
    borderline_reason = borderline_reason,
    total_reads = borderline_iq$total.reads,
    used_reads = borderline_iq$used.reads,
    rpc = borderline_iq$rpc,
    ploidy = borderline_iq$ploidy,
    cycling_activity = borderline_iq$cycling_activity,
    cycling_cutoff = borderline_iq$cycling_cutoff_sd,
    mapd = borderline_iq$mapd,
    mapd_z = borderline_iq$dmapd,
    gini_normalized = borderline_iq$gini_normalized,
    gini_z = borderline_iq$dgini,
    alpha = borderline_iq$alpha.select,
    segments = borderline_segments,
    dropout_bins = colSums(borderline_cn == 0, na.rm = TRUE),
    percent_cn2 = colMeans(borderline_cn == 2, na.rm = TRUE) * 100,
    metric_flag_count = borderline_iq$metric_flag_count,
    manual_decision = rep("", nrow(borderline_iq)),
    review_notes = rep("", nrow(borderline_iq)),
    stringsAsFactors = FALSE
  )
  write.csv(borderline_review, file.path(output_dir, "borderline_cells.csv"), row.names = FALSE)

  tryCatch(
    plot_cn_profiles_pdf(
      borderline_object,
      file.path(figures_dir, "cn_profiles_borderline.pdf"),
      title_lookup = qc_title_lookup,
      document_title = paste0("SLX-", sample_id, " - Borderline copy-number profiles"),
      qc_criteria = res$qc_criteria
    ),
    error = function(e) warning("CN profiles (borderline) failed for ", sample_id, ": ", e$message)
  )

  # 5. QC metric plots (alpha distribution, gini, mapd, rpc, cell cycle) — individual files
  tryCatch({
    iq_plot <- res$iq
    if (is.data.frame(iq_plot$dmapd.outlier)) iq_plot$dmapd.outlier <- iq_plot$dmapd.outlier[, 1]
    if (is.data.frame(iq_plot$dgini.outlier)) iq_plot$dgini.outlier <- iq_plot$dgini.outlier[, 1]
    if (is.data.frame(iq_plot$alpha.outlier)) iq_plot$alpha.outlier <- iq_plot$alpha.outlier[, 1]
    general_analysis_plots(res$df, iq_plot, figures_dir)
  }, error = function(e) warning("QC metric plots failed for ", sample_id, ": ", e$message))

  # 6. Save QC metric data so step 6 can regenerate the panel with dropout heatmap + CN heatmap
  saveRDS(
    list(df = res$df, iq = res$iq, summary_df = res$summary_df,
         cn_heatmap_png = cn_heatmap_png),
    file.path(output_dir, "qc_metric_data.rds")
  )

  # 7. Compact publication-quality QC summary panel (CN heatmap included; dropout added in step 6)
  tryCatch(
    plot_qc_summary_panel(
      df              = res$df,
      iq              = res$iq,
      summary_df      = res$summary_df,
      sample_id       = sample_id,
      outfile         = file.path(figures_dir, "qc_summary_panel.pdf"),
      cn_heatmap_path = cn_heatmap_png
    ),
    error = function(e) warning("QC summary panel failed for ", sample_id, ": ", e$message)
  )

  # 8. Per-condition QC panels (optional — requires 'condition_patterns' column in samples CSV)
  condition_patterns_str <- if ("condition_patterns" %in% names(samples_tbl))
    samples_tbl$condition_patterns[i] else NA
  cond_patterns <- parse_condition_patterns(condition_patterns_str)

  if (length(cond_patterns) > 0) {
    orig_names <- Biobase::pData(object)$name
    cond_df    <- extract_cell_conditions(orig_names, cond_patterns)
    conditions <- sort(setdiff(unique(cond_df$condition), "all"))

    message("  Condition patterns [", paste(cond_patterns, collapse=";"), "] → ",
            length(conditions), " condition(s): ", paste(conditions, collapse=", "))

    for (cond in conditions) {
      safe_cond  <- gsub("[^A-Za-z0-9_-]", "_", cond)
      cond_cells <- cond_df$cell[cond_df$condition == cond]
      cond_mask  <- Biobase::pData(object)$name %in% cond_cells

      if (sum(cond_mask) < 3) {
        message("  Skipping condition '", cond, "' — only ", sum(cond_mask), " cells")
        next
      }
      message("  Condition '", cond, "': ", sum(cond_mask), " cells")

      cond_object <- object[, cond_mask]

      res_cond <- tryCatch(
        get_summary_of_outliers(
          object                 = cond_object,
          SLX                    = paste0("SLX-", sample_id),
          UID                    = if (!is.na(category) && nzchar(category)) category else sample_id,
          rpc_cutoff             = QC_PARAMS$rpc_cutoff,
          mapd_cutoff            = QC_PARAMS$mapd_cutoff,
          mapd_density_control   = QC_PARAMS$mapd_density_control,
          gini_norm_cutoff       = QC_PARAMS$gini_norm_cutoff,
          alpha_cutoff           = QC_PARAMS$alpha_cutoff,
          alpha_hard_cutoff      = QC_PARAMS$alpha_hard_cutoff,
          borderline_extreme_cutoff = QC_PARAMS$borderline_extreme_cutoff,
          borderline_max_metric_flags = QC_PARAMS$borderline_max_metric_flags,
          gini_density_control   = QC_PARAMS$gini_density_control,
          densityCutoff          = QC_PARAMS$densityCutoff,
          cutoff_percentile      = QC_PARAMS$cutoff_percentile,
          replicating_cutoff_value = QC_PARAMS$replicating_cutoff_value,
          replicating_iqr_value  = QC_PARAMS$replicating_iqr_value,
          normal_threshold       = QC_PARAMS$normal_threshold,
          save_results           = FALSE
        ),
        error = function(e) { warning("QC failed for condition '", cond, "': ", e$message); NULL }
      )
      if (is.null(res_cond)) next

      # Per-condition CN heatmap
      cn_heatmap_cond_png <- NULL
      cond_cn_mat <- tryCatch(
        res_cond$non_outlier_object@assayData$copynumber, error = function(e) NULL
      )
      if (!is.null(cond_cn_mat) && ncol(cond_cn_mat) > 0 && any(!is.na(cond_cn_mat))) {
        cn_ht_cond <- tryCatch(
          build_cn_heatmap(
            res_cond$non_outlier_object,
            row_title = paste0(sample_id, " [", cond, "]")
          ),
          error = function(e) { warning("Heatmap failed for condition '", cond, "': ", e$message); NULL }
        )
        if (!is.null(cn_ht_cond)) {
          png_cond_path <- file.path(figures_dir, paste0("heatmap_clustered_", safe_cond, ".png"))
          tryCatch({
            png(png_cond_path, width = 4000, height = 1275, res = 300)
            ComplexHeatmap::draw(cn_ht_cond)
            dev.off()
            cn_heatmap_cond_png <- png_cond_path
          }, error = function(e) { try(dev.off(), silent = TRUE) })
        }
      }

      # Save per-condition qc_metric_data for step 6 to use when regenerating full panel
      saveRDS(
        list(df = res_cond$df, iq = res_cond$iq, summary_df = res_cond$summary_df,
             cn_heatmap_png = cn_heatmap_cond_png),
        file.path(output_dir, paste0("qc_metric_data_", safe_cond, ".rds"))
      )

      # Per-condition QC summary panel
      tryCatch(
        plot_qc_summary_panel(
          df              = res_cond$df,
          iq              = res_cond$iq,
          summary_df      = res_cond$summary_df,
          sample_id       = sample_id,
          condition       = cond,
          outfile         = file.path(figures_dir, paste0("qc_summary_panel_", safe_cond, ".pdf")),
          cn_heatmap_path = cn_heatmap_cond_png
        ),
        error = function(e) warning("QC panel failed for condition '", cond, "': ", e$message)
      )
    }
  }

  message("Done: ", sample_id)
}

message("\nAll samples processed.")
message("QC config used: ", qc_config)
