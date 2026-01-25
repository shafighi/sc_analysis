#!/usr/bin/env Rscript
# Usage: Rscript scripts/01_generate_outlier_summaries.R <samples_csv> [obj_base] [out_base] [bin_size]

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
obj_base    <- if (length(args) >= 2) args[[2]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/scAboslute-obj"
out_base    <- if (length(args) >= 3) args[[3]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute"
bin_size    <- if (length(args) >= 4) args[[4]] else "100"

library(Biobase)
library(dplyr)
source("R/core.R")
source("R/visualization_helpers.R")
source("R/summary_helpers.R")

# Load samples from CSV or use fallback
if (!is.null(samples_csv) && file.exists(samples_csv)) {
  samples_tbl <- read.csv(samples_csv)
  if (!"sample" %in% colnames(samples_tbl)) stop("CSV must contain 'sample' column")
} else {
  samples_tbl <- data.frame(sample = "25393", category = "PEO1 parent")
  message("No CSV provided - using fallback sample. Provide CSV with 'sample' and optional 'category' columns.")
}

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

  # Run QC pipeline
  res <- tryCatch({
    get_summary_of_outliers(
      object = object,
      SLX = paste0("SLX-", sample_id),
      UID = if (!is.na(category) && nzchar(category)) category else sample_id,
      rpc_cutoff = 25, mapd_cutoff = 2, mapd_density_control = FALSE,
      gini_norm_cutoff = 2, alpha_cutoff = 1.5, alpha_hard_cutoff = 0.05,
      gini_density_control = TRUE, densityCutoff = 0.1, cutoff_percentile = 0.05,
      save_results = TRUE,
      outlier_cellbase_path = file.path(output_dir, "cellbased_outliers.rds"),
      non_outlier_obj_path  = file.path(output_dir, paste0(sample_id, "_non_outlier.rds")),
      outlier_summary_path  = file.path(output_dir, "outlier_summary.rds"),
      normals_path          = file.path(output_dir, "normals.rds")
    )
  }, error = function(e) { warning("QC failed for ", sample_id, ": ", e$message); NULL })

  if (is.null(res)) next

  # Generate heatmap if valid data exists
  if (!is.null(res$non_outlier_object)) {
    cn_mat <- tryCatch(res$non_outlier_object@assayData$copynumber, error = function(e) NULL)
    if (!is.null(cn_mat) && ncol(cn_mat) > 0 && any(!is.na(cn_mat))) {
      tryCatch({
        plotCopynumberHeatmap(
          res$non_outlier_object,
          file = file.path(output_dir, "heatmap_clustered.pdf"),
          cluster_rows = "ploidy", row_split = NULL, cutoff = 10,
          show_unobserved_states = TRUE, har = NULL, useCopynumber = TRUE,
          show_cell_names = FALSE, abbreviate_cell_names = TRUE,
          show_chromosome_names = FALSE,
          fontsize_row = 9, fontsize_col = 9, fontsize_chr = 12,
          fontsize_leg_title = 18, fontsize_leg_label = 14,
          row_title = sample_id
        )
      }, error = function(e) warning("Heatmap failed for ", sample_id, ": ", e$message))
    }
  }
  message("Done: ", sample_id)
}
