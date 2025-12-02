#!/usr/bin/env Rscript
## Runner: generate_summary_of_outliers_scAbsolute.R
## Usage:
# Rscript scripts/generate_summary_of_outliers_scAbsolute.R [samples_csv]
# If no CSV is provided, the script falls back to a built-in sample vector.

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL

library(Biobase)
library(dplyr)

# Source project-local helpers
source("R/core.R")
source("R/visualization_helpers.R")
source("R/summary_helpers.R")

# default locations (adjust if needed)
obj_base <- "/Volumes/Fl/sc_analysis/scAboslute-obj"
out_base <- "/Volumes/Fl/sc_analysis/post-scAbsolute"
bin_size <- "100"

if (!is.null(samples_csv) && file.exists(samples_csv)) {
  samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE)
  # expect columns: sample, category (optional)
  if (!"sample" %in% colnames(samples_tbl)) stop("samples CSV must contain a 'sample' column")
} else {
  # fallback: small example list (edit as needed)
  samples_tbl <- data.frame(sample = c("25393"), category = c("PEO1 parent"), stringsAsFactors = FALSE)
  message("No samples CSV provided or file not found — using fallback sample list. To run many samples, provide a CSV with 'sample' and optional 'category' columns.")
}

for (i in seq_len(nrow(samples_tbl))) {
  SAMPLENAME <- samples_tbl$sample[i]
  CATEGORY <- if ("category" %in% colnames(samples_tbl)) samples_tbl$category[i] else NA
  SAMPLEONJ <- paste0("SLX-", SAMPLENAME, "_", bin_size)
  ALLCELLS <- file.path(obj_base, paste0(SAMPLEONJ, '.rds'))
  OUTPUT <- file.path(out_base, SAMPLEONJ)
  if (!dir.exists(OUTPUT)) dir.create(OUTPUT, recursive = TRUE)

  OUTLIERS_CELLBASE_PATH <- file.path(OUTPUT, "cellbased_outliers.rds")
  NON_OUTLIER_OBJ_PATH <- file.path(OUTPUT, paste0(SAMPLENAME, "_non_outlier.rds"))
  OUTLIER_SUMMARY_PATH <- file.path(OUTPUT, "outlier_summary.rds")
  NORMALS_PATH <- file.path(OUTPUT, "normals.rds")

  if (!file.exists(ALLCELLS)) {
    warning(paste0("Object file not found: ", ALLCELLS, " — skipping sample ", SAMPLENAME))
    next
  }

  message("Processing sample: ", SAMPLENAME)
  object <- tryCatch({ readRDS(ALLCELLS) }, error = function(e){
    warning(paste0("Failed to read RDS for sample ", SAMPLENAME, ": ", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(object)) next

  # Call the refactored helper (returns list). Enable saving by passing paths.
  res <- tryCatch({
    # if no category/UID is provided in the samples CSV, fall back to the sample name
    uid_to_pass <- if (!is.na(CATEGORY) && nzchar(CATEGORY)) CATEGORY else SAMPLENAME
    get_summary_of_outliers(
      object = object,
      SLX = paste0("SLX-", SAMPLENAME),
      UID = uid_to_pass,
      rpc_cutoff = 25,
      mapd_cutoff = 2,
      mapd_density_control = FALSE,
      gini_norm_cutoff = 2,
      alpha_cutoff = 1.5,
      alpha_hard_cutoff = 0.05,
      gini_density_control = TRUE,
      densityCutoff = 0.1,
      cutoff_percentile = 0.05,
      save_results = TRUE,
      outlier_cellbase_path = OUTLIERS_CELLBASE_PATH,
      non_outlier_obj_path = NON_OUTLIER_OBJ_PATH,
      outlier_summary_path = OUTLIER_SUMMARY_PATH,
      normals_path = NORMALS_PATH
    )
  }, error = function(e) {
    warning("Error running get_summary_of_outliers for sample ", SAMPLENAME, ": ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(res)) next

  # Plot heatmap of non-outlier object if available
  if (!is.null(res$non_outlier_object)) {
    # guard: ensure there is valid copy-number data to plot
    cn_mat <- tryCatch({ res$non_outlier_object@assayData$copynumber }, error = function(e) NULL)
    if (!is.null(cn_mat) && ncol(cn_mat) > 0 && any(!is.na(cn_mat))) {
      tryCatch({
        plotCopynumberHeatmap(res$non_outlier_object,
                              file = file.path(OUTPUT, "heatmap_clustered.pdf"),
                              cluster_rows = "ploidy",
                              row_split = NULL,
                              cutoff = 10,
                              show_unobserved_states = TRUE,
                              har = NULL,
                              useCopynumber = TRUE,
                              show_cell_names = FALSE,
                              abbreviate_cell_names = TRUE,
                              show_chromosome_names = FALSE,
                              fontsize_row = 9,
                              fontsize_col = 9,
                              fontsize_chr = 12,
                              fontsize_leg_title = 18,
                              fontsize_leg_label = 14,
                              row_title = SAMPLENAME)
      }, error = function(e) {
        warning("Failed to plot heatmap for sample ", SAMPLENAME, ": ", conditionMessage(e))
      })
    } else {
      message("Skipping heatmap: no valid copy-number data for sample ", SAMPLENAME)
    }
  }
  message("Done: ", SAMPLENAME)
}


