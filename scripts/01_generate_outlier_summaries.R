#!/usr/bin/env Rscript
# ==============================================================================
# Generate Outlier Summaries for scAbsolute Objects
# ==============================================================================
#
# Usage:
#   Rscript scripts/01_generate_outlier_summaries.R <samples_csv> [obj_base] [out_base] [bin_size] [qc_config]
#
# Arguments:
#   samples_csv - CSV file with 'sample' column (and optional 'category')
#   obj_base    - Directory containing scAbsolute RDS objects
#   out_base    - Output directory for results
#   bin_size    - Bin size (default: 100)
#   qc_config   - QC parameters config file (default: config/qc_params_default.csv)
#
# Output per sample:
#   - outlier_summary.rds      : Summary statistics (RDS format)
#   - outlier_summary.csv      : Summary with QC criteria columns
#   - qc_criteria_used.csv     : Full QC parameters used
#   - cellbased_outliers.rds   : Per-cell outlier status
#   - {sample}_non_outlier.rds : Filtered object with good quality cells
#   - normals.rds              : Normal (diploid) cells
#   - heatmap_clustered.pdf    : Copy number heatmap
#
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
obj_base    <- if (length(args) >= 2) args[[2]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/scAboslute-obj"
out_base    <- if (length(args) >= 3) args[[3]] else "/Volumes/LenovoPS8/FI backup/sc_analysis/analysis_per_sample"
bin_size    <- if (length(args) >= 4) args[[4]] else "100"
qc_config   <- if (length(args) >= 5) args[[5]] else "config/qc_params_default.csv"

library(Biobase)
library(dplyr)
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

message("\n=== QC Parameters ===")
message("RPC cutoff: >= ", QC_PARAMS$rpc_cutoff)
message("MAPD cutoff: <= ", QC_PARAMS$mapd_cutoff, " (density_control: ", QC_PARAMS$mapd_density_control, ")")
message("Gini cutoff: <= ", QC_PARAMS$gini_norm_cutoff, " (density_control: ", QC_PARAMS$gini_density_control, ")")
message("Alpha cutoff: <= ", QC_PARAMS$alpha_cutoff, " SD (hard: ", QC_PARAMS$alpha_hard_cutoff, ")")
message("Replicating cutoff: ", QC_PARAMS$replicating_cutoff_value, " SD")
message("Normal threshold: > ", QC_PARAMS$normal_threshold, "% bins with CN=2")
message("=====================\n")

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
      gini_density_control = QC_PARAMS$gini_density_control,
      densityCutoff = QC_PARAMS$densityCutoff,
      cutoff_percentile = QC_PARAMS$cutoff_percentile,
      replicating_cutoff_value = QC_PARAMS$replicating_cutoff_value,
      replicating_iqr_value = QC_PARAMS$replicating_iqr_value,
      normal_threshold = QC_PARAMS$normal_threshold,
      # Output paths
      save_results = TRUE,
      outlier_cellbase_path = file.path(output_dir, "cellbased_outliers.rds"),
      non_outlier_obj_path  = file.path(output_dir, paste0(sample_id, "_non_outlier.rds")),
      outlier_summary_path  = file.path(output_dir, "outlier_summary.rds"),
      normals_path          = file.path(output_dir, "normals.rds")
    )
  }, error = function(e) { warning("QC failed for ", sample_id, ": ", e$message); NULL })

  if (is.null(res)) next

  # Save QC criteria in separate file (for reference)
  qc_criteria_df <- as.data.frame(res$qc_criteria)
  qc_criteria_df$config_file <- basename(qc_config)
  write.csv(qc_criteria_df, file.path(output_dir, "qc_criteria_used.csv"), row.names = FALSE)

  # Save summary WITHOUT QC columns (QC params go in header of combined output)
  write.csv(res$summary_df, file.path(output_dir, "outlier_summary.csv"), row.names = FALSE)

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

message("\nAll samples processed.")
message("QC config used: ", qc_config)
