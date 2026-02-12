# ==============================================================================
# Dropout Analysis Helper Functions
# ==============================================================================
# Shared functions for analyzing copy number dropouts (bins with segVal == 0)
# across cells and chromosomes in single-cell copy number data.
#
# Used by: scripts/04_generate_dropout_summaries.R
#          scripts/05_visualize_dropout_heatmap.R
# ==============================================================================

library(dplyr)
library(tidyr)

#' Parse "chr:start-end" row names into a data frame
#'
#' @param df A data frame whose row names follow the format "chromosome:start-end"
#' @return A data frame with columns: chromosome, start, end
split_row_names <- function(df) {
  parts <- strsplit(rownames(df), "[:-]")
  data.frame(
    chromosome = sapply(parts, "[", 1),
    start      = as.integer(sapply(parts, "[", 2)),
    end        = as.integer(sapply(parts, "[", 3)),
    stringsAsFactors = FALSE
  )
}

#' Compute dropout counts per cell and per chromosome from a QDNAseq object
#'
#' @param object A QDNAseq copy number object (from scAbsolute)
#' @return A list with:
#'   - cn_binned: long-format data frame (chromosome, start, end, segVal, cell)
#'   - cell_chr_dropout: data frame with columns cell, chromosome, dropout_count
get_dropout_by_chromosome <- function(object) {
  cn <- object@assayData$copynumber
  coords <- split_row_names(cn)

  # Convert to long format
  cn_df <- as.data.frame(cn, check.names = FALSE)
  cn_df <- cbind(coords, cn_df)

  cn_binned <- cn_df %>%
    pivot_longer(
      cols      = -c(chromosome, start, end),
      names_to  = "cell",
      values_to = "segVal"
    ) %>%
    filter(!is.na(segVal))

  # Count dropouts (segVal == 0) per cell per chromosome
  cell_chr_dropout <- cn_binned %>%
    filter(segVal == 0) %>%
    count(cell, chromosome, name = "dropout_count")

  # Ensure all cell x chromosome combinations exist (fill missing with 0)
  all_cells <- unique(cn_binned$cell)
  all_chroms <- unique(cn_binned$chromosome)
  complete_grid <- expand.grid(
    cell       = all_cells,
    chromosome = all_chroms,
    stringsAsFactors = FALSE
  )
  cell_chr_dropout <- complete_grid %>%
    left_join(cell_chr_dropout, by = c("cell", "chromosome")) %>%
    mutate(dropout_count = replace_na(dropout_count, 0L))

  list(
    cn_binned        = cn_binned,
    cell_chr_dropout = cell_chr_dropout
  )
}

#' Assign QC status to cells based on cellbased_outliers flags
#'
#' @param cell_names Character vector of cell identifiers
#' @param cellbased_outliers Data frame with columns: cellid, replicating,
#'   rpc_outlier, alpha_outlier, gini_outlier, mapd_outlier, na_alpha_outliers
#' @return Character vector of status labels (same length as cell_names)
assign_cell_status <- function(cell_names, cellbased_outliers) {
  # Build lookup: first matching flag wins (priority order)
  flag_cols <- c(
    replicating      = "replicating",
    rpc_outlier      = "rpc_outlier",
    alpha_outlier    = "alpha_outlier",
    gini_outlier     = "gini_outlier",
    mapd_outlier     = "mapd_outlier",
    na_alpha_outlier = "na_alpha_outliers"
  )

  # Only use columns that actually exist in the data
  available <- flag_cols[flag_cols %in% names(cellbased_outliers)]

  status <- rep("pass", length(cell_names))
  for (i in seq_along(cell_names)) {
    cn <- cell_names[i]
    row <- cellbased_outliers[cellbased_outliers$cellid == cn, , drop = FALSE]
    if (nrow(row) == 0) next
    for (label in names(available)) {
      col <- available[[label]]
      if (!is.na(row[[col]][1]) && row[[col]][1] == 1) {
        status[i] <- label
        break
      }
    }
  }
  status
}
