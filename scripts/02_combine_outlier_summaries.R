#!/usr/bin/env Rscript

# Improved visualization with metadata support
# - Accepts optional CSV manifest as first argument (columns: sample, category, HRD_HRP, Resistance)
# - If no manifest provided, falls back to the hard-coded vectors below
# - Configurable `obj_base` and `out_base` via args 2/3

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- if (length(args) >= 1) args[[1]] else NULL
out_base <- if (length(args) >= 2) args[[2]] else "/Volumes/Fl/sc_analysis/all_samples_23july2024"
obj_base <- if (length(args) >= 3) args[[3]] else "/Volumes/Fl/sc_analysis/post-scAbsolute"
bin_size <- if (length(args) >= 4) args[[4]] else "100"

library(dplyr)
library(knitr)

# Hard-coded metadata vectors (kept for backward compatibility)
ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24491","24518","24835")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116 BRCA2 -/-","PEO14","UWB1.289 BRCA","PEO1","PEO1","UWB1.289","PEO1 MISSENSE","PEO1 STOP","PEO1")
HRD_HRP <- c("HRD","HRD","HRD","HRP","HRP","HRP","HRP","HRP","HRD","HRP","HRP","HRD","HRD","HRD","HRD","HRP","HRD","HRD","HRD","HRD","HRD","HRD")
Resistance <- c("Cisplatin sensitive","Cisplatin sensitive","Cisplatin sensitive","Cisplatin Resistant","Cisplatin Resistant","Cisplatin Resistant","-","Cisplatin Resistant","-","Cisplatin Resistant","BRCA2 WT","-","Cisplatin Sensitive","BRCA2 null","Cisplatin Sensitive","BRCA1 restored","Cisplatin Sensitive","Cisplatin Sensitive","BCRA1 mutation","-","-","Cisplatin Sensitive")

# If a samples CSV is provided, use it. Accept either the old columns
# (`category`, `HRD_HRP`, `Resistance`) or the new ones
# (`Cell line`/`Cell_line`, `feature1`, `feature2`). We normalize to
# internal names `category`, `feature1`, `feature2`.
if (!is.null(samples_csv) && file.exists(samples_csv)) {
  samples_tbl <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% colnames(samples_tbl)) stop("samples CSV must contain a 'sample' column")

  # prefer new-style names if present
  if (all(c("Cell line", "feature1", "feature2") %in% colnames(samples_tbl))) {
    samples_meta <- samples_tbl %>% dplyr::select(sample, `Cell line`, feature1, feature2)
    colnames(samples_meta) <- c("sample", "category", "feature1", "feature2")
  } else if (all(c("Cell_line", "feature1", "feature2") %in% colnames(samples_tbl))) {
    samples_meta <- samples_tbl %>% dplyr::select(sample, Cell_line, feature1, feature2)
    colnames(samples_meta) <- c("sample", "category", "feature1", "feature2")
  } else if (all(c("category","HRD_HRP","Resistance") %in% colnames(samples_tbl))) {
    # backward compatible: map legacy names
    samples_meta <- samples_tbl %>% dplyr::select(sample, category = category, feature1 = HRD_HRP, feature2 = Resistance)
  } else if ("category" %in% colnames(samples_tbl)) {
    samples_meta <- samples_tbl %>% dplyr::mutate(HRD_HRP=NA, Resistance=NA) %>% dplyr::select(sample, category, HRD_HRP, Resistance)
    colnames(samples_meta) <- c("sample", "category", "feature1", "feature2")
  } else {
    samples_meta <- NULL
  }
} else {
  samples_meta <- NULL
}

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Build working list from either CSV or hard-coded vectors
# Build working list from either CSV or hard-coded vectors
if (!is.null(samples_meta)) {
  ALLSAMPLES_USE <- samples_meta$sample
  categories_use <- samples_meta$category
  feature1_use <- samples_meta$feature1
  feature2_use <- samples_meta$feature2
} else {
  ALLSAMPLES_USE <- ALLSAMPLES
  categories_use <- categories
  feature1_use <- HRD_HRP
  feature2_use <- Resistance
}

# Validate lengths
if (length(ALLSAMPLES_USE) != length(categories_use) || length(ALLSAMPLES_USE) != length(feature1_use) || length(ALLSAMPLES_USE) != length(feature2_use)) {
  stop("Metadata vectors (sample, category, feature1, feature2) must have the same length")
}

# Initialize
combined_outlier_summary <- tibble()

for (i in seq_along(ALLSAMPLES_USE)) {
  SAMPLENAME <- ALLSAMPLES_USE[i]
  SAMPLEONJ <- paste0("SLX-", SAMPLENAME, "_", bin_size)
  OUTPUT <- file.path(obj_base, SAMPLEONJ)
  out_rds <- file.path(OUTPUT, "outlier_summary.rds")

  if (!file.exists(out_rds)) {
    warning("Missing outlier summary: ", out_rds, " — skipping sample ", SAMPLENAME)
    next
  }

  outlier_summary <- tryCatch(readRDS(out_rds), error = function(e){ warning("Failed to read: ", out_rds); return(NULL) })
  if (is.null(outlier_summary)) next

  outlier_summary <- as.data.frame(outlier_summary)
  outlier_summary$Sample <- as.character(SAMPLENAME)
  # write the new column names: `Cell line`, `feature1`, `feature2`
  outlier_summary[["Cell line"]] <- as.character(categories_use[i])
  outlier_summary$feature1 <- as.character(feature1_use[i])
  outlier_summary$feature2 <- as.character(feature2_use[i])

  combined_outlier_summary <- bind_rows(combined_outlier_summary, outlier_summary)
}

if (nrow(combined_outlier_summary) == 0) stop("No outlier summaries found — nothing to visualize.")

# Define order and arrange
category_order <- c("PEO1", "PEO4", "PEO6", "PEO1 MISSENSE", "PEO1 STOP", "PEO14", "PEO23", "CIOV3", "CIOV6", "HCT116", "HCT116 BRCA2 -/-", "UWB1.289", "UWB1.289 BRCA")
if ("Category" %in% colnames(combined_outlier_summary)) {
  combined_outlier_summary$Category <- factor(combined_outlier_summary$Category, levels = category_order)
}

sorted_data <- combined_outlier_summary %>% arrange(`Cell line`) %>% select(Sample, `Cell line`, feature1, feature2, everything())

# Ensure Sample/Category are character and add totals
sorted_data$Sample <- as.character(sorted_data$Sample)
sorted_data$`Cell line` <- as.character(sorted_data$`Cell line`)
sum_row <- sorted_data %>% summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% mutate(Sample = "Total", `Cell line` = "", feature1 = "", feature2 = "")
sorted_data <- bind_rows(sorted_data, sum_row)

# Always write CSV; optional rendering below
outfile_csv <- file.path(out_base, "outlier_summary_table_meta_combined.csv")
write.csv(sorted_data, outfile_csv, row.names = FALSE)
message("Wrote: ", outfile_csv)

# Try to render PNG if packages available
if (requireNamespace("kableExtra", quietly = TRUE) && requireNamespace("webshot2", quietly = TRUE)) {
  outfile_png <- file.path(out_base, "outlier_summary_table_meta.png")
  kable(sorted_data, format = "html") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    kableExtra::add_header_above(c(" " = 4, "Datasets" = ncol(sorted_data) - 4)) %>%
    kableExtra::row_spec(nrow(sorted_data), bold = TRUE, background = "#EFEFEF") %>%
    kableExtra::column_spec(4, width = "150px", extra_css = "white-space: nowrap;") %>%
    kableExtra::save_kable(file = outfile_png, zoom = 2)
  message("Wrote: ", outfile_png)
} else {
  message("Optional rendering packages not available — CSV was written to: ", outfile_csv)
}
