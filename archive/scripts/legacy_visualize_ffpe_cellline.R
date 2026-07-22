#!/usr/bin/env Rscript

# Robust visualizer for summary CSVs produced by the generator.
# Usage: Rscript scripts/visualize_summary.R [combined_csv] [ffpe_csv] [out_base]

args <- commandArgs(trailingOnly = TRUE)
combined_csv <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else "/Volumes/Fl/sc_analysis/all_samples_23july2024/outlier_summary_table_meta_combined.csv"
ffpe_csv <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "/Volumes/Fl/sc_analysis/all_samples_23july2024/FFPE_95.csv"
out_base <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else "/Volumes/Fl/sc_analysis/all_samples_23july2024"

library(ggplot2)
library(tidyr)
library(dplyr)

# helper: return first matching column name or NULL
pick_col <- function(df, candidates) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  NULL
}

if (!file.exists(combined_csv)) stop("Combined CSV not found: ", combined_csv)
df <- read.csv(combined_csv, stringsAsFactors = FALSE, check.names = FALSE)

# Detect metadata / cell-line column
cell_line_col <- pick_col(df, c("Cell line", "Cell_line", "category", "Category", "sample", "Sample"))
if (!is.null(cell_line_col)) df$Cell_line <- as.character(df[[cell_line_col]]) else df$Cell_line <- NA_character_

# Detect numeric columns with fallbacks
rep_col <- pick_col(df, c("Replicating", "Replicating...RPC", "Replicating..RPC", "replicating"))
good_col <- pick_col(df, c("Good.Quality.Cells", "Good.Quality"))
normal_col <- pick_col(df, c("Normal", "Normal.Cells"))
seq_col <- pick_col(df, c("Sequenced", "Sequenced.Cells", "Sequenced"))

# If Sequenced missing, use a sensible default vector padded/truncated to nrows
default_sequenced <- c(288,288,288,134,384,384,384,384,384,384,384,384,384,384,384,288,384,384,384,384,288,384)
n <- nrow(df)
if (is.null(seq_col)) {
  seq_vec <- rep(NA_real_, n)
  m <- min(length(default_sequenced), n)
  seq_vec[seq_len(m)] <- default_sequenced[seq_len(m)]
  if (n > length(default_sequenced)) seq_vec[(length(default_sequenced)+1):n] <- default_sequenced[length(default_sequenced)]
  df$Sequenced <- seq_vec
} else {
  df$Sequenced <- suppressWarnings(as.numeric(df[[seq_col]]))
}

df$Replicating <- if (!is.null(rep_col)) suppressWarnings(as.numeric(df[[rep_col]])) else NA_real_
df$Good.Quality.Cells <- if (!is.null(good_col)) suppressWarnings(as.numeric(df[[good_col]])) else NA_real_
df$Normal.Cells <- if (!is.null(normal_col)) suppressWarnings(as.numeric(df[[normal_col]])) else NA_real_

# Compute high-quality percentage safely
df <- df %>% mutate(High_Quality_Percentage = round((coalesce(Good.Quality.Cells, 0) + coalesce(Replicating, 0)) * 100 / pmax(Sequenced, 1)))

# Ensure out_base exists
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# Save normalized CSV
write.csv(df, file.path(out_base, "Cellines_95_2.csv"), row.names = FALSE)

# Density plot (main cohort)
p_main <- ggplot(df, aes(x = High_Quality_Percentage)) +
  geom_density(fill = "skyblue", color = "navy", alpha = 0.7, linewidth = 0.6) +
  theme_minimal() + labs(title = "Distribution of High-Quality Cell Percentage", x = "High-Quality (%)", y = "Density")
ggsave(file.path(out_base, "Good_Quality_Distribution.pdf"), p_main, width = 6, height = 4)
message("Wrote: ", file.path(out_base, "Good_Quality_Distribution.pdf"))

# Helper: safe pivot_longer that returns empty tibble if no cols present
safe_pivot <- function(dat, rows, cols_vec) {
  dat_sub <- dat[rows, , drop = FALSE]
  if (!any(cols_vec %in% colnames(dat_sub))) return(tibble::tibble())
  pivot_longer(dat_sub, cols = any_of(cols_vec), names_to = "Kind", values_to = "Count")
}

# Stacked composition (top N cell lines)
df$Filtered <- df$Good.Quality.Cells - df$Normal.Cells
df$Cell_line <- ifelse(is.na(df$Cell_line) & "Sample" %in% colnames(df), as.character(df$Sample), df$Cell_line)
N <- min(22, nrow(df))
df_long <- safe_pivot(df, seq_len(N), c("Replicating", "Normal.Cells", "Filtered"))
cols <- c("Replicating" = "#81C784", "Filtered" = "#D1C4E9", "Normal.Cells" = "#F0E68C")
if (nrow(df_long) > 0) {
  p_comp <- ggplot(df_long, aes(x = factor(Cell_line, levels = unique(Cell_line)), y = Count, fill = Kind)) +
    geom_bar(stat = "identity") + scale_fill_manual(values = cols) +
    geom_text(data = df[seq_len(N), ], aes(x = factor(Cell_line, levels = unique(Cell_line)), y = Sequenced + 5, label = paste0(High_Quality_Percentage, "%")), inherit.aes = FALSE, size = 3, angle = 90, vjust = 0.5) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) + labs(title = "Cell Composition by Cell line", x = "Cell line", y = "Count")
  ggsave(file.path(out_base, "cell_composition_plot.pdf"), p_comp, width = 8, height = 6)
  message("Wrote: ", file.path(out_base, "cell_composition_plot.pdf"))
} else {
  message("Skipping cell composition plot: no composition columns present in combined CSV")
}

# FFPE optional processing
if (file.exists(ffpe_csv)) {
  df_ffpe <- read.csv(ffpe_csv, stringsAsFactors = FALSE, check.names = FALSE)
  rep_col_f <- pick_col(df_ffpe, c("Replicating", "Replicating...RPC", "Replicating..RPC"))
  good_col_f <- pick_col(df_ffpe, c("Good.Quality.Cells", "Good.Quality"))
  normal_col_f <- pick_col(df_ffpe, c("Normal", "Normal.Cells"))
  seq_col_f <- pick_col(df_ffpe, c("Sequenced", "Sequenced.Cells", "Sequenced"))

  df_ffpe$Sequenced <- if (!is.null(seq_col_f)) suppressWarnings(as.numeric(df_ffpe[[seq_col_f]])) else NA_real_
  df_ffpe$Replicating <- if (!is.null(rep_col_f)) suppressWarnings(as.numeric(df_ffpe[[rep_col_f]])) else NA_real_
  df_ffpe$Good.Quality.Cells <- if (!is.null(good_col_f)) suppressWarnings(as.numeric(df_ffpe[[good_col_f]])) else NA_real_
  df_ffpe$Normal.Cells <- if (!is.null(normal_col_f)) suppressWarnings(as.numeric(df_ffpe[[normal_col_f]])) else NA_real_

  df_ffpe <- df_ffpe %>% mutate(High_Quality_Percentage = round((coalesce(Good.Quality.Cells, 0) + coalesce(Replicating, 0)) * 100 / pmax(Sequenced, 1)))
  write.csv(df_ffpe, file.path(out_base, "FFPE_95_2.csv"), row.names = FALSE)

  df_ffpe$Filtered <- df_ffpe$Good.Quality.Cells - df_ffpe$Normal.Cells
  if (!"Sample" %in% colnames(df_ffpe)) df_ffpe$Sample <- seq_len(nrow(df_ffpe))
  M <- min(19, nrow(df_ffpe))
  df_long_ffpe <- safe_pivot(df_ffpe, seq_len(M), c("Replicating", "Normal.Cells", "Filtered"))

  p_ffpe_dist <- ggplot(df_ffpe, aes(x = High_Quality_Percentage)) + geom_density(fill = "skyblue", color = "navy", alpha = 0.7, linewidth = 0.6) + theme_minimal() + labs(title = "High-Quality % (FFPE)")
  ggsave(file.path(out_base, "Good_Quality_Distribution_FFPE.pdf"), p_ffpe_dist, width = 6, height = 4)
  message("Wrote: ", file.path(out_base, "Good_Quality_Distribution_FFPE.pdf"))

  if (nrow(df_long_ffpe) > 0) {
    p_ffpe_comp <- ggplot(df_long_ffpe, aes(x = factor(Sample, levels = unique(as.character(Sample))), y = Count, fill = Kind)) + geom_bar(stat = "identity") + scale_fill_manual(values = cols) + theme_minimal() + labs(title = "Cell Composition in FFPE Samples", x = "Sample", y = "Count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    ggsave(file.path(out_base, "cell_composition_plot_FFPE.pdf"), p_ffpe_comp, width = 8, height = 6)
    message("Wrote: ", file.path(out_base, "cell_composition_plot_FFPE.pdf"))
  } else {
    message("Skipping FFPE composition plot: no composition columns present in FFPE CSV")
  }

  # Combined PDFs (guarded)
  library(gridExtra)
  tryCatch({
    pdf(file.path(out_base, "Good_Quality_Distribution_Combined.pdf"), width = 10, height = 4)
    grid.arrange(p_main, p_ffpe_dist, ncol = 2)
    dev.off()
    message("Wrote: ", file.path(out_base, "Cell_Composition_Combined.pdf"))
  }, error = function(e) message("Skipping combined composition PDF: ", e$message))
} else {
  message("FFPE CSV not found - skipping FFPE plots: ", ffpe_csv)
}
ggsave(file.path(out_base, "Good_Quality_Distribution.pdf"), p_main, width = 6, height = 4)
message("Wrote: ", file.path(out_base, "Good_Quality_Distribution.pdf"))

# Helper: safe pivot_longer that returns empty tibble if no cols present
safe_pivot <- function(dat, rows, cols_vec) {
  dat_sub <- dat[rows, , drop = FALSE]
  if (!any(cols_vec %in% colnames(dat_sub))) return(tibble::tibble())
  pivot_longer(dat_sub, cols = any_of(cols_vec), names_to = "Kind", values_to = "Count")
}

# Stacked composition (top N cell lines)
df$Filtered <- df$Good.Quality.Cells - df$Normal.Cells
df$Cell_line <- ifelse(is.na(df$Cell_line) & "Sample" %in% colnames(df), as.character(df$Sample), df$Cell_line)
N <- min(22, nrow(df))
