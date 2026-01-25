#!/usr/bin/env Rscript

# Cluster Classification Analysis
# --------------------------------
# Performs machine learning classification and clustering on scAbsolute copy number data:
# 1. LASSO regression for feature selection (identifying important genomic positions)
# 2. t-SNE dimensionality reduction for visualization
# 3. PCA + hierarchical clustering for sample grouping
#
# Usage:
#   Rscript scripts/cluster_classification_analysis.R <input_dir> <output_dir> <bin_size> <sample1,sample2,...> <labels>
#
# Example:
#   Rscript scripts/cluster_classification_analysis.R /path/to/scAbsolute-obj /path/to/output 100 "24490,24489,24174,24175" "1,0,1,0"
#
# Arguments:
#   input_dir  - Directory containing scAbsolute RDS objects
#   output_dir - Directory to save results
#   bin_size   - Bin size used (e.g., 100)
#   samples    - Comma-separated sample IDs
#   labels     - Comma-separated labels for each sample (for classification)

args <- commandArgs(trailingOnly = TRUE)

# Default values for interactive use
if (length(args) < 5) {
  message("Usage: Rscript cluster_classification_analysis.R <input_dir> <output_dir> <bin_size> <samples> <labels>")
  message("Running with default example values...")
  input_dir <- "/Volumes/Fl/sc_analysis/scAbsolute-obj"
  output_dir <- "/Volumes/Fl/sc_analysis/cluster_analysis"
  bin_size <- "100"
  samples <- c("24490", "24489", "24174", "24175")
  labels <- c(1, 0, 1, 0)  # IHR=1, HRP=0
} else {
  input_dir <- args[1]
  output_dir <- args[2]
  bin_size <- args[3]
  samples <- unlist(strsplit(args[4], ","))
  labels <- as.numeric(unlist(strsplit(args[5], ",")))
}

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load required packages
library(glmnet)
library(ggplot2)
library(reshape2)
library(Rtsne)
library(stats)
library(dendextend)

set.seed(2020)

# ============================================================================
# Load Data
# ============================================================================
message("Loading scAbsolute objects...")

data_list <- list()
label_vector <- c()

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  rds_path <- file.path(input_dir, paste0("SLX-", sample_id, "_", bin_size, ".rds"))

  if (!file.exists(rds_path)) {
    warning("File not found: ", rds_path, " - skipping")
    next
  }

  obj <- readRDS(rds_path)
  df <- as.data.frame(obj@assayData$copynumber)
  data_list[[sample_id]] <- df
  label_vector <- c(label_vector, rep(labels[i], ncol(df)))
  message("  Loaded: ", sample_id, " (", ncol(df), " cells)")
}

if (length(data_list) == 0) {
  stop("No data files could be loaded.")
}

# Combine all samples
X <- do.call(cbind, data_list)
y <- label_vector

# Remove rows with NA values
X <- X[complete.cases(X), ]
message("Combined matrix: ", nrow(X), " bins x ", ncol(X), " cells")

# ============================================================================
# 1. LASSO Feature Selection
# ============================================================================
message("Running LASSO feature selection...")

lasso_model <- cv.glmnet(t(as.matrix(X)), y, alpha = 1, family = "binomial")
lasso_coefficients <- coef(lasso_model, s = "lambda.min")

# Identify important genomic positions
important_positions <- which(lasso_coefficients != 0)
important_features <- rownames(X)[important_positions]

message("Found ", length(important_positions), " important genomic positions")

# Save important positions
write.csv(
  data.frame(position = important_features),
  file.path(output_dir, "lasso_important_positions.csv"),
  row.names = FALSE
)

# ============================================================================
# 2. t-SNE Visualization
# ============================================================================
message("Running t-SNE...")

perplexity_value <- min(30, floor((ncol(X) - 1) / 3))
tsne_result <- Rtsne(as.matrix(t(X)), perplexity = perplexity_value, check_duplicates = FALSE)

tsne_data <- data.frame(
  Dim1 = tsne_result$Y[, 1],
  Dim2 = tsne_result$Y[, 2],
  Label = as.factor(y)
)

p_tsne <- ggplot(tsne_data, aes(x = Dim1, y = Dim2, color = Label)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(
    title = "t-SNE of Copy Number Profiles",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    color = "Group"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "tsne_plot.pdf"), p_tsne, width = 8, height = 6)
ggsave(file.path(output_dir, "tsne_plot.png"), p_tsne, width = 8, height = 6, dpi = 300)

# ============================================================================
# 3. PCA + Hierarchical Clustering
# ============================================================================
message("Running PCA and hierarchical clustering...")

# Normalize data
data_scaled <- scale(t(X))

# PCA
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Keep components explaining 90% variance
variance_explained <- summary(pca_result)$importance[2, ]
num_pcs <- which(cumsum(variance_explained) >= 0.90)[1]
message("Using ", num_pcs, " PCs (90% variance explained)")

pca_data <- pca_result$x[, 1:num_pcs]

# Hierarchical clustering
dist_matrix <- dist(pca_data, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

# Save dendrogram
pdf(file.path(output_dir, "dendrogram.pdf"), width = 12, height = 6)
plot(hc, labels = y, main = "Hierarchical Clustering Dendrogram", cex = 0.6)
dev.off()

# Cut tree into clusters
n_clusters <- length(unique(y))
clusters <- cutree(hc, k = n_clusters)

# Save cluster assignments
final_data <- data.frame(
  Cell_ID = colnames(X),
  True_Label = y,
  Cluster = clusters
)
write.csv(final_data, file.path(output_dir, "cluster_results.csv"), row.names = FALSE)

# PCA scatter plot with clusters
if (num_pcs >= 2) {
  pca_plot_data <- data.frame(
    PC1 = pca_data[, 1],
    PC2 = pca_data[, 2],
    Cluster = factor(clusters),
    True_Label = factor(y)
  )

  p_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Cluster, shape = True_Label)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = "PCA with Hierarchical Clustering",
      x = paste0("PC1 (", round(variance_explained[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(variance_explained[2] * 100, 1), "%)")
    ) +
    theme_minimal()

  ggsave(file.path(output_dir, "pca_clusters.pdf"), p_pca, width = 8, height = 6)
  ggsave(file.path(output_dir, "pca_clusters.png"), p_pca, width = 8, height = 6, dpi = 300)
}

message("Results saved to: ", output_dir)
message("Done!")
