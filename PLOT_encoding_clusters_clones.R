library(ggpubr)
library(grid)
library(dendextend)
library(pheatmap)

samples <- c("SLX-25394_B4-15-qc_100", "SLX-25394_C9-15-qc_100", "SLX-25394_D5-15-qc_100", "SLX-25394_D7-15-qc_100", "SLX-25394_D11-15-qc_100", "SLX-25394_G3-15-qc_100", "SLX-25393-15-qc_100")
base_path <- "../../Volumes/Fl/scUnique_results/"

get_uniqueEventEncoding <- function(sample) {
  object <- readRDS(paste0(base_path, sample, "/", sample, ".finalCN.RDS"))
  background <- readRDS(paste0(base_path, sample, "/", sample, ".profiles_background_freq.RDS"))
  
  # Compute unique event encoding
  uniqueEventEncoding <- matrix(-1, nrow=dim(object)[1], ncol=dim(object)[2])
  uniqueEventIndices <- object@assayData$copynumber != background & !is.na(object@assayData$copynumber)
  uniqueEventEncoding[uniqueEventIndices] <- object@assayData$copynumber[uniqueEventIndices]
  
  return(uniqueEventEncoding)
}


# Get unique event encodings for all samples
event_encodings <- lapply(samples, get_uniqueEventEncoding)

# Merge matrices by column
merged_matrix <- do.call(cbind, event_encodings)

sample_names <- c("B4", "C9", "D5", "D7", "D11", "G3", "Parent")

# Assign the sample names as column names (repeat for all cells in each sample)
colnames(merged_matrix) <- rep(sample_names, sapply(event_encodings, ncol))




# Perform clustering
distance_matrix <- dist(t(merged_matrix))  # Distance between sub-samples
clustering <- hclust(distance_matrix)

# Plot dendrogram
plot(clustering, main = "Clustering of Samples by Unique Event Encoding", xlab = "Samples", sub = "")

# Heatmap to visualize clustering with sample names
pheatmap(merged_matrix, clustering_distance_cols="euclidean", clustering_method="ward.D2",
         main = "Heatmap of Unique Event Encodings", labels_col = colnames(merged_matrix))




# Perform K-means clustering
set.seed(123)  # Ensure reproducibility
num_clusters <- length(samples)  # Choose k equal to the number of samples

non_constant_cols <- apply(t(merged_matrix), 2, var) != 0
filtered_matrix <- t(merged_matrix)[, non_constant_cols, drop = FALSE]


kmeans_result <- kmeans(t(filtered_matrix), centers = num_clusters)


# Perform PCA for visualization
pca_result <- prcomp(t(filtered_matrix), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Cluster = as.factor(kmeans_result$cluster),
  Sample = colnames(filtered_matrix)
)

# Scatter plot of PCA with colors based on clusters
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 0.5) +
  ggtitle("K-means Clustering of Samples (PCA Projection)") +
  theme_minimal()