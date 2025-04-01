# Install and load necessary packages
install.packages("glmnet")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("Rtsne")
library(glmnet)
library(ggplot2)
library(reshape2)
library(Rtsne)

IHR1 <- readRDS("../../Volumes/Fl/sc_analysis/scAboslute-obj/SLX-24490_100.rds")
HRP1 <- readRDS("../../Volumes/Fl/sc_analysis/scAboslute-obj/SLX-24489_100.rds")

IHR2 <- readRDS("../../Volumes/Fl/sc_analysis/scAboslute-obj/SLX-24174_100.rds")
HRP2 <- readRDS("../../Volumes/Fl/sc_analysis/scAboslute-obj/SLX-24175_100.rds")


HRP_df1 <- as.data.frame(HRP1@assayData$copynumber)
IHR_df1 <- as.data.frame(IHR1@assayData$copynumber)
HRP_df2 <- as.data.frame(HRP2@assayData$copynumber)
IHR_df2 <- as.data.frame(IHR2@assayData$copynumber)
# Combine the matrices and create a response vector
X <- cbind(HRP_df1, IHR_df1,HRP_df2, IHR_df2)
y <- c(rep(0, ncol(HRP_df1)), rep(1, ncol(IHR_df1)),rep(0, ncol(HRP_df2)), rep(2, ncol(IHR_df2)))

# Remove rows with NA values
X <- X[complete.cases(X), ]

# Fit the lasso model
lasso_model <- cv.glmnet(t(as.matrix(X)), y, alpha = 1, family = "binomial")

# Get the coefficients of the model
lasso_coefficients <- coef(lasso_model, s = "lambda.min")

# Identify the important positions in the genome
important_positions <- which(lasso_coefficients != 0)


# Print the important positions
print(important_positions)
rownames(X)[important_positions]




# Subset the data to include only the important features
X_important <- X[ important_positions,]


# Perform t-SNE
tsne_result <- Rtsne(as.matrix(t(X)), perplexity = 7, check_duplicates = FALSE)

# Create a data frame for t-SNE results
tsne_data <- data.frame(tsne_result$Y)
tsne_data$label <- y
colnames(tsne_data) <- c("Dim1", "Dim2", "label")

# Create the t-SNE plot
ggplot(tsne_data, aes(x = Dim1, y = Dim2, color = as.factor(label))) +
  geom_point() +
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2", color = "Label") +
  theme_minimal()












# Load necessary libraries
library(stats)     # For PCA & clustering
library(ggplot2)   # For visualization
library(dendextend) # For dendrogram customization

# ---- Load Your Data ----
# Assuming `data` is a matrix/dataframe of shape (1470, 23607)
# You can replace this with actual file loading
# data <- read.csv("your_file.csv", row.names = 1)  # Example for CSV

# ---- Step 1: Normalize the Data ----
data_scaled <- scale(t(X))  # Standardize the data (mean=0, sd=1)

# ---- Step 2: Perform PCA ----
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Choose number of principal components (e.g., keeping 90% variance)
variance_explained <- summary(pca_result)$importance[2, ]
num_pcs <- which(cumsum(variance_explained) >= 0.90)[1]  # Keep 90% variance

# Transform data using selected PCs
pca_data <- pca_result$x[, 1:num_pcs]

# ---- Step 3: Compute Distance Matrix ----
dist_matrix <- dist(pca_data, method = "euclidean")  # Euclidean distance

# ---- Step 4: Apply Hierarchical Clustering ----
hc <- hclust(dist_matrix, method = "ward.D2")  # Ward's method minimizes variance

# ---- Step 5: Visualize the Dendrogram ----
plot(hc, labels = y, main = "Hierarchical Clustering Dendrogram", cex = 0.6)

# Optionally, cut the tree into clusters (e.g., k=5)
clusters <- cutree(hc, k = 5)

# ---- Step 6: Add Cluster Labels to Data ----
final_data <- data.frame(Sample_ID = rownames(t(X)), Cluster = clusters)

# ---- Step 7: Visualize Clusters (if PCA reduced to 2D) ----
if (num_pcs >= 2) {
  ggplot(data.frame(PC1 = pca_data[, 1], PC2 = pca_data[, 2], Cluster = factor(clusters)), 
         aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 2) +
    labs(title = "Clusters after PCA and Hierarchical Clustering") +
    theme_minimal()
}

# Save clusters
write.csv(final_data, "cluster_results.csv", row.names = FALSE)

