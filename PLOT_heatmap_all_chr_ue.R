library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(viridis)

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO23","HCT116\nBRCA2 -/-")

all_samples_ue_chr <- list()
# Loop to create sublists and add them to the main list
for (i in 1:length(samples)) {
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[i],'/')
  chromosomes <- c(paste0("chr", 1:22))
  UE_cells_PATH  <- paste0(OUTPUT, 'unique_events_cells.rds')
  ue_df <- readRDS(UE_cells_PATH)
  
  # Ensure the data is numeric
  ue_df <- data.frame(lapply(ue_df[chromosomes], as.numeric))
  
  # Check for NA values and handle them (e.g., replace NAs with zeros)
  ue_df[is.na(ue_df)] <- 0
  
  all_samples_ue_chr[[i]] <- ue_df
}

# Create a dataframe with sample categories
sample_categories <- data.frame(
  sample = samples,
  category = factor(category, levels = category_order)  # Use the category_order here
)

# Order samples by category
ordered_samples <- sample_categories %>%
  arrange(category) %>%
  pull(sample)

# Reorder all_samples_ue_chr and samples
all_samples_ue_chr <- all_samples_ue_chr[match(ordered_samples, samples)]
samples <- ordered_samples

# Define a color palette with breaks and colors
color_breaks <- seq(0, max(unlist(all_samples_ue_chr), na.rm = TRUE), length.out = 100)
color_palette <- viridis(100)
col_fun <- colorRamp2(color_breaks, color_palette)

# Function to create a heatmap for a single sample
create_sample_heatmap <- function(data, sample_name, category) {
  # Perform clustering on rows
  row_clusters <- hclust(dist(data))
  row_order <- row_clusters$order
  
  # Create the heatmap
  heatmap <- Heatmap(
    data[row_order, ], 
    name = sample_name,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = paste(sample_name, "\n(", category, ")", sep = ""),
    col = col_fun
  )
  
  return(heatmap)
}

# Create a list of heatmaps
heatmap_list <- pmap(
  list(
    data = all_samples_ue_chr, 
    sample_name = samples, 
    category = sample_categories$category[match(samples, sample_categories$sample)]
  ),
  create_sample_heatmap
)

# Create row annotations with sample names
row_anno <- rowAnnotation(
  Sample = anno_text(samples, gp = gpar(fontsize = 10))
)

# Combine heatmaps into a single plot using the HeatmapList
ht_list <- Reduce(`%v%`, heatmap_list)

# Add the row annotations to the combined heatmap
ht_list <- ht_list + row_anno

# Save the plot
png("Documents/sc_analysis/all_samples_23july2024/combined_heatmap.png", width = 12 * 300, height = 4 * length(all_samples_ue_chr) * 300, res = 300)
draw(ht_list, heatmap_legend_side = "right", column_title = "Combined Heatmap of Samples")
dev.off()
