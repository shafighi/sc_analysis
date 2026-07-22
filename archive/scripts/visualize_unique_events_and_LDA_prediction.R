library(ggplot2)
library(dplyr)
library(forcats)

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA")  # Adjust this to your preferred order

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24831","24832","24833","25146","25147","25148","25149")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA","FFPE")  # Adjust this to your preferred order


samples <- c("23003","23303","23359","23526","24077","24130","24532","24007","23965")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO4","PEO6","PEO1","PEO1","PEO1")
category_order <- c("PEO1", "PEO4","PEO6")  # Adjust this to your preferred order


samples <- c("23303","23526","24077","24130")
category <- c("PEO1","PEO4","PEO4","PEO6")
category_order <- c("PEO1", "PEO4","PEO6")  # Adjust this to your preferred order

samples <- c("24175","24174")
category <- c("HCT116", "HCT116\nBRCA2 -/-")
category_order <- c("HCT116", "HCT116\nBRCA2 -/-")  # Adjust this to your preferred order


samples <- c("23961","24078","24441")
category <- c("CIOV3", "CIOV6","CIOV6")
category_order <- c("CIOV3", "CIOV6")  # Adjust this to your preferred order


all_samples_ue <- list()

# Loop to create sublists and add them to the main list
for (i in 1:length(samples)) {
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[i],'/')
  UE_cells_PATH  <- paste0(OUTPUT,'unique_events_cells_compact.rds')
  ue_df <- readRDS(UE_cells_PATH)
  ue_df$sample <- samples[i]  # Add sample ID to each dataframe
  ue_df$category <- category[i] 
  all_samples_ue[[i]] <- ue_df  # Store dataframe in the list
}

# Concatenate all dataframes into one
combined_df <- bind_rows(all_samples_ue)

# Summarize the data: count rows per node and sample, calculate the mean of 'length'
result_df <- combined_df %>%
  group_by(node, sample,category) %>%
  summarise(
    count = n(),  # Count rows per node and sample
    mean_length = mean(length, na.rm = TRUE)  # Calculate mean length
  ) %>%
  ungroup()

# View the result
print(result_df)

df_scaled <- result_df %>%
  mutate(feature1_scaled = scale(count),
         feature2_scaled = scale(mean_length))

library(viridis)  # For a colorblind-friendly palette

# Create the scatter plot with viridis colors
#p <- ggplot(result_df[result_df$count<150 & result_df$mean_length<150,], aes(x = count, y = mean_length, color = category)) +
#p <- ggplot(result_df, aes(x = count, y = mean_length, color = category)) +
p <- ggplot(df_scaled, aes(x = feature1_scaled[,1], y = feature2_scaled[,1], color = category)) +
  geom_point(size = 1) +  # Scatter plot with points
  labs(
    title = "Unique Events vs. Length by Sample",
    x = "Number of Unique Events",
    y = "Mean Length of Unique Events"
  ) +
  theme_minimal() +  # Clean theme
  theme(
    legend.title = element_text(size = 10),  # Adjust legend title size
    legend.text = element_text(size = 8)     # Adjust legend text size
  ) +
  scale_color_viridis_d(option = "plasma")  # Use distinct viridis color palette

p
# Save the plot to a file
ggsave("Documents/sc_analysis/all_samples_23july2024/unique_events_plot.png", plot = p, width = 8, height = 6, dpi = 300)


library(MASS)

# Perform LDA
lda_model <- lda(category ~ count + mean_length, data = result_df)

# Get the LDA predictions
lda_pred <- predict(lda_model)

# Add the LDA predictions to the dataframe
result_df$LD1 <- lda_pred$x[, 1]  # First linear discriminant
result_df$LD2 <- lda_pred$x[, 2]  # Second linear discriminant (if available)

# Plot the LDA results
ggplot(result_df, aes(x = LD1, y = LD2, color = category)) +
  geom_point(size = 1) +  # Scatter plot with points
  labs(
    title = "LDA Plot: Separation of Samples",
    x = "LD1 (First Linear Discriminant)",
    y = "LD2 (Second Linear Discriminant)"
  ) +
  theme_minimal() +  # Clean theme
  scale_color_brewer(palette = "Set1")  # Optional color palette
