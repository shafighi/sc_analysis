library(ggplot2)
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
  
  ue_df$Sample <- samples[i]
  all_samples_ue_chr[[i]] <- ue_df
}

# Combine all data frames into one
combined_data <- bind_rows(all_samples_ue_chr)

# Reshape the data into long format
long_data <- combined_data %>%
  pivot_longer(cols = -Sample, names_to = "Chromosome", values_to = "Value")

# Assign unique colors to each sample
sample_colors <- viridis(length(samples))
names(sample_colors) <- samples

# Create boxplots
boxplot_gg <- ggplot(long_data, aes(x = Chromosome, y = Value, fill = Sample)) +
  geom_boxplot() +
  scale_fill_manual(values = sample_colors) +
  theme_minimal() +
  labs(title = "Boxplots of Chromosomes Across Samples",
       x = "Chromosome",
       y = "Value",
       fill = "Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/combined_boxplot.png", plot = boxplot_gg, width = 25, height = 8)







