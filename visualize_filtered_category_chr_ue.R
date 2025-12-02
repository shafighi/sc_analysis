library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(viridis)

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1")
category_order <- c("PEO1", "PEO4", "PEO6", "CIOV3", "CIOV6", "PEO23", "HCT116\nBRCA2 -/-")

# Define desired categories
desired_categories <- c("UWB1.289", "UWB1.289\nBRCA")

# Filter samples and categories based on the desired categories
sample_categories <- data.frame(
  sample = samples,
  category = category
) %>% 
  filter(category %in% desired_categories)


# Update all_samples_ue_chr to only include desired categories
filtered_samples <- sample_categories$sample
all_samples_ue_chr <- all_samples_ue_chr[match(filtered_samples, samples)]
samples <- filtered_samples
category <- sample_categories$category

# Combine all data frames into one
all_samples_ue_chr_combined <- bind_rows(all_samples_ue_chr)
all_samples_ue_chr_combined <- all_samples_ue_chr_combined %>%
  filter(Sample %in% filtered_samples)

# Reshape the data into long format
long_data <- all_samples_ue_chr_combined %>%
  pivot_longer(cols = -Sample, names_to = "Chromosome", values_to = "Value") %>%
  left_join(sample_categories, by = c("Sample" = "sample"))

# Define colors for each category
category_colors <- viridis(length(desired_categories))
names(category_colors) <- desired_categories

# Create boxplots
boxplot_gg <- ggplot(long_data, aes(x = Chromosome, y = Value, fill = category)) +
  geom_boxplot() +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  labs(title = "Boxplots of Chromosomes for Selected Categories",
       x = "Chromosome",
       y = "Value",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/combined_boxplot_filtered_UWB1.289.png", plot = boxplot_gg, width = 15, height = 10)

