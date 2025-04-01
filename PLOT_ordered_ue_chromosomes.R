
library(ggplot2)
library(dplyr)
library(forcats)

samples <- c("23003","23303","23359","23526","24077","24130","24532","24007","23965")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO4","PEO6","PEO1","PEO1","PEO1")
category_order <- c("PEO1", "PEO4","PEO6")  # Adjust this to your preferred order
chromosomes <- c(paste0("chr", 1:22))

all_samples_ue <- list()

# Loop to collect data for all chromosomes
for (i in 1:length(samples)) {
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[i],'/')
  UE_cells_PATH  <- paste0(OUTPUT,'unique_events_cells.rds')
  ue_df <- readRDS(UE_cells_PATH)
  
  for (chr in chromosomes) {
    ue_df_sum <- rowSums(ue_df[chr])  # Summing across the chromosome
    all_samples_ue[[length(all_samples_ue) + 1]] <- data.frame(
      Sample = samples[i],
      Chromosome = chr,
      Value = ue_df_sum
    )
  }
}

# Combine the list into a single dataframe
df <- do.call(rbind, all_samples_ue)

# Create a mapping of Sample to Category
sample_categories <- data.frame(
  Sample = unique(df$Sample),
  Category = category
)

# Join the category information to the main dataframe
df <- df %>%
  left_join(sample_categories, by = "Sample")

df <- df %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)),
         group = interaction(Category, Sample, drop = TRUE, sep = "_"))

# Arrange the data by group and order values within each group
df_ordered <- df %>%
  group_by(group, Chromosome) %>%
  arrange(Value) %>%
  mutate(Cell_Order = row_number())  # Create a new column for the order of cells

# Plot the data with dots, keeping the same y-axis across all subplots
p=ggplot(df_ordered, aes(x = Cell_Order, y = Value, color = Category)) +
  geom_point(size = 0.5) +  # Use geom_point to plot dots, and adjust size as needed
  labs(x = "Cell Order (from lowest to highest value)", y = "Number of Unique Events") +
  theme_minimal() +
  theme(legend.title = element_blank()) +  # Remove legend title
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~ Chromosome)  # Keep the same y-axis across all subplots

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/ue_ordered_chr_PEO.pdf", 
       plot = p, width = 6, height = 6, dpi = 300)
