library(ggplot2)
library(dplyr)
library(forcats)


samples <- c("23003","23303","23359","23526","24077","24130","24532","24007","23965")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO4","PEO6","PEO1","PEO1","PEO1")
category_order <- c("PEO1", "PEO4","PEO6")  # Adjust this to your preferred order



all_samples_ue <- list()

# Loop to create sublists and add them to the main list
for (i in 1:length(samples)) {
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[i],'/')
  chromosomes <- c(paste0("chr", 1:22))
  UE_cells_PATH  <- paste0(OUTPUT,'unique_events_cells.rds')
  ue_df <- readRDS(UE_cells_PATH)
  ue_df_sum <- rowSums(ue_df["chr22"])
  all_samples_ue[[i]] <- ue_df_sum
}

# Assuming all_samples_ue is a list, convert it to a data frame
df <- data.frame(
  Sample = rep(samples, sapply(all_samples_ue, length)),
  Value = unlist(all_samples_ue)
)

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

# Assuming df has columns: Cell, Value, and Group (representing different groups)
custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
                   "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB")

# Step 1: Arrange the data by group and order values within each group
df_ordered <- df %>%
  group_by(group) %>%
  arrange(Value) %>%
  mutate(Cell_Order = row_number())  # Create a new column for the order of cells


# Step 2: Plot the data with dots
ggplot(df_ordered, aes(x = Cell_Order, y = Value, color = Category)) +
  geom_point(size = 2) +  # Use geom_point to plot dots, and adjust size as needed
  labs(x = "Cell Order (from lowest to highest value)", y = "Value") +
  theme_minimal() +
  theme(legend.title = element_blank()) +  # Remove legend title
  scale_color_brewer(palette = "Set1")  

