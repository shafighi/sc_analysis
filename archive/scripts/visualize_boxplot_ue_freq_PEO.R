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
  ue_df_sum <- rowSums(ue_df[chromosomes])
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


custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
                   "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB")
# Create the plot
p <- ggplot(df[df$Value<200,], aes(x = group, y = Value, fill = Category)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5) +
  theme_classic(base_size = 7) +
  labs(title = "Frequency of Unique Events",
       x = "Samples",
       y = "Unique events in each sample") +
  scale_fill_manual(values = custom_colors) +
  #scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels = function(x) sub("^[^_]+_", "", x),
                   expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    axis.text.x.bottom = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0))
  ) +
  # Add category labels below sample names
  annotate("text", 
           x = tapply(1:nlevels(df$group), df$Category[match(levels(df$group), df$group)], mean),
           y = min(df$Value) - 0.01 * diff(range(df$Value)),
           label = levels(df$Category),
           size = 1.5)

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/ue_boxplot_PEO146.png", 
       plot = p, width = 3, height = 3, dpi = 300)

