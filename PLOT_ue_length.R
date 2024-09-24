library(ggplot2)
library(dplyr)
library(forcats)

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA")  # Adjust this to your preferred order

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24831","24832","24833","25146","25147","25148","25149")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA","FFPE")  # Adjust this to your preferred order



all_samples_ue <- list()

# Loop to create sublists and add them to the main list
for (i in 1:length(samples)) {
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[i],'/')
  chromosomes <- c(paste0("chr", 1:22))
  UE_cells_PATH  <- paste0(OUTPUT,'unique_events_cells_compact.rds')
  ue_df <- readRDS(UE_cells_PATH)
  #ue_df_sum <- rowSums(ue_df[chromosomes])
  all_samples_ue[[i]] <- ue_df$length
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
threshold=500


# Create the plot
p <- ggplot(df, aes(x = group, y = Value, fill = Category)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5) +
  theme_classic(base_size = 7) +
  labs(title = "Unique events",
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
           y = min(df$Value) - 0.1 * diff(range(df$Value)),
           label = levels(df$Category),
           size = 3,
           fontface = "bold")

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/ue_lenght_boxplot_ffpe_included.png", 
       plot = p, width = 14, height = 6, dpi = 300)



# Create the plot
p <- ggplot(df[df$Value>0,], aes(x = group, y = Value, fill = Category)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5) +
  theme_classic(base_size = 7) +
  labs(title = "Unique events",
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
           y = min(df$Value) - 0.1 * diff(range(df$Value)),
           label = levels(df$Category),
           size = 3,
           fontface = "bold")

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/ue_lenght_boxplot_ffpe_nonzero.png", 
       plot = p, width = 14, height = 6, dpi = 300)




threshold=100
# Create the plot
p <- ggplot(df[df$Value<threshold,], aes(x = group, y = Value, fill = Category)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5) +
  theme_classic(base_size = 7) +
  labs(title = "Unique events",
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
           y = min(df$Value) - 0.1 * diff(range(df$Value)),
           label = levels(df$Category),
           size = 3,
           fontface = "bold")

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/ue_lenght_boxplot_ffpe_included_small_500.png", 
       plot = p, width = 14, height = 6, dpi = 300)


p <- ggplot(df[df$Value>threshold,], aes(x = group, y = Value, fill = Category)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5) +
  theme_classic(base_size = 7) +
  labs(title = "Unique events",
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
           y = min(df$Value) - 0.1 * diff(range(df$Value)),
           label = levels(df$Category),
           size = 3,
           fontface = "bold")

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/ue_lenght_boxplot_ffpe_included_big_500.png", 
       plot = p, width = 14, height = 6, dpi = 300)

