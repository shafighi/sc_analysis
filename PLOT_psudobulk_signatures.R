


library(dplyr)
library(tibble)
library(ggplot2)
library(ggforce)
library(reshape2)
signature <- readRDS("../../Volumes/Fl/sc_pseudo_rerun_30kb/sc_pseudo_30kb_ds_absCopyNumber_acts.rds")

if(any(duplicated(rownames(signature)))) {
  stop("Duplicates found in rownames of signature")
}
signature <- data.frame(signature)
signature$Sample <- rownames(signature)




# Assuming your data is in a data frame called signature
# Add the group column to the DataFrame
group_mapping <- c("PEO1_2N" = "PEO1/4/6", "PEO1_FUCCI" = "PEO1/4/6", "PEO1_missense" = "PEO1/4/6", "PEO1_stop" = "PEO1/4/6", "PEO4_smallbam" = "PEO1/4/6", "PEO4_bigbam" = "PEO1/4/6", "PEO6_smallbam" = "PEO1/4/6", "PEO6_bigbam" ="PEO1/4/6", "CIOV3" = "CIOV", "CIOV6" = "CIOV","PEO14"="PEO14/23","PEO23"="PEO14/23")
signature$Group <- group_mapping[rownames(signature)]



# Define the mapping of old names to new names
name_mapping <- c("PEO1_2N" = "PEO1 p135", "PEO1_FUCCI" = "PEO1 FUCCI px+27", "PEO1_missense" = "PEO1 Missense Px+5", "PEO1_stop" = "PEO1 Stop Px+9", 
                  "PEO4_smallbam" = "PEO4 Px+8", "PEO4_bigbam" = "PEO4 (Px+8) + (Px+5)", "PEO6_smallbam" = "PEO6 px+11", "PEO6_bigbam" = "PEO6 2(px+11)", 
                  "CIOV3" = "CIOV3 Px+48", "CIOV6" = "CIOV6 px+9", "PEO14" = "PEO14 2(P56)", "PEO23" = "PEO23 2(Px+29)")

# Replace the sample names using the mapping
signature <- signature %>%
  mutate(Sample = recode(Sample, !!!name_mapping))



# Melt the data for ggplot2
signature_melted <- melt(signature, id.vars = c("Sample", "Group"), measure.vars = paste0("CX", 1:17),
                         variable.name = "Signatures")

# Assuming sample_order is a numeric vector indicating the order of samples
sample_order <- c(11,12 , 1, 2, 3, 4,9,10, 5, 6, 7, 8)

# Add sample_order to the data frame
signature_melted$sample_order <- sample_order

# Convert the Sample column to a factor with the specified order
signature_melted <- signature_melted %>%
  arrange(sample_order) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Calculate the number of bars in each group
bar_counts <- signature_melted %>%
  group_by(Group) %>%
  summarise(bar_count = n())

# Merge the bar counts back into the original data
signature_melted <- signature_melted %>%
  left_join(bar_counts, by = "Group")

pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", 
                   "#D4A5A5", "#A5D4A5", "#A5A5D4", "#D4A5D4", "#A5D4D4",
                   "#D4D4A5", "#F3C1C6", "#C1F3C6", "#C1C6F3", "#F3C6C1",
                   "#C6F3C1", "#C6C1F3")

# Define the dark colors
dark_colors <- c("#8B0000", "#006400", "#00008B", "#FF8C00", "#4B0082", 
                 "#008B8B", "#8B008B", "#556B2F", "#483D8B", "#2F4F4F",
                 "#B8860B", "#BDB76B", "#8FBC8F", "#00CED1", "#9400D3",
                 "#E9967A", "#9932CC")

light_colors <- c("#F08080", "#FFA07A", "#FFB6C1", "#90EE90", "#87CEFA", 
                  "#FAFAD2", "#B0C4DE", "#E0FFFF", "#20B2AA", "#778899",
                  "#ADD8E6", "#FFFFE0", "#D3D3D3", "#FF9999", "#C5CBE1",
                  "#E6E6FA", "#FFDAB9")
fun_colors <- c("#7DF9FF", "#39FF14", "#FF69B4", "#FFA500", "#32CD32", 
                "#8F00FF", "#FFD700", "#40E0D0", "#FF00FF", "#FF0000",
                "#00FFFF", "#7FFF00", "#FF7F50", "#FF00FF", "#FFD700",
                "#00FF7F", "#DC143C")

lighter_fun_colors <- c("#ADD8E6", "#98FF98", "#FFB6C1", "#FFDAB9", "#E6E6FA", 
                        "#87CEEB", "#FFFACD", "#00FFFF", "#F08080", "#AFEEEE",
                        "#FFA07A", "#FAFAD2", "#98FB98", "#FFB6C1", "#E0FFFF",
                        "#20B2AA", "#87CEFA")

# Define the random colors
random_colors <- c("#DC143C", "#00BFFF", "#32CD32", "#FF69B4", "#FFD700", 
                   "#BA55D3", "#FF6347", "#1E90FF", "#00FF7F", "#FF4500",
                   "#40E0D0", "#EE82EE", "#7FFF00", "#FF7F50", "#DA70D6",
                   "#9ACD32", "#4682B4")


# Create the stacked bar plot
p <- ggplot(signature_melted, aes(x = Sample, y = value, fill = Signatures)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free") +
  coord_flip() +
  theme_minimal() +
  labs(title = "CN Signatures (pseudo-bulk)", x = "Sample", y = "Value") +
  scale_fill_manual(values = random_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the plot
print(p)

# Save the plot with adjusted size
ggsave("../../Volumes/Fl/sc_analysis/all_samples_23july2024/stacked_bar_plot_psudo.png", plot = p, width = 10, height = 5)

