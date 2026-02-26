
library(dplyr)
library(tibble)
library(ggplot2)
library(ggforce)
library(reshape2)
signature <- read.csv("/Volumes/LenovoPS8/FIbackup/bulk_sc_siganture_quantification_thresholded_activities.copied.tsv",sep="\t")
annotation <- read.csv("/Volumes/LenovoPS8/FIbackup/Paired Cell Lines Summary Sheet - Sheet1.csv",sep=",")
annotation <- annotation[,c("SLX","Cell.Line","Cell.Line.Passage","SINCEL")]
# Remove duplicates in df based on SINCEL, keeping the first occurrence
annotation <- annotation %>% distinct(SINCEL, .keep_all = TRUE)
annotation$sample =  paste0(annotation$Cell.Line," ", annotation$Cell.Line.Passage)
# Load necessary libraries


# Assuming your matrix is stored in a variable called 'signature'
# Rename samples
#rownames(signature)
#[1] "SINCEL-193" "SINCEL-194" "SINCEL-206" "SINCEL-211" "SINCEL-212" "SINCEL-218" "SINCEL-219" "SINCEL-220" "SINCEL-221" "SINCEL-245"
#[11] "SINCEL-246" "SINCEL-249" "SINCEL-250" "SINCEL-251" "SINCEL-252" "SINCEL-257"


#SINCEL-194 FUCCI
#SINCEL-221 CIOV1

# Check for duplicates in SINCEL column
if(any(duplicated(annotation$SINCEL))) {
  stop("Duplicates found in SINCEL column of df")
}
if(any(duplicated(rownames(signature)))) {
  stop("Duplicates found in rownames of signature")
}
if(any(duplicated(annotation$sample))) {
  stop("Duplicates found in samples")
}


# Join the dataframes
joined_df <- signature %>%
  rownames_to_column("SINCEL") %>%
  left_join(annotation, by = "SINCEL")

joined_df <- na.omit(joined_df)
rownames(joined_df)<-NULL
# Update row names
signature <- joined_df %>%
  column_to_rownames("sample")

# View the updated signature dataframe
print(signature)

signature$Sample <- rownames(signature)



# Assuming your data is in a data frame called signature
# Add the group column to the DataFrame
group_mapping <- c("PEO1 P90" = "PEO1/4/6", "PEO1 P112" = "PEO1/4/6", "PEO1 P135" = "PEO1/4/6", "PEO1 STOP Px+9" = "PEO1/4/6", "PEO1 MISSENSE Px+5" = "PEO1/4/6", "PEO4 Px+5" = "PEO1/4/6", "PEO6 Px+11" = "PEO1/4/6", "CIOV3 Px+48" = "CIOV", "CIOV6 px+9" = "CIOV", "UWB1.289 Px+5" = "UWB1.289", "UWB1.289 BRCA Px+5" = "UWB1.289", "HCT116 BRCA2 -/- Px+9" = "HCT", "PEO14 P56" = "PEO14/23", "PEO23 Px+29" = "PEO14/23")
signature$Group <- group_mapping[rownames(signature)]

# Melt the data for ggplot2
signature_melted <- melt(signature, id.vars = c("Sample", "Group", "Position"), measure.vars = paste0("CX", 1:17),
                         variable.name = "Signatures")

# Assuming sample_order is a numeric vector indicating the order of samples
sample_order <- c(3, 8, 2, 5, 12, 4, 13, 10, 11, 6, 9, 14, 7, 1)

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
  labs(title = "CN Signatures", x = "Sample", y = "Value") +
  scale_fill_manual(values = random_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the plot
print(p)

# Save the plot with adjusted size
ggsave("../../Volumes/Fl/sc_analysis/all_samples_23july2024/stacked_bar_plot.png", plot = p, width = 10, height = 5)
