# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)


df<- read.csv("../../Volumes/Fl/sc_analysis/all_samples_23july2024/Cellines_95.csv")
df <- df %>%
  mutate(Sequenced = c(288,288, 288, 134, 384,384,384 ,384, 384,384,384 ,384,384,384,384,288,384,384,384,384,288,384,7718))  # Adding a new column


df <- df %>%
  mutate(High_Quality_Percentage = (Good.Quality.Cells+Replicating)*100/Sequenced)  # Adding a new column


ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24491","24518","24835")
categories <- c("PEO1(23003)","PEO1(23303)","PEO1(23359)","PEO4(23526)","PEO23(23527)","PEO23(23528)","CIOV3(23961)","PEO4(24077)","CIOV6(24078)","PEO6(24130)","HCT116(24175)","CIOV6(24441)","PEO1(24532)","HCT116\nBRCA2 -/-(24174)","PEO14(24173)","UWB1.289\nBRCA(24489)","PEO1(24007)","PEO1(23965)","UWB1.289(24490)","PEO1 MISSENSE(24491)","PEO1 STOP(24518)","PEO1(24835)")

# Creating a named vector from ALLSAMPLES and categories
sample_category_map <- setNames(categories, ALLSAMPLES)

# Adding the category to the dataframe
df$Cell_line <- sample_category_map[df$Sample]

write.csv(df,"../../Volumes/Fl/sc_analysis/all_samples_23july2024/Cellines_95_2.csv")
# Create the plot
p <- ggplot(df, aes(x = High_Quality_Percentage)) +
  geom_density(fill = "skyblue", color = "navy", alpha = 0.7, linewidth = 1.5) +
  theme_minimal(base_size = 12) +  # Reduce base font size
  labs(title = "Distribution of Good-Quality Cell Percentage\n in Cell Lines", 
       x = "Percentage of Good Quality Cells", 
       y = "Density") +
  scale_x_continuous(breaks = seq(10, 100, by = 5)) +  # More x-axis ticks
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01)) +  # More y-axis ticks
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12),   # Smaller tick labels
        axis.title = element_text(face = "bold", size = 12))  # Smaller axis labels

p
# Save as a small PDF (adjust size as needed)
ggsave("../../Volumes/Fl/sc_analysis/all_samples_23july2024/Good_Quality_Distribution.pdf", plot = p, width = 5, height = 4, units = "in") 


df$High_Quality_Percentage = round(df$High_Quality_Percentage)
df$Filtered = df$Good.Quality.Cells - df$Normal.Cells
colnames(df) = c( "Processed.Cells","Replicating","Replicating...RPC" ,"RPC.Outliers","Alpha.Mapd.Gini","Good.Quality","Normal","Sample","Category","Sequenced","High_Quality_Percentage","Cell_line","Filtered")
# Reshape the data to long format for ggplot
df_long <- df[1:22,] %>%
  pivot_longer(cols = c(Replicating, Normal, Filtered),
               names_to = "Good.Quality.Cells", values_to = "Count")

colors <- c(
  "Replicating" = "#81C784",   # Soft green
  "Filtered" = "#D1C4E9",      # Soft lavender
  "Normal" = "#F0E68C"         # A softer, muted yellow (Khaki)
)
p2 <- ggplot(df_long, aes(x = Cell_line, y = Count, fill = Good.Quality.Cells)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  scale_fill_manual(values = colors) +  # Custom colors
  geom_text(data = df[1:22,], aes(x = Cell_line, y = Sequenced + 5, label = paste0(High_Quality_Percentage, "%")),
            size = 3, fontface = "bold", angle = 90, inherit.aes = FALSE) +  # Increase font size & rotate
  theme_minimal(base_size = 12) +  # Clean white background with slightly bigger base font
  labs(title = "Cell Composition by Cell line", x = "Cell line", y = "Cell Count", fill = "Good Quality Cells") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 1,size = 8))
p2
# Save the figure with a smaller height (e.g., 5) and a reasonable width (e.g., 8)
ggsave("../../Volumes/Fl/sc_analysis/all_samples_23july2024/cell_composition_plot.pdf", p2, width = 6, height = 7)



# FFPE


df_ffpe<- read.csv("../../Volumes/Fl/sc_analysis/all_samples_23july2024/FFPE_95.csv")
df_ffpe <- df_ffpe %>%
  mutate(Sequenced = c(384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384,7296))  # Adding a new column


df_ffpe <- df_ffpe %>%
  mutate(High_Quality_Percentage = round((Good.Quality.Cells+Replicating)*100/Sequenced))  # Adding a new column

write.csv(df_ffpe,"../../Volumes/Fl/sc_analysis/all_samples_23july2024/FFPE_95_2.csv")

# Create the plot
p3 <- ggplot(df_ffpe, aes(x = High_Quality_Percentage)) +
  geom_density(fill = "skyblue", color = "navy", alpha = 0.7, linewidth = 1.5) +
  theme_minimal(base_size = 12) +  # Reduce base font size
  labs(title = "Distribution of Good-Quality Cell Percentage\n in FFPE", 
       x = "Percentage of Good Quality Cells", 
       y = "Density") +
  scale_x_continuous(breaks = seq(10, 100, by = 5)) +  # More x-axis ticks
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01)) +  # More y-axis ticks
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12),   # Smaller tick labels
        axis.title = element_text(face = "bold", size = 12))  # Smaller axis labels

p3
# Save as a small PDF (adjust size as needed)
ggsave("../../Volumes/Fl/sc_analysis/all_samples_23july2024/Good_Quality_Distribution_FFPE.pdf", plot = p3, width = 5, height = 4, units = "in") 


df_ffpe$High_Quality_Percentage = round(df_ffpe$High_Quality_Percentage)
df_ffpe$Filtered = df_ffpe$Good.Quality.Cells - df_ffpe$Normal.Cells
colnames(df_ffpe) = c( "Processed.Cells","Replicating","Replicating...RPC","RPC","Alpha.Mapd.Gini","Good.Quality","Normal","Sample","Category","Sequenced","High_Quality_Percentage","Filtered")
# Reshape the data to long format for ggplot
df_long <- df_ffpe[1:19,] %>%
  pivot_longer(cols = c(Replicating, Normal, Filtered),
               names_to = "Kind", values_to = "Count")

colors <- c(
  "Replicating" = "#81C784",   # Soft green
  "Filtered" = "#D1C4E9",      # Soft lavender
  "Normal" = "#F0E68C"         # A softer, muted yellow (Khaki)
)
p4 <- ggplot(df_long, aes(x = Sample, y = Count, fill = Kind)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  scale_fill_manual(values = colors) +  # Custom colors
  geom_text(data = df_ffpe[1:19,], aes(x = Sample, y = Sequenced + 5, label = paste0(High_Quality_Percentage, "%")),
            size = 3, fontface = "bold", angle = 90, inherit.aes = FALSE) +  # Increase font size & rotate
  theme_minimal(base_size = 12) +  # Clean white background with slightly bigger base font
  labs(title = "Cell Composition in FFPE Samples", x = "FFPE", y = "Cell Count", fill = "Good Quality Cells") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 1,size = 8))
p4
# Save the figure with a smaller height (e.g., 5) and a reasonable width (e.g., 8)
ggsave("../../Volumes/Fl/sc_analysis/all_samples_23july2024/cell_composition_plot_FFPE.pdf", p4, width = 6, height = 7)




library(gridExtra)

pdf("../../Volumes/Fl/sc_analysis/all_samples_23july2024/Good_Quality_Distribution_Combined.pdf", width = 10, height = 4)
grid.arrange(p, p3, ncol = 2)
dev.off()

pdf("../../Volumes/Fl/sc_analysis/all_samples_23july2024/Cell_Composition_Combined.pdf", width = 10, height = 5)
grid.arrange(p2, p4, ncol = 2)
dev.off()


