library(dplyr)
library(knitr)
library(gt)
library(webshot2)
library(kableExtra)
library(ggplot2)
library(tidyr)
library(viridis)

# Define the samples and categories
ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24835")
sample_categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","PEO1")
category_order <- c("PEO1", "PEO4", "PEO6", "PEO14", "PEO23", "CIOV3", "CIOV6", "HCT116", "HCT116\nBRCA2 -/-", "UWB1.289", "UWB1.289\nBRCA")


ALLSAMPLES=c("25146","25147","25149","24831","24832","24833")
category_order <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
sample_categories <- c("FFPE")

ALLSAMPLES=c("24911","24912","24913")
category_order <- c("FFPE","FFPE","FFPE")
sample_categories <- c("FFPE")

bin_size <- "100"

# Initialize an empty data frame to store the combined outlier summaries
dropout_summary <- data.frame()

for (i in 1:length(ALLSAMPLES)) {
  SAMPLENAME <- ALLSAMPLES[i]
  SAMPLEONJ <- paste0("SLX-",SAMPLENAME,"_",bin_size)
  OUTPUT <- file.path("Documents/sc_analysis/post-scAbsolute", SAMPLEONJ)
  
  # Read the outlier_summary for the current sample
  cn_binned <- readRDS(paste0(OUTPUT, "/cn_binned.rds"))
  dropout <- as.data.frame(table(cn_binned[cn_binned$segVal==0,]$chromosome))
  dropout$Freq_per_cell <- dropout$Freq/length(unique(cn_binned$sample))
  dropout$Sample <- SAMPLENAME
  
  # Combine the current outlier summary with the accumulated data
  dropout_summary <- rbind(dropout_summary, dropout)
}

# Create a complete set of chromosomes for each sample
all_combinations <- expand.grid(Sample = ALLSAMPLES, Var1 = as.character(1:22))

# Merge with the existing data and fill in missing values
dropout_summary_complete <- all_combinations %>%
  left_join(dropout_summary, by = c("Sample", "Var1")) %>%
  mutate(
    Freq = ifelse(is.na(Freq), 0, Freq),
    Freq_per_cell = ifelse(is.na(Freq_per_cell), 0, Freq_per_cell)
  )

# Create a dataframe linking samples to categories
sample_category_df <- data.frame(Sample = ALLSAMPLES, Category = sample_categories)

# Merge categories with dropout summary and prepare for plotting
heatmap_data <- dropout_summary_complete %>%
  left_join(sample_category_df, by = "Sample") %>%
  mutate(
    #Category = factor(Category, levels = category_order),
    Sample = factor(Sample, levels = ALLSAMPLES),
    Var1 = factor(Var1, levels = as.character(1:22))
  ) %>%
  arrange(Category, Sample, Var1)

# Create the heatmap
ggplot(heatmap_data, aes(x = Var1, y = Sample, fill = Freq_per_cell)) +
  geom_tile() +
  scale_fill_viridis(option = "D", name = "Frequency\nper cell") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Chromosome",
    y = "Sample",
    title = "Dropout Frequency Heatmap",
    subtitle = "Frequency per cell across chromosomes and samples, grouped by category"
  )

# Save the plot
ggsave("Documents/sc_analysis/all_samples_23july2024/summary_dropout_heatmap_FFPE.png", width = 6, height = 12, dpi = 300)

