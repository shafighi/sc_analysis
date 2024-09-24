#install.packages("kableExtra")
#install.packages("webshot")
#install.packages("htmltools")
#install.packages("rmarkdown")
library(kableExtra)
library(webshot)
library(htmltools)

library(dplyr)
library(knitr)
library(gt)
library(webshot2)
library(kableExtra)

# Define the samples and categories
ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24491","24518","24835")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","PEO1\nMISSENSE","PEO1\nSTOP","PEO1")


# Define the samples and categories
ALLSAMPLES <- c("23355","25146","24831","24832","24833")
ALLSAMPLES=c("25146","25147","25149","24831","24832","24833")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
ALLSAMPLES=c("24518","24491")
categories <- c("PEO STOP","PEO MISSENSE")

ALLSAMPLES=c("24911","24912","24913")
categories <- c("FFPE","FFPE","FFPE")


bin_size <- "100"

# Initialize an empty data frame to store the combined outlier summaries
combined_outlier_summary <- data.frame()

for (i in 1:length(ALLSAMPLES)) {
  SAMPLENAME <- ALLSAMPLES[i]
  SAMPLEONJ <- paste0("SLX-",SAMPLENAME,"_",bin_size)
  OUTPUT <- file.path("Documents/sc_analysis/post-scAbsolute", SAMPLEONJ)
  
  # Read the outlier_summary for the current sample
  outlier_summary <- readRDS(paste0(OUTPUT, "/outlier_summary.rds"))
  
  # Add columns for the sample name and category
  outlier_summary$Sample <- SAMPLENAME
  outlier_summary$Category <- categories[i]
  
  # Combine the current outlier summary with the accumulated data
  combined_outlier_summary <- rbind(combined_outlier_summary, outlier_summary)
}

# Optionally, set the row names to the sample names if desired
# rownames(combined_outlier_summary) <- combined_outlier_summary$Sample

# Define the custom order for categories
#categories <- c("PEO1", "PEO4", "PEO6", "PEO1\nMISSENSE", "PEO1\nSTOP", "PEO14", "PEO23", "CIOV3", "CIOV6", "HCT116", "HCT116\nBRCA2 -/-", "UWB1.289", "UWB1.289\nBRCA")

# Convert the Category column to a factor with specified levels
####xf$Category <- factor(combined_outlier_summary$Category, levels = categories)

# Order the dataframe by Category
sorted_data <- combined_outlier_summary %>%
  arrange(Category) %>%
  select(Sample, Category, everything())

# Calculate the sum for numeric columns
sum_row <- sorted_data %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(Sample = "Total", Category = "")

# Add the sum row to the dataframe
sorted_data <- bind_rows(sorted_data, sum_row)

# Create and save the table as PNG
kable(sorted_data, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  add_header_above(c(" " = 2, "Datasets" = ncol(sorted_data) - 2)) %>%
  row_spec(nrow(sorted_data), bold = TRUE, background = "#EFEFEF") %>%
  save_kable(file = "Documents/sc_analysis/all_samples_23july2024/outlier_summary_table_FFPE.png", zoom = 2)

