library(dplyr)
library(knitr)
library(gt)
library(webshot2)
library(kableExtra)

# Define the samples and categories
ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24491","24518","24835")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116 BRCA2 -/-","PEO14","UWB1.289 BRCA","PEO1","PEO1","UWB1.289","PEO1 MISSENSE","PEO1 STOP","PEO1")
HRD_HRP <- c("HRD","HRD","HRD","HRP","HRP","HRP","HRP","HRP","HRD","HRP","HRP","HRD","HRD","HRD","HRD","HRP","HRD","HRD","HRD","HRD","HRD","HRD")
Resistance <- c("Cisplatin sensitive","Cisplatin sensitive","Cisplatin sensitive","Cisplatin Resistant","Cisplatin Resistant","Cisplatin Resistant","-","Cisplatin Resistant","-","Cisplatin Resistant","BRCA2 WT","-","Cisplatin Sensitive","BRCA2 null","Cisplatin Sensitive","BRCA1 restored","Cisplatin Sensitive","Cisplatin Sensitive","BCRA1 mutation","-","-","Cisplatin Sensitive")

bin_size <- "100"

# Initialize an empty data frame to store the combined outlier summaries
combined_outlier_summary <- data.frame()

for (i in 1:length(ALLSAMPLES)) {
  SAMPLENAME <- ALLSAMPLES[i]
  SAMPLEONJ <- paste0("SLX-",SAMPLENAME,"_",bin_size)
  OUTPUT <- file.path("Documents/sc_analysis/post-scAbsolute", SAMPLEONJ)
  
  # Read the outlier_summary for the current sample
  outlier_summary <- readRDS(paste0(OUTPUT, "/outlier_summary.rds"))
  
  # Add columns for the sample name, category, HRD_HRP, and Resistance
  outlier_summary$Sample <- SAMPLENAME
  outlier_summary$Category <- categories[i]
  outlier_summary$HRD_HRP <- HRD_HRP[i]
  outlier_summary$Resistance <- Resistance[i]
  
  # Combine the current outlier summary with the accumulated data
  combined_outlier_summary <- rbind(combined_outlier_summary, outlier_summary)
}

# Define the custom order for categories
category_order <- c("PEO1", "PEO4", "PEO6", "PEO1 MISSENSE", "PEO1 STOP", "PEO14", "PEO23", "CIOV3", "CIOV6", "HCT116", "HCT116 BRCA2 -/-", "UWB1.289", "UWB1.289 BRCA")

# Convert the Category column to a factor with specified levels
combined_outlier_summary$Category <- factor(combined_outlier_summary$Category, levels = category_order)

# Order the dataframe by Category
sorted_data <- combined_outlier_summary %>%
  arrange(Category) %>%
  select(Sample, Category, HRD_HRP, Resistance, everything())

# Calculate the sum for numeric columns
sum_row <- sorted_data %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(Sample = "Total", Category = "", HRD_HRP = "", Resistance = "")

# Add the sum row to the dataframe
sorted_data <- bind_rows(sorted_data, sum_row)

kable(sorted_data, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  add_header_above(c(" " = 4, "Datasets" = ncol(sorted_data) - 4)) %>%
  row_spec(nrow(sorted_data), bold = TRUE, background = "#EFEFEF") %>%
  column_spec(4, width = "150px", extra_css = "white-space: nowrap;") %>%  # Adjust the Resistance column
  save_kable(file = "Documents/sc_analysis/all_samples_23july2024/outlier_summary_table_meta.png", zoom = 2)
