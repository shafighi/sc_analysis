library("reshape")
library(dplyr)
library(ComplexHeatmap)
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("Documents/scUnique/R/visualize.R")


CORES=7
RMNORM=TRUE
samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24491","24518","24835")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","PEO1\nMISSENSE","PEO1\nSTOP","PEO1")

bin_size <- "100"
s=1

# Initialize an empty DataFrame for storing sample information
sample_summary <- data.frame(
  SLXNumber = character(),
  SampleName = character(),
  MeanReads = numeric(),
  MedianReads = numeric(),
  NumCells = numeric(),
  stringsAsFactors = FALSE
)

# Initialize an empty array for total reads
combined_reads <- c()

# Loop through samples
for (s in 1:length(samples)) {
  print(paste0("SAMPLE NAME: ", samples[s]))
  print("------------------ SETTING FILE NAMES ------------------")
  
  # File paths
  OUTPUT <- paste0('/Volumes/Fl/sc_analysis1/scUnique-obj/scAbsolute/SLX-', samples[s], '/')
  OBJ_PATH <- paste0('/Volumes/Fl/sc_analysis1/scAboslute-obj/SLX-', samples[s], "_100.rds")
  
  # Create directory if not exists
  if (file.exists(OUTPUT)) {
    cat("Directory is already there!", OUTPUT, "\n")
  } else {
    dir.create(OUTPUT, recursive = TRUE)
    cat("Directory is created:", OUTPUT, "\n")
  }
  
  # Read the RDS object
  object <- readRDS(OBJ_PATH)
  object_df <- Biobase::pData(object)
  
  # Calculate statistics for this sample
  mean_reads <- mean(object_df$total.reads, na.rm = TRUE)
  median_reads <- median(object_df$total.reads, na.rm = TRUE)
  num_cells <- nrow(object_df)
  
  # Add this sample's stats to the DataFrame
  sample_summary <- rbind(
    sample_summary,
    data.frame(
      SLXNumber = samples[s],
      SampleName = categories[s],
      MeanReads = mean_reads,
      MedianReads = median_reads,
      NumCells = num_cells
    )
  )
  
  # Add this sample's reads to the combined array
  combined_reads <- c(combined_reads, object_df$total.reads)
}

# Print the resulting DataFrame and combined array
print(sample_summary)
