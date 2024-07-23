
library("reshape")
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
library(ggplot2)
library(viridis)
library(dendextend)

source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/scUnique-main/R/visualize.R")


# Function to split row names into chromosome, start, and end
split_row_names <- function(df) {
  parts <- strsplit(rownames(df), '[:-]')
  data.frame(
    chromosome = sapply(parts, '[', 1),
    start = as.integer(sapply(parts, '[', 2)),
    end = as.integer(sapply(parts, '[', 3))
  )
}

get_summary_of_dropouts <- function(object,cellbased_outliers,CN_BINNED_PATH,CELL_POS_DROPOUT_PATH,CN_SEGMENTED_PATH,CELL_CHR_DROPOUT_LONG_PATH){
  cn <- object@assayData$copynumber
  
  # Split row names and combine with original dataframe
  cn <- cbind(split_row_names(cn), cn)
  
  # Melt to long format
  cn_binned <- melt(cn, id.vars = c("chromosome", "start", "end"))
  
  # Extract segVal and sample 
  cn_binned$segVal <- cn_binned$value
  cn_binned$sample <- cn_binned$variable
  
  # Reorder columns
  cn_binned <- cn_binned[c("chromosome", "start", "end", "segVal", "sample")]
  cn_binned <- na.omit(cn_binned)
  saveRDS(cn_binned,CN_BINNED_PATH)
  
  # Group by Sample_Name, Chromosome, and Value, and then summarize the data
  
  cn_segmented <- cn_binned %>%
    group_by(sample, chromosome, segVal) %>%
    summarize(start = min(start), end = max(end)) %>%
    ungroup()
  
  saveRDS(cn_binned,CN_SEGMENTED_PATH)
  
  # Split the sample column into three columns: SINCEL, Plate, and AZ
  cn_binned_split <- cn_binned %>%
    separate(sample, into = c("SINCEL","210","SC", "Plate", "Plate_number", "AZNumber"), sep = "_|(?<=[0-9])(?=[A-Z])")
  
  
  
  
  
  # Split the AZ column into A and Number
  #cn_binned_split <- cn_binned_split %>%
  #  separate(AZNumber, into = c("AZ", "Number"), sep = "(?<=[A-Za-z])(?=[0-9])")
  
  # Function to split AZNumber
  split_AZNumber <- function(x) {
    az_part <- gsub("([A-Za-z]+).*", "\\1", x)
    number_part <- gsub("[A-Za-z]+", "", x)
    return(list(AZ = az_part, Number = number_part))
  }
  
  # Apply the function
  result <- split_AZNumber(cn_binned_split$AZNumber)
  
  # Add new columns to the dataframe
  cn_binned_split$AZ <- result$AZ
  cn_binned_split$Number <- result$Number
  
  
  cn_binned_split$dropout=(cn_binned_split$segVal==0)*1
  
  # Group by the newly created columns (A and Number) and summarize segVal
  cn_binned_summary <- cn_binned_split %>%
    group_by(AZ, Number) %>%
    summarize(total_dropout = sum(dropout), .groups = 'drop')
  
  
  # Set factor levels again in df_long to ensure correct order
  cn_binned_summary$Number <- factor(cn_binned_summary$Number, levels = as.character(1:24)) # Explicit order
  cn_binned_summary$AZ <- factor(cn_binned_summary$AZ, levels = rev(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")))
  
  saveRDS(cn_binned_summary,CELL_POS_DROPOUT_PATH)
  
  cn_binned_split$cell_name <- paste(cn_binned_split$AZ, cn_binned_split$Number, sep = "")
  
  
  # First, let's create a function to assign the status
  assign_status <- function(cell_name,cellbased_outliers) {
    if (cell_name %in% cellbased_outliers$cellid[cellbased_outliers$replicating == 1]) {
      return("replicating")
    } else if (cell_name %in% cellbased_outliers$cellid[cellbased_outliers$rpc_outlier == 1]) {
      return("rpc_outlier")
    } else if (cell_name %in% cellbased_outliers$cellid[cellbased_outliers$alpha_outlier == 1]) {
      return("alpha_outlier")
    } else if (cell_name %in% cellbased_outliers$cellid[cellbased_outliers$alpha_outlier == 1]) {
      return("gini_outlier")
    } else if (cell_name %in% cellbased_outliers$cellid[cellbased_outliers$alpha_outlier == 1]) {
      return("dmap_outlier")
    }else if (cell_name %in% cellbased_outliers$cellid[cellbased_outliers$na_alpha_outliers == 1]) {
      return("na_alpha_outlier")
    }else {
      return("pass")
    }
  }
  
  # Step 2: Count occurrences where segVal == 0
  cell_chr_dropout <- cn_binned_split %>%
    filter(segVal == 0) %>%
    count(chromosome, cell_name) %>%
    spread(chromosome, n, fill = 0)
  
  # Step 3: Create the heatmap
  cell_chr_dropout_long <- cell_chr_dropout %>%
    gather(key = "chromosome", value = "count", -cell_name)
  
  
  # Now, let's apply this function to heatmap_long
  cell_chr_dropout_long <- cell_chr_dropout_long %>%
    mutate(status = sapply(cell_name, assign_status, cellbased_outliers = cellbased_outliers))
  
  # Ensure the status column in heatmap_long has these four categories
  cell_chr_dropout_long$status <- factor(cell_chr_dropout_long$status, 
                                         levels = c("pass", "rpc_outlier", "alpha_outlier", "replicating","na_alpha_outlier","dmap_outlier","gini_outlier"))
  
  saveRDS(cell_chr_dropout_long,CELL_CHR_DROPOUT_LONG_PATH)

}

ALLSAMPLES=c("23003","24078","24533","23359","23526","24173","24534","24174","24535","23527","24175","24536","23961","24441","24835","23965","24489","24007","24490","24076","24518")
bin_size <- "100"

for (i in 1:length(ALLSAMPLES)){
  SAMPLENAME = ALLSAMPLES[i]
  SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
  ALLCELLS=file.path("Documents/sc_analysis/scAboslute-obj",paste0(SAMPLEONJ,'.rds')) 
  OUTPUT= file.path("Documents/sc_analysis/post-scAbsolute", SAMPLEONJ)
  CELL_POS_DROPOUT_PATH <- paste0(OUTPUT,'/cell_pos_dropout.rds')
  CN_BINNED_PATH <- paste0(OUTPUT,'/cn_binned.rds')
  CN_SEGMENTED_PATH <- paste0(OUTPUT,'/cn_segmented.rds')
  CELL_CHR_DROPOUT_LONG_PATH <- paste0(OUTPUT,'/cell_chr_dropout.rds')
  
  cellbased_outliers = readRDS(paste0(OUTPUT,"/cellbased_outliers.rds"))#summary_outliers.rds
  object<- readRDS(ALLCELLS)
  
  get_summary_of_dropouts(object,cellbased_outliers,CN_BINNED_PATH,CELL_POS_DROPOUT_PATH,CN_SEGMENTED_PATH,CELL_CHR_DROPOUT_LONG_PATH)
  
}
