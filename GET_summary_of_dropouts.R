# For running this, you need to run GET_summary_of_outliers first and have cellbased_outliers object
library("reshape")
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
library(ggplot2)
library(viridis)
library(dendextend)
library(stringr)

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
  ##cn_binned_split <- cn_binned %>%
    ##separate(sample, into = c("SINCEL", "Plate", "Plate_number", "AZNumber"), sep = "_|(?<=[0-9])(?=[A-Z])")
    ###separate(sample, into = c("SINCEL","210","SC", "Plate", "Plate_number", "AZNumber"), sep = "_|(?<=[0-9])(?=[A-Z])")
  
 # cn_binned_split <- cn_binned %>%
#    rowwise() %>%
 #   mutate(parts_count = str_count(sample, "_") + 1) %>%
  #  ungroup() %>%
  #  mutate(
#      sample_split = case_when(
 #       parts_count == 4 ~ list(separate(data.frame(sample = sample), sample, into = c("SINCEL", "Plate", "Plate_number", "AZNumber"), sep = "_|(?<=[0-9])(?=[A-Z])")),
  #      parts_count == 6 ~ list(separate(data.frame(sample = sample), sample, into = c("SINCEL", "210", "SC", "Plate", "Plate_number", "AZNumber"), sep = "_|(?<=[0-9])(?=[A-Z])")),
   #     parts_count == 7 ~ list(separate(data.frame(sample = sample), sample, into = c("SINCEL", "210", "SC", "GEMNEG", "Plate", "Plate_number", "AZNumber"), sep = "_|(?<=[0-9])(?=[A-Z])")),
    #    TRUE ~ list(data.frame(SINCEL = NA, Plate = NA, Plate_number = NA, AZNumber = NA)) # Default case to handle unexpected inputs
#      )
 #   ) %>%
  #  unnest(cols = sample_split, keep_empty = TRUE) %>%
   # select(-parts_count)
  
  # Convert 'sample' column to character to avoid issues with factors
  cn_binned <- cn_binned %>%
    mutate(sample = as.character(sample))
  
  # Extracting the last part (AZNumber) from the sample
  #cn_binned_split <- cn_binned %>%
  #  mutate(AZNumber = sapply(strsplit(sample, "_"), function(x) tail(x, 1)))

  
  cn_binned_split <- cn_binned %>%
    mutate(AZNumber = sapply(gsub("(_[123456789]?)$", "", sample), function(x) {
      last_part <- tail(strsplit(x, "_")[[1]], 1)
      if (last_part == "") "" else last_part
    }))
  
  
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
  
  print(cell_chr_dropout_long$status)
  # Ensure the status column in heatmap_long has these four categories
  cell_chr_dropout_long$status <- factor(cell_chr_dropout_long$status, 
                                         levels = c("pass", "rpc_outlier", "alpha_outlier", "replicating","na_alpha_outlier","dmap_outlier","gini_outlier"))
  
  saveRDS(cell_chr_dropout_long,CELL_CHR_DROPOUT_LONG_PATH)
  
  
  # Assuming your original dataframe is named 'df'
  cell_dropout_status <- cell_chr_dropout_long %>%
    select(-chromosome) %>%  # Remove the chromosome column
    group_by(cell_name, status) %>%  # Group by cell_name and status
    summarise(count = sum(count, na.rm = TRUE)) %>%  # Sum the counts
    ungroup()  # Ungroup to remove the grouping attribute
  
  # Write the new dataframe to a CSV file
  write.csv(cell_dropout_status, CELLS_DROPOUT_STATUS_PATH, row.names = FALSE)

}

ALLSAMPLES=c("23003","24078","24533","23359","23526","24173","24534","24174","24535","23527")
ALLSAMPLES=c("24175","24536","23961","24441","24835","23965","24489","24007","24490","24076","24518")
ALLSAMPLES <- c("24173","24489","24007","23965","24490","24174")
ALLSAMPLES=c("23003","23359","23961","23965","24007")
ALLSAMPLES=c("24078","24518")
ALLSAMPLES=c("24077","24130","24491","23303","24532","23528")
ALLSAMPLES=c("23355","25146","25148","24831","24832","24833")

ALLSAMPLES=c("25146","25147","25149","24831","24832","24833")
ALLSAMPLES=c("24831","24832","24833")

ALLSAMPLES=c("24518","24491")

ALLSAMPLES=c("24911","24912","24913")

ALLSAMPLES=c("24757")

ALLSAMPLES=c("25393","25394","25394_B4","25394_C9","25394_D5","25394_D7","25394_D11","25394_G3")
categories <- c("PEO1 parent","PEO1 children","PEO1 B4","PEO1 C9","PEO1 D5","PEO1 B7","PEO1 D11","PEO1 G3")

bin_size <- "100"

for (i in 1:length(ALLSAMPLES)){
  SAMPLENAME = ALLSAMPLES[i]
  SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
  ALLCELLS=file.path("../../Volumes/Fl/sc_analysis/scAboslute-obj",paste0(SAMPLEONJ,'.rds')) 
  OUTPUT= file.path("../../Volumes/Fl/sc_analysis/post-scAbsolute", SAMPLEONJ)
  CELL_POS_DROPOUT_PATH <- paste0(OUTPUT,'/cell_pos_dropout.rds')
  CN_BINNED_PATH <- paste0(OUTPUT,'/cn_binned.rds')
  CN_SEGMENTED_PATH <- paste0(OUTPUT,'/cn_segmented.rds')
  CELL_CHR_DROPOUT_LONG_PATH <- paste0(OUTPUT,'/cell_chr_dropout.rds')
  CELLS_DROPOUT_STATUS_PATH <- paste0(OUTPUT,"/cells_dropout_status.csv")
  cellbased_outliers = readRDS(paste0(OUTPUT,"/cellbased_outliers.rds"))#summary_outliers.rds
  object<- readRDS(ALLCELLS)
  
  
  get_summary_of_dropouts(object,cellbased_outliers,CN_BINNED_PATH,CELL_POS_DROPOUT_PATH,CN_SEGMENTED_PATH,CELL_CHR_DROPOUT_LONG_PATH)
  
}

