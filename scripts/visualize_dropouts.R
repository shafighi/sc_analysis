# For running this, you need to run GET_summary_of_dropouts first and have cell_pos_dropout,
# cn_binned, cn_segmented, and cell_chr_dropout

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

ALLSAMPLES=c("23003","24078","24533","23359","23526","24173","24534","24174","24535","23527","24175","24536","23961","24441","24835","23965","24489","24007","24490","24076","24518")
ALLSAMPLES=c("23003","23359","23961","23965","24007")
ALLSAMPLES=c("24078","24518")#,"24835"
ALLSAMPLES=c("24078")
ALLSAMPLES=c("24077","24130","24491","23303","24532","23528")

ALLSAMPLES <- c("23355","25146","24831","24832","24833","25146")
ALLSAMPLES=c("25146","25147","25148","25149","24831","24832","24833")
ALLSAMPLES=c("24831","24832","24833")
ALLSAMPLES=c("24518")
ALLSAMPLES=c("24911","24912","24913")
ALLSAMPLES=c("24757")

ALLSAMPLES=c("25393","25394","25394_B4","25394_C9","25394_D5","25394_D7","25394_D11","25394_G3")

ALLSAMPLES=c("25393")
bin_size <- "100"

for (i in 1:length(ALLSAMPLES)){
  SAMPLENAME = ALLSAMPLES[i]
  SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
  ALLCELLS=file.path("../../Volumes/Fl/sc_analysis/scAboslute-obj",paste0(SAMPLEONJ,'.rds')) 
  OUTPUT= file.path("../../Volumes/Fl/sc_analysis/post-scAbsolute", SAMPLEONJ)
  if (file.exists(paste0(OUTPUT,'/figures'))) {
    cat("Directory is already there!", paste0(OUTPUT,'/figures'), "\n")
  } else {
    dir.create(paste0(OUTPUT,'/figures'))
    cat("Directory is created:", paste0(OUTPUT,'/figures'), "\n")
  }
  DROPOUT_PATH_SEGMENT <- paste0(OUTPUT,'/figures/dropout_segment.pdf')
  DROPOUT_PATH_BIN <- paste0(OUTPUT,'/figures/dropout_bin.pdf')
  CN_BINNED_HEATMAP_DROPOUT <- paste0(OUTPUT,'/figures/dropout_heatmap.png')
  CN_BINNED_HEATMAP_DROPOUT2 <- paste0(OUTPUT,'/figures/dropout_heatmap2.png')
  CN_BINNED_HEATMAP_DROPOUT_PASSED <- paste0(OUTPUT,'/figures/dropout_heatmap2_passed.png')
  
  CELL_POS_DROPOUT_PATH <- paste0(OUTPUT,'/cell_pos_dropout.rds')
  CN_BINNED_PATH <- paste0(OUTPUT,'/cn_binned.rds')
  CN_SEGMENTED_PATH <- paste0(OUTPUT,'/cn_segmented.rds')
  CELL_CHR_DROPOUT_LONG_PATH <- paste0(OUTPUT,'/cell_chr_dropout.rds')
  
  object<- readRDS(ALLCELLS)
  
  plot_dropouts(object,CELL_POS_DROPOUT_PATH,CN_BINNED_PATH,CN_SEGMENTED_PATH,CELL_CHR_DROPOUT_LONG_PATH,DROPOUT_PATH_SEGMENT,DROPOUT_PATH_BIN,CN_BINNED_HEATMAP_DROPOUT,CN_BINNED_HEATMAP_DROPOUT2)
}


plot_dropouts <- function(object,CELL_POS_DROPOUT_PATH,CN_BINNED_PATH,CN_SEGMENTED_PATH,CELL_CHR_DROPOUT_LONG_PATH,DROPOUT_PATH_SEGMENT,DROPOUT_PATH_BIN,CN_BINNED_HEATMAP_DROPOUT,CN_BINNED_HEATMAP_DROPOUT2){
  cn_segmented <- readRDS(CN_SEGMENTED_PATH)
  cn_binned <- readRDS(CN_BINNED_PATH)
  
  
  dropout <- as.data.frame(table(cn_segmented[cn_segmented$segVal==0,]$chromosome) )
  pdf(DROPOUT_PATH_SEGMENT, width = 6, height = 3)
  barplot(dropout$Freq, names.arg = dropout$Var1, xlab = "", ylab = "Dropout Frequency (segments)", col = "skyblue",
          main = SAMPLENAME, cex.axis = 0.4, cex.lab = 0.4, cex.names = 0.4)
  dev.off()
  
  
  dropout <- as.data.frame(table(cn_binned[cn_binned$segVal==0,]$chromosome) )
  pdf(DROPOUT_PATH_BIN, width = 6, height = 3)
  barplot(dropout$Freq, names.arg = dropout$Var1, xlab = "", ylab = "Dropout Frequency (bins)", col = "skyblue",
          main = SAMPLENAME, cex.axis = 0.4, cex.lab = 0.4, cex.names = 0.4)
  dev.off()
  
  
  cell_pos_dropout <- readRDS(CELL_POS_DROPOUT_PATH)
  
  # Create the heatmap
  heatmap_plot <- ggplot(cell_pos_dropout, aes(x = Number, y = AZ, fill = total_dropout)) +
    geom_tile(color = "white") +
    geom_text(aes(label = total_dropout), color = "black", size = 3) + # Add text labels
    #scale_fill_gradient2(low = "lightgreen", mid = "lightyellow", high = "plum", midpoint = max(cell_pos_dropout$total_dropout)/2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = max(cell_pos_dropout$total_dropout)/2) +
      labs(title = paste0("Heatmap of total number of bins with 0 reads for each cell: ",SAMPLENAME), x = "Number", y = "AZ") +
    theme_minimal() +
    theme(
      text = element_text(size = 12), # Adjust the base font size
      axis.title = element_text(size = 12), # Adjust axis title font size
      axis.text = element_text(size = 12), # Adjust axis text font size
      plot.title = element_text(size = 12) # Adjust plot title font size
    )
  
  # Print the heatmap to the console (optional)
  print(heatmap_plot)
  
  # Save the heatmap to a file
  ggsave(CN_BINNED_HEATMAP_DROPOUT, plot = heatmap_plot, width = 10, height = 6, dpi = 300)

  
  
  # Create the heatmap
  heatmap_plot_passed <- ggplot(cell_pos_dropout, aes(x = Number, y = AZ, fill = total_dropout)) +
    geom_tile(color = "white") +
    geom_text(aes(label = total_dropout), color = "black", size = 3) + # Add text labels
    #scale_fill_gradient2(low = "lightgreen", mid = "lightyellow", high = "plum", midpoint = max(cell_pos_dropout$total_dropout)/2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = max(cell_pos_dropout$total_dropout)/2) +
    labs(title = paste0("Heatmap of total number of bins with 0 reads for each cell: ",SAMPLENAME), x = "Number", y = "AZ") +
    theme_minimal() +
    theme(
      text = element_text(size = 12), # Adjust the base font size
      axis.title = element_text(size = 12), # Adjust axis title font size
      axis.text = element_text(size = 12), # Adjust axis text font size
      plot.title = element_text(size = 12) # Adjust plot title font size
    )
  
  # Print the heatmap to the console (optional)
  print(heatmap_plot_passed)
  
  # Save the heatmap to a file
  ggsave(CN_BINNED_HEATMAP_DROPOUT_PASSED, plot = heatmap_plot, width = 10, height = 6, dpi = 300)  
  
  cell_chr_dropout <- readRDS(CELL_CHR_DROPOUT_LONG_PATH)
  
  # Create the heatmap with clustered cells within each facet and increased space between facets
  heatmap_plot_2 <- ggplot(cell_chr_dropout, aes(x = chromosome, y = cell_name, fill = count)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = max(cell_chr_dropout$count)/2) +
    facet_grid(status ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text.y = element_text(angle = 0),
      strip.background = element_rect(fill = "lightgrey", colour = "black", size = 1),
      panel.spacing = unit(2, "cm"),  # Increase space between facets
      strip.text = element_text(size = 12, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    labs(title = "Heatmap of Counts by Cell Name and Chromosome",
         x = "Chromosome",
         y = "Cell Name",
         fill = "Count")
  
  # Save the plot
  ggsave(CN_BINNED_HEATMAP_DROPOUT2, plot = heatmap_plot_2, width = 8, height = 40, dpi = 300)
  
}
