
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

BASE="OneDrive - CRUK Cambridge Institute/scCNV_analysis"
BASEOUT="OneDrive - CRUK Cambridge Institute/scCNV_analysis/post-scAbsolute"
SAMPLENAME = "24532"
bin_size <- "100"
SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
ALLCELLS=file.path(BASE,"scAboslute-obj",paste0(SAMPLEONJ,'.rds')) 
OUTPUT= file.path(BASEOUT, SAMPLEONJ)

if (file.exists(OUTPUT)) {
  cat("Directory is already there!", OUTPUT, "\n")
} else {
  dir.create(OUTPUT)
  cat("Directory is created:", OUTPUT, "\n")
}

DROPOUT_PATH_SEGMENT <- paste0(OUTPUT,'/dropout_segment.pdf')
DROPOUT_PATH_BIN <- paste0(OUTPUT,'/dropout_bin.pdf')
CN_BINNED_PATH <- paste0(OUTPUT,'/cn_binned.rds')
CN_BINNED_HEATMAP_DROPOUT <- paste0(OUTPUT,'/dropout_heatmap.png')
CN_BINNED_HEATMAP_DROPOUT2 <- paste0(OUTPUT,'/dropout_heatmap2.png')

CELL_POS_DROPOUT_PATH <- paste0(OUTPUT,'/cell_pos_dropout.rds')
CN_BINNED_PATH <- paste0(OUTPUT,'/cn_binned.rds')
CN_SEGMENTED_PATH <- paste0(OUTPUT,'/cn_segmented.rds')
CELL_CHR_DROPOUT_LONG_PATH <- paste0(OUTPUT,'/cell_chr_dropout.rds')
object<- readRDS(ALLCELLS)
cn <- object@assayData$copynumber
summary_df_outliers <- readRDS(paste0(OUTPUT,"/summary_outliers.rds"))


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
  geom_text(aes(label = total_dropout), color = "black", size = 2) + # Add text labels
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = max(cell_pos_dropout$total_dropout)/2) +
  labs(title = paste0("Heatmap of total number of bins with 0 reads for each cell: ",SAMPLENAME), x = "Number", y = "AZ") +
  theme_minimal() +
  theme(
    text = element_text(size = 7), # Adjust the base font size
    axis.title = element_text(size = 7), # Adjust axis title font size
    axis.text = element_text(size = 7), # Adjust axis text font size
    plot.title = element_text(size = 7) # Adjust plot title font size
  )

# Print the heatmap to the console (optional)
print(heatmap_plot)

# Save the heatmap to a file
ggsave(CN_BINNED_HEATMAP_DROPOUT, plot = heatmap_plot, width = 10, height = 6, dpi = 300)


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



