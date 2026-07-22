library("reshape")
library(dplyr)
library(ComplexHeatmap)
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/scUnique-main/R/visualize.R")

CORES=7
RMNORM=TRUE
PREPATH="OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/data/metadata/"
samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532")
samples <- c("24173","24489","24007","23965","24490","24174")
samples <- c("24518", "24491")
samples <- c("23355")
samples <- c("24831", "24832" ,"24833")
samples <- c("25147","25148")
samples <- c("25149")
samples <- c("24911","24912","24913")
samples <- c("23964","24831","24832","24833","25146","25147","25149","24911","24912","24913","24915","24806","24807","24808","24809")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")

samples <- c("23962","23963","23964","24831","24832","24833","25146","25147","25148","25149","24911","24912","24913","24914","24915","24806","24807","24808","24809")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2_null","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289")


bin_size <- "100"


for (s in 1:length(samples)){
  print(paste0("SAMPLE NAME: ", samples[s]))
  print("------------------ SETTING FILE NAMES ------------------")
  OUTPUT <- paste0('../../Volumes/Fl/sc_analysis/scUnique-obj/scAbsolute/FFPE_heatmaps/')
  OBJ_PATH <- paste0('../../Volumes/Fl/sc_analysis/post-scAbsolute/SLX-',samples[s],'_100/',samples[s],'_non_outlier.rds')
  if (file.exists(OUTPUT)) {
    cat("Directory is already there!", OUTPUT, "\n")
  } else {
    dir.create(OUTPUT, recursive = TRUE)
    cat("Directory is created:", OUTPUT, "\n")
  }  
  
  
  HEATMAP_PATH <- paste0(OUTPUT,categories[s],"_",samples[s],"_scAbsolute_copynumber_non_outlier.pdf")
  object<- readRDS(OBJ_PATH)
  
  plotCopynumberHeatmap(object,file =HEATMAP_PATH ,cluster_rows = FALSE, row_split=NULL,
                        cutoff=10, show_unobserved_states=TRUE,
                        har=NULL, useCopynumber=TRUE,
                        show_cell_names=FALSE, abbreviate_cell_names=TRUE, show_chromosome_names=FALSE,
                        use_cell_names=NULL, fontsize_row=9, fontsize_col=9, fontsize_chr=12,
                        fontsize_leg_title=18, fontsize_leg_label=14, 
                        row_title_gp=13, column_title_gp=13,
                        column_title="", bottom_annotation=NULL,
                        show_heatmap_legend=TRUE, scale_copynumber=1.0,
                        raster_device="tiff",
                        colorMap="MSK",row_title=paste0(categories[s],": ",samples[s]))
  
}
