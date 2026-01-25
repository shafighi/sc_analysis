
library("reshape")
library(dplyr)
library(ComplexHeatmap)
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/scUnique-main/R/visualize.R")


CORES=7
RMNORM=TRUE
PREPATH="OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/data/metadata/"

ALLSAMPLES <- c("23962","23963","23964","24831","24832","24833","25146","25147","25148","25149","24911","24912","24913","24914","24915","24806","24807","24808")

ALLSAMPLES <- c("23964","24831","24832","24833","25146","25147","25149","24911","24912","24913","24915","24806","24807","24808")
samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1")

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24831","24832","24833","25146","25147","25148","25149")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA","FFPE")  # Adjust this to your preferred order



bin_size <- "100"



ht_list = list()

  SAMPLENAME = "24077"
  SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
  ALLCELLS = paste0("Documents/sc_analysis/scUnique-obj/SLX-",SAMPLENAME,"/jointSegmentationObject.rds")
  OUTPUT = paste0("Documents/sc_analysis/scUnique-obj/SLX-",SAMPLENAME)
  if (file.exists(OUTPUT)) {
    cat("Directory is already there!", OUTPUT, "\n")
  } else {
    dir.create(OUTPUT)
    cat("Directory is created:", OUTPUT, "\n")
  }
  NON_OUTLIER_OBJ_PATH = paste0(OUTPUT,"/",SAMPLENAME,"_non_outlier_scUnique.rds")
  jointSegmentationObject = readRDS(ALLCELLS)
  object <- jointSegmentationObject$newObject 
  
  HEATMAP_PATH <- paste0(OUTPUT,"/",SAMPLENAME,"_scUnique_non_outlier_copynumber.pdf")
  

  #HEATMAP_PATH
  plotCopynumberHeatmap(object,file =HEATMAP_PATH ,cluster_rows = FALSE, row_split=NULL,
                             cutoff=10, show_unobserved_states=TRUE,
                             har=NULL, useCopynumber=TRUE,
                             show_cell_names=FALSE, abbreviate_cell_names=TRUE, show_chromosome_names=FALSE,
                             use_cell_names="", fontsize_row=9, fontsize_col=9, fontsize_chr=12,
                             fontsize_leg_title=18, fontsize_leg_label=14, 
                             row_title_gp=13, column_title_gp=13,
                             column_title="", bottom_annotation=NULL,
                             show_heatmap_legend=TRUE, scale_copynumber=1.0,
                             raster_device="tiff",
                             colorMap="MSK",row_title=paste0("PEO4: ",SAMPLENAME))
  
