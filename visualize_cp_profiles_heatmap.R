
library("reshape")
library(dplyr)
library(ComplexHeatmap)
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/scUnique-main/R/visualize.R")


CORES=7
RMNORM=TRUE
PREPATH="OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/data/metadata/"

samples <- c("23962","23963","23964","24831","24832","24833","25146","25147","25148","25149","24911","24912","24913","24914","24915","24806","24807","24808")

samples <- c("23964","24831","24832","24833","25146","25147","25149","24911","24912","24913","24915","24806","24807","24808")

bin_size <- "100"


ht_list = list()
for (s in 1:length(samples)){
  print(paste0("SAMPLE NAME: ", samples[s]))
  print("------------------ SETTING FILE NAMES ------------------")
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/scAbsolute/SLX-',samples[s],'/')
  OBJ_PATH <- paste0('Documents/sc_analysis/scAboslute-obj/SLX-',samples[s],"_100.rds")
  if (file.exists(OUTPUT)) {
    cat("Directory is already there!", OUTPUT, "\n")
  } else {
    dir.create(OUTPUT, recursive = TRUE)
    cat("Directory is created:", OUTPUT, "\n")
  }  
  
  HEATMAP_PATH <- paste0(OUTPUT,samples[s],"_scAbsolute_copynumber.pdf")

  
  
  print("------------------ LOADING OBJECTS ------------------")
  object<- readRDS(OBJ_PATH)

#HEATMAP_PATH
  ht = plotCopynumberHeatmap(object,file =NULL ,cluster_rows = FALSE, row_split=NULL,
                      cutoff=10, show_unobserved_states=TRUE,
                      har=NULL, useCopynumber=TRUE,
                      show_cell_names=FALSE, abbreviate_cell_names=TRUE, show_chromosome_names=FALSE,
                      use_cell_names="", fontsize_row=9, fontsize_col=9, fontsize_chr=12,
                      fontsize_leg_title=18, fontsize_leg_label=14, 
                      row_title_gp=13, column_title_gp=13,
                      column_title="", bottom_annotation=NULL,
                      show_heatmap_legend=TRUE, scale_copynumber=1.0,
                      raster_device="tiff",
                      colorMap="MSK",row_title=paste0("FFPE: ",samples[s]))

  ht_list[[s]] = ht
  
}

final_ht = Reduce(`%v%`, ht_list)

pdf('Documents/sc_analysis/all.pdf', width = 11, height = 50)
draw(final_ht)
dev.off()

