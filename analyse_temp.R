library("reshape")
library(dplyr)
library(ComplexHeatmap)
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/scUnique-main/R/visualize.R")

OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[s],'/')
HEATMAP_PATH <- paste0(OUTPUT,"scUnique_copynumber.pdf")

bin_size=100
ALLSAMPLES=c("24518","24491")
categories <- c("PEO STOP","PEO MISSENSE")
SAMPLENAME = ALLSAMPLES[2]
SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
OUTPUT= file.path("Documents/sc_analysis/post-scAbsolute", SAMPLEONJ)
CELLS_DROPOUT_STATUS_PATH <- paste0(OUTPUT,"/cells_dropout_status.csv")
cell_dropout_status <- read.csv(CELLS_DROPOUT_STATUS_PATH)
passed = cell_dropout_status[cell_dropout_status$status=="pass",]
passed[passed$count>400,]

cellbased_outliers = readRDS(paste0(OUTPUT,"/cellbased_outliers.rds"))#summary_outliers.rds

plotCopynumberHeatmap(new_object,file = HEATMAP_PATH,show_chromosome_names=TRUE)
