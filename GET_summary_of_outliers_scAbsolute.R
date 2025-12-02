
library(gridExtra)
library(pheatmap) 
library(Biobase)
library(gridExtra)
library(ggplot2)
library(dplyr)
#source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/sc_visualization.R")
source("R/core.R")
source("R/visualization_helpers.R")

current_directory <- getwd()
print(paste("Current Working Directory:", current_directory))



ALLSAMPLES=c("23003","24078","24533","23359","23526","24173","24534","24174","24535","23527","24175","24536","23961","24441","24835","23965","24489","24007","24490","24076","24518")

ALLSAMPLES <- c("24173","24489","24007","23965","24490","24174")
ALLSAMPLES=c("24077","24130","24491","23303","24532","23528")
ALLSAMPLES=c("23355","25146","25148","24831","24832","24833")

ALLSAMPLES=c("25146","25147","25149","24831","24832","24833")
ALLSAMPLES=c("24518","24491")
categories <- c("PEO STOP","PEO MISSENSE")

ALLSAMPLES=c("24911","24912","24913")
categories <- c("FFPE","FFPE","FFPE")

ALLSAMPLES <- c("23964","24831","24832","24833","25146","25147","25149","24911","24912","24913","24915","24806","24807","24808")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")

ALLSAMPLES <- c("23962","23963","23964","24831","24832","24833","25146","25147","25148","25149","24911","24912","24913","24914","24915","24806","24807","24808","24809")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")

ALLSAMPLES=c("24809")
categories <- c("FFPE")

ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24174","24173","24489","24007","23965","24490","24532")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","PEO1")

ALLSAMPLES=c("24518","24491")
categories <- c("PEO STOP","PEO MISSENSE")
ALLSAMPLES=c("24757")
categories <- c("Organoid")

ALLSAMPLES=c("25393","25394","25394_B4","25394_C9","25394_D5","25394_D7","25394_D11","25394_G3")
categories <- c("PEO1 parent","PEO1 children","PEO1 B4","PEO1 C9","PEO1 D5","PEO1 B7","PEO1 D11","PEO1 G3")
ALLSAMPLES=c("25393")
categories <- c("PEO1 parent")

bin_size <- "100"

for (i in 1:length(ALLSAMPLES)){
  
  SAMPLENAME = ALLSAMPLES[i]
  SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
  ALLCELLS=file.path("../../Volumes/Fl/sc_analysis/scAboslute-obj",paste0(SAMPLEONJ,'.rds')) 
  OUTPUT= file.path("../../Volumes/Fl/sc_analysis/post-scAbsolute", SAMPLEONJ)
  if (file.exists(OUTPUT)) {
    cat("Directory is already there!", OUTPUT, "\n")
  } else {
    dir.create(OUTPUT)
    cat("Directory is created:", OUTPUT, "\n")
  }
  OUTLIERS_CELLBASE_PATH = paste0(OUTPUT,"/cellbased_outliers.rds")#summary_outliers.rds
  NON_OUTLIER_OBJ_PATH = paste0(OUTPUT,"/",SAMPLENAME,"_non_outlier.rds")
  OUTLIER_SUMMARY_PATH = paste0(OUTPUT,"/outlier_summary.rds")#summary_df.rds
  NORMALS_PATH = paste0(OUTPUT,"/normals.rds")#summary_df.rds
  object = readRDS(ALLCELLS)
  SLX_value <- paste0("SLX-",SAMPLENAME)
  UID_value <- categories[i]
  get_summary_of_outliers(object,OUTLIERS_CELLBASE_PATH,NON_OUTLIER_OBJ_PATH,OUTLIER_SUMMARY_PATH,SLX_value,UID_value,NORMALS_PATH)
  

  normal_rows = rownames(readRDS(NORMALS_PATH))
  non_outlier_object <- readRDS(NON_OUTLIER_OBJ_PATH)
  
  #Biobase::pData(object)$IsNormal <- seq_len(nrow(df)) %in% normal_rows

  
  plotCopynumberHeatmap(non_outlier_object,file =paste0(OUTPUT,"/heatmap_clustered.pdf") ,cluster_rows = "ploidy", row_split=NULL,
                        cutoff=10, show_unobserved_states=TRUE,
                        har=NULL, useCopynumber=TRUE,
                        show_cell_names=FALSE, abbreviate_cell_names=TRUE, show_chromosome_names=FALSE,
                        use_cell_names=NULL, fontsize_row=9, fontsize_col=9, fontsize_chr=12,
                        fontsize_leg_title=18, fontsize_leg_label=14, 
                        row_title_gp=13, column_title_gp=13,
                        column_title="", bottom_annotation=NULL,
                        show_heatmap_legend=TRUE, scale_copynumber=1.0,
                        raster_device="tiff",
                        colorMap="MSK",row_title=SAMPLENAME)
}

