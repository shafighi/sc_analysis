library(caret)
library(randomForest)
library(e1071)
library(class)
library(ggplot2)
library("flexmix")
library(NMF)
library(tictoc)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)


CORES=15
SEED=1990
NMFALG="brunet"
ITER=15


OUTPUT <- paste0('Documents/sc_analysis/all_samples_23july2024/')
allMats <- readRDS(paste0(OUTPUT,"sum_of_posterior.rds"))
#,"changepoint1", "changepoint2", "changepoint3", "changepoint4"
selected_cols <- c("segsize1","segsize2","segsize3","segsize4","segsize5","segsize6","segsize7","segsize8","bpchrarm1","bpchrarm2","bpchrarm3","uechr1","uechr2", "uechr3","changepoint1", "changepoint2" ,"changepoint3" , "uechr4" , "uechr5" )      
#sample_names <- c("24007","23526","24489","24490","24077","23965") 
sample_names <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24173","24489","24007","23965","24490","24174")


for (RANK in c(3,4,5)){
  
  nmfObject = NMF::nmf(allMats[,selected_cols], rank = RANK, seed = SEED, nrun = ITER, method = NMFALG, .opt = paste0("p", CORES) )
  BASE=paste0(OUTPUT,"SCSIG_",RANK,"/")
    
  if (file.exists(BASE)) {
    cat("Directory is already there!", BASE, "\n")
  } else {
    dir.create(BASE, recursive = TRUE)
    cat("Directory is created:", BASE, "\n")
  } 
  saveRDS(object = nmfObject, file = paste0(BASE, "nmfObject.rds"))
  
  # Normalise both and save as txt and rds
  if(nrow(allMats) > ncol(allMats)) {
    ## Sample by component input
    sigs = nmfObject@fit@H
    exp = nmfObject@fit@W
  } else {
    ## Component by sample input
    sigs = t( nmfObject@fit@W )
    exp = t( nmfObject@fit@H )
  }
  
  rownames(sigs) = paste0("S", 1:nrow(sigs))
  write.table(sigs, file = paste0(BASE, "Signatures.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  sigs_normalised = t( apply(sigs, 1, function(x) x/sum(x)) )
  write.table(sigs_normalised, file = paste0(BASE, "Signatures_normalised.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  colnames(exp) = paste0("S", 1:ncol(exp))
  exp_normalised = t( apply(exp, 1, function(x) x/sum(x)) )
  write.table(exp_normalised, file = paste0(BASE, "Exposures_normalised.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

}


pdf(file = paste0(OUTPUT, "Signatures.pdf"), width = 6, height = 2.5)
for (RANK in c(3,4,5)){
  
  BASE=paste0(OUTPUT,"SCSIG_",RANK,"/")
  sigs_normalised <- read.table( file = paste0( BASE,"Signatures_normalised.txt"), sep = "\t", header = TRUE, row.names = 1)
  exp_normalised <- read.table(file = paste0( BASE,"Exposures_normalised.txt"), sep = "\t", header = TRUE, row.names = 1)
  
  print(Heatmap(sigs_normalised, column_title = "Normalised signatures", cluster_rows = FALSE, cluster_columns = FALSE))
  sample_sincel <- sub("^(SINCEL[-_]\\d+)_((SC_GEMNEG_Plate|SC_Plate|Plate)_\\d+_.+)$", "\\1", rownames(exp_normalised))

  sample_sincel <- gsub("SINCEL-194", "SINCEL-194 (PEO1)", sample_sincel)
  sample_sincel <- gsub("SINCEL_210", "SINCEL_210 (PEO4)", sample_sincel)
  sample_sincel <- gsub("SINCEL-245", "SINCEL-245 (UWB1.289 BRCA)", sample_sincel)
  sample_sincel <- gsub("SINCEL-246", "SINCEL-246 (UWB1.289)", sample_sincel)
  sample_sincel <- gsub("SINCEL-219", "SINCEL-219 (PEO4)", sample_sincel)
  sample_sincel <- gsub("SINCEL-211", "SINCEL-211 (PEO1)", sample_sincel)
  
  color_palette <- brewer.pal(length( unique(sample_sincel) ), "Set1")
  names(color_palette) <- unique(sample_sincel)  
  column_ha = HeatmapAnnotation(Samples = sample_sincel,col = list(Samples = color_palette))
}
dev.off()
  


pdf(file = paste0(OUTPUT, "Normalised_Exp.pdf"), width =4, height =6)
for (RANK in c(3,4,5)){
  
  BASE=paste0(OUTPUT,"SCSIG_",RANK,"/")
  sigs_normalised <- read.table( file = paste0( BASE,"Signatures_normalised.txt"), sep = "\t", header = TRUE, row.names = 1)
  exp_normalised <- read.table(file = paste0( BASE,"Exposures_normalised.txt"), sep = "\t", header = TRUE, row.names = 1)

  sample_sincel <- sub("^(SINCEL[-_]\\d+)_((SC_GEMNEG_Plate|SC_Plate|Plate)_\\d+_.+)$", "\\1", rownames(exp_normalised))
  
  sample_sincel <- gsub("SINCEL-194", "SINCEL-194 (PEO1)", sample_sincel)
  sample_sincel <- gsub("SINCEL_210", "SINCEL_210 (PEO4)", sample_sincel)
  sample_sincel <- gsub("SINCEL-245", "SINCEL-245 (UWB1.289 BRCA)", sample_sincel)
  sample_sincel <- gsub("SINCEL-246", "SINCEL-246 (UWB1.289)", sample_sincel)
  sample_sincel <- gsub("SINCEL-219", "SINCEL-219 (PEO4)", sample_sincel)
  sample_sincel <- gsub("SINCEL-211", "SINCEL-211 (PEO1)", sample_sincel)
  
  color_palette <- brewer.pal(length( unique(sample_sincel) ), "Set1")
  names(color_palette) <- unique(sample_sincel)  
  column_ha = rowAnnotation(Samples = sample_sincel,col = list(Samples = color_palette))
  
  print(Heatmap(
    matrix = (exp_normalised),
    name = "Exp",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_dend_side = "left",
    left_annotation = column_ha
  ))
}
dev.off()

pdf(file = paste0(OUTPUT, "Normalised_Exp_Clustered.pdf"), width = 12, height = 8)
for (RANK in c(3,4,5,6,7,8,9,10,11,12,13,14,15)){
  
  BASE=paste0(OUTPUT,"SCSIG_",RANK,"/")
  sigs_normalised <- read.table( file = paste0( BASE,"Signatures_normalised.txt"), sep = "\t", header = TRUE, row.names = 1)
  exp_normalised <- read.table(file = paste0( BASE,"Exposures_normalised.txt"), sep = "\t", header = TRUE, row.names = 1)
  
  sample_sincel <- sub("^(SINCEL[-_]\\d+)_((SC_GEMNEG_Plate|SC_Plate|Plate)_\\d+_.+)$", "\\1", rownames(exp_normalised))
  sample_sincel <- gsub("SINCEL-194", "SINCEL-194 (PEO1)", sample_sincel)
  sample_sincel <- gsub("SINCEL_210", "SINCEL_210 (PEO4)", sample_sincel)
  sample_sincel <- gsub("SINCEL-245", "SINCEL-245 (UWB1.289 BRCA)", sample_sincel)
  sample_sincel <- gsub("SINCEL-246", "SINCEL-246 (UWB1.289)", sample_sincel)
  sample_sincel <- gsub("SINCEL-219", "SINCEL-219 (PEO4)", sample_sincel)
  sample_sincel <- gsub("SINCEL-211", "SINCEL-211 (PEO1)", sample_sincel)
  
  
  color_palette <- brewer.pal(length( unique(sample_sincel) ), "Set1")
  names(color_palette) <- unique(sample_sincel)  
  column_ha = HeatmapAnnotation(Samples = sample_sincel,col = list(Samples = color_palette))
  
  print(Heatmap(
    matrix = t(exp_normalised),
    name = "Exp",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_dend_side = "left",
    top_annotation = column_ha
  ))
}
dev.off()





pdf(file = paste0(OUTPUT, "Mixture_components.pdf"), width = 6, height = 6)

sample_sincel <- sub("^(SINCEL[-_]\\d+)_((SC_GEMNEG_Plate|SC_Plate|Plate)_\\d+_.+)$", "\\1", rownames(allMats[,selected_cols]))

sample_sincel <- gsub("SINCEL-194", "SINCEL-194 (PEO1)", sample_sincel)
sample_sincel <- gsub("SINCEL_210", "SINCEL_210 (PEO4)", sample_sincel)
sample_sincel <- gsub("SINCEL-245", "SINCEL-245 (UWB1.289 BRCA)", sample_sincel)
sample_sincel <- gsub("SINCEL-246", "SINCEL-246 (UWB1.289)", sample_sincel)
sample_sincel <- gsub("SINCEL-219", "SINCEL-219 (PEO4)", sample_sincel)
sample_sincel <- gsub("SINCEL-211", "SINCEL-211 (PEO1)", sample_sincel)

color_palette <- brewer.pal(length( unique(sample_sincel) ), "Set1")
names(color_palette) <- unique(sample_sincel)  
column_ha = rowAnnotation(Samples = sample_sincel,col = list(Samples = color_palette))

print(Heatmap(
  matrix = allMats[,selected_cols],
  name = "Mixture components",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_dend_side = "left",
  left_annotation = column_ha
))
dev.off()

