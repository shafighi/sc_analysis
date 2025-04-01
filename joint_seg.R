library(dplyr)
library(QDNAseq, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(ggbeeswarm, quietly = TRUE)
library("ggpubr", quietly = TRUE)
library(robustbase, quietly = TRUE)
library(Biobase, quietly = TRUE)

BASEU="/home/schnei01/scUnique"
source(file.path(BASEU,"R/core.R"))
BASEDIR="/home/schnei01/scAbsolute"
#ALLSAMPLES="out_23525_100.rds"
source(file.path(BASEU, "sc_visualization.R"))

# Check if commandArgs is not empty
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  # Get the value of ALLSAMPLES from command line
  ALLSAMPLES <- commandArgs(trailingOnly = TRUE)[1]
  OUTPUT <- commandArgs(trailingOnly = TRUE)[2]
  HMMPATH <- commandArgs(trailingOnly = TRUE)[3]
  # Print the value for verification (optional)
  cat("ALLSAMPLES:", ALLSAMPLES, "\n")

  # Rest of your script using ALLSAMPLES
  # ...
} else {

  cat("Error: Please provide the value for ALLSAMPLES and OUTPUT as a command-line argument.\n")
  # Optionally, print instructions on how to run the script
}

SAMPLENAME = "out"

# Output folder path
#output_folder <- file.path(BASEU, OUTPUT)

# Check if the folder exists, and create it if not
if (!dir.exists(OUTPUT)) {
  dir.create(OUTPUT, recursive = TRUE)
}

source(file.path(BASEDIR, "R/core.R"))
source(file.path(BASEDIR, "R/visualization.R"))


scabsolute = readRDS(ALLSAMPLES)
df = predict_replicating(Biobase::pData(scabsolute) %>%
                           tidyr::separate(name, sep="_", into=c("UID", "SLX", "cellid", "celltag"), remove=FALSE)
                         ,cutoff_value=1,iqr_value = 1.5)

all_cellids <- unique(df$cellid)
na_alpha_cells <- df[!complete.cases(df$hmm.alpha),]
df <- df[complete.cases(df$hmm.alpha), ]
 
#scabsolute <- scabsolute[, complete.cases(df$hmm.alpha)]
rep_cells = df[df$replicating==TRUE,]$name
condition_to_remove <- pData(scabsolute)$name %in% rep_cells
non_rep_object <- scabsolute[, !condition_to_remove]


mapd_cutoff <- 1.5
mapd_density_control <- TRUE
gini_norm_cutoff <- 1.5
alpha_cutoff <- 2
alpha_hard_cutoff <- 0.05
gini_density_control <- TRUE
rpc_cutoff <- 15
densityCutoff <- 0.1
cutoff_percentile <-0.05

#filtered_rpc <- df %>% dplyr::filter(!replicating, rpc >= rpc_cutoff)
#filtered_mapd <- qc_mapd(filtered_rpc, cutoff_value = mapd_cutoff, densityControl=mapd_density_control)
#filtered_gini <- qc_gini_norm(filtered_mapd, cutoff_value = gini_norm_cutoff, densityControl = gini_density_control)
#iq = qc_alpha(filtered_gini, cutoff_value=alpha_cutoff, hard_cutoff = alpha_hard_cutoff)

filtered_rpc <- df[df$replicating==FALSE,] %>% dplyr::filter(!replicating, rpc >= rpc_cutoff)    
filtered_mapd <- qc_mapd(filtered_rpc,cutoff_percentile=cutoff_percentile, cutoff_value = mapd_cutoff, symmetric=FALSE, densityControl=mapd_density_control, densityCutoff=densityCutoff) 
filtered_mapd$dmapd.outlier <- filtered_mapd$dmapd.outlier[,"outlier"] 
filtered_gini <- qc_gini_norm(filtered_mapd, cutoff_value = gini_norm_cutoff, densityControl = gini_density_control, densityCutoff=densityCutoff) 
filtered_gini$dgini.outlier <- filtered_gini$dgini.outlier[,"outlier"] 
iq = qc_alpha(filtered_gini, cutoff_percentile=cutoff_percentile, cutoff_value=alpha_cutoff, hard_cutoff = alpha_hard_cutoff) 
iq$outlier = iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier


print(iq$dmapd.outlier)
print(iq$dgini.outlier)
print(iq$alpha.outlier)
print(iq$outlier)


report <- paste0("Out of ",length(df$replicating)," cells, ",sum((df$replicating)), " are replicating. Out of the ",sum((!df$replicating))," non-replicating cells, mapd outliers: ", sum(iq$dmapd.outlier)," gini outliers: ", sum(iq$dgini.outlier) , ", and alpha outliers: ",sum(iq$alpha.outlier), " and cells that are outlier at least in one of these categories: ",sum(iq$outlier))
writeLines(report, paste0(OUTPUT, "/report.txt"))

general_analysis_plots(df,iq,OUTPUT)

non_outlier_cells = iq[iq$outlier==FALSE,]$name
condition_to_stay <- pData(non_rep_object)$name %in% non_outlier_cells
non_outlier_object <- non_rep_object[, condition_to_stay]
#saveRDS(non_outlier_object@assayData$copynumber,paste0(OUTPUT,"/",SAMPLENAME,"_copynumber.rds"))
saveRDS(non_outlier_object,paste0(OUTPUT,"/",SAMPLENAME,"_non_outlier.rds"))



source(file.path(BASEU,"R/scUnique.R"))
print("started")	

#print(Biobase::pData(scabsolute)$outlier)

#condition_to_remove <- iq$outlier==TRUE 
#non_outlier_scabsolute <- scabsolute[, !condition_to_remove]

#print(ncol(non_outlier_scabsolute))
jointSegmentationObject = jointSegmentationUpdateMEDICC(non_outlier_object, change_prob = 1e-1, windowSize=5,
                                                        mediccSizeCutoff=5, mediccDiffAmplitude=1.5,
                                                        splitPerChromosome=TRUE, ncores=4, HMMPATH=HMMPATH, BASEDIR="/home/schnei01/scUnique")

#HMMPATH='/home/schnei01/scUnique/scAbsolute/results_SLX-23525/'
#jointSegmentationObject = jointSegmentationUpdateMEDICC(scabsolute, change_prob = 1e-1, windowSize=5,
#							mediccSizeCutoff=5, mediccDiffAmplitude=1.5,
#							splitPerChromosome=TRUE, ncores=4, HMMPATH=NULL, BASEDIR="/opt")

print("Segmentation is done!")
saveRDS(jointSegmentationObject,file.path(OUTPUT,"jointSegmentationObject.rds"))





pdf(paste0(OUTPUT, "/dmapd_quantile.pdf"), width = 3, height = 2, pointsize=7) # Adjust width and height as needed
cutoff_quantile <- quantile(abs(iq$dmapd), 1 - cutoff_percentile)
hist(abs(iq$dmapd), breaks = 100, col = "skyblue", main = "Filtering based on MAPD",
     xlab = "Linear Regression Residual ", ylab = "Frequency")
abline(v = cutoff_quantile, col = "red", lty = 2, lwd = 2)
dev.off()


pdf(paste0(OUTPUT, "/dgini_quantile.pdf"), width = 3, height = 2, pointsize=7) # Adjust width and height as needed
cutoff_quantile <- quantile(abs(iq$dgini), 1 - cutoff_percentile)
hist(abs(iq$dgini), breaks = 100, col = "skyblue", main = "Filtering based on Gini",
     xlab = "Linear Regression Residual ", ylab = "Frequency")
abline(v = cutoff_quantile, col = "red", lty = 2, lwd = 2)
dev.off()

densityCutoff=0.1
kde2d_filter(df$rpc, df$gini_normalized, q=densityCutoff,n=1000, plt=TRUE,paste0(OUTPUT,"/",SAMPLENAME,"_contur_rpc&gini_normalized.pdf"),'RPC', 'Normalized gini')
kde2d_filter(df$rpc, df$mapd, q=densityCutoff,n=1000, plt=TRUE,paste0(OUTPUT,"/",SAMPLENAME,"_contur_rpc&mapd.pdf"),'RPC', 'MAPD')
kde2d_filter(df$rpc, df$hmm.alpha, q=densityCutoff,n=1000, plt=TRUE,paste0(OUTPUT,"/",SAMPLENAME,"_contur_alpha.pdf"),'RPC', 'Alpha')


plotDataset(object=scabsolute, file=paste0(OUTPUT,"/out.pdf"), f=plotCopynumber,verbose=FALSE)
plotDataset(object=non_rep_object, file=paste0(OUTPUT,"/non_replicating_out.pdf"), f=plotCopynumber,verbose=FALSE)
plotDataset(object=non_outlier_object, file=paste0(OUTPUT,"/non_outlier_out.pdf"), f=plotCopynumber,verbose=FALSE)

#plot_heatmap(scabsolute,OUTPUT,SAMPLENAME,"all",col_cluster=FALSE)
#plot_heatmap(scabsolute,OUTPUT,SAMPLENAME,"all_samples",col_cluster=TRUE)

#plot_heatmap(non_rep_object,OUTPUT,SAMPLENAME,"nreplicating",col_cluster=FALSE)
#plot_heatmap(non_rep_object,OUTPUT,SAMPLENAME,"nreplicating_samples",col_cluster=TRUE)

#plot_heatmap(non_outlier_object,OUTPUT,SAMPLENAME,"noutlier",col_cluster=FALSE)
#plot_heatmap(non_outlier_object,OUTPUT,SAMPLENAME,"noutlier_samples",col_cluster=TRUE)


