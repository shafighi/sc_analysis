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
bin_size <- "100"


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
  

  CN_PATH <- paste0(OUTPUT,'cn_',samples[s],'_',bin_size,'.rds')
  CN_BINNED_PATH <- paste0(OUTPUT,'cn_',samples[s],'_binned_',bin_size,'.rds')
  FEATURE_PATH = file.path(OUTPUT, paste0("cn_",samples[s],'_',bin_size,"_features.rds"))
  HEATMAP_PATH <- paste0(OUTPUT,"scAbsolute_copynumber.pdf")
  USEDREAD_RPC_PATH <- paste0(OUTPUT,"scUnique_rpc_usedreads_alpha_replicatings.jpg")
  USED_TOTAL_READ_PATH <- paste0(OUTPUT,"scUnique_usedreads_totalreads_alpha_replicatings.jpg")
  
  USEDREAD_RPC_HIST_PATH <- paste0(OUTPUT,"histogram_plot.pdf")
  DROPOUT_PATH_SEGMENT <- paste0(OUTPUT,'dropout_segment.pdf')
  DROPOUT_PATH_BIN <- paste0(OUTPUT,'dropout_bin.pdf')
  CN_FREQ_PATH <- paste0(OUTPUT,"cn_freq_scUnique.pdf")
  BP10MB_FREQ_PATH <- paste0(OUTPUT,"cn_bp10MB.pdf")
  osCN_FREQ_PATH <- paste0(OUTPUT,"cn_osCN.pdf")
  CHANGEPOINT_PATH <- paste0(OUTPUT,"cn_changepoint.pdf")
  SEGSIZE_PATH <- paste0(OUTPUT,"cn_segsize.pdf")
  UE_FREQ_PATH <- paste0(OUTPUT,"ue_freq.pdf")
  UE_HEATMAP_PATH <- paste0(OUTPUT,"ue_heatmap.pdf")
  UE_CHANGES_PATH <- paste0(OUTPUT,"ue_changes.pdf")
  
  
  
  
  print("------------------ LOADING OBJECTS ------------------")
  object<- readRDS(OBJ_PATH)
  cn <- object@assayData$copynumber
  object_df <- Biobase::pData(object)
  
  
  print("------------------ GENERATING COPYNUMBER DATAFRAME ------------------")
  
  # Function to split row names into chromosome, start, and end
  split_row_names <- function(df) {
    parts <- strsplit(rownames(df), '[:-]')
    data.frame(
      chromosome = sapply(parts, '[', 1),
      start = as.integer(sapply(parts, '[', 2)),
      end = as.integer(sapply(parts, '[', 3))
    )
  }
  
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
  
  
  cn_segmented <- cn_segmented[, c("chromosome","start","end","segVal", "sample" )]
  saveRDS(cn_segmented,CN_PATH)
  
  
  
  
  
  
  
  print("------------------ EXTRACTING FEATURES FROM COPYNUMBER PROFILE ------------------")
  
  # Set factors again so no empty samples are present. To be sure.
  cn_segmented$sample = factor(cn_segmented$sample)
  cn_segmented_list = split( cn_segmented, cn_segmented$sample )
  cn_segmented_list <- lapply(cn_segmented_list, as.data.frame)
  cn_features = extractCopynumberFeatures( cn_segmented_list, cores = CORES, prePath=PREPATH, rmNorm = RMNORM )
  saveRDS(cn_features, file = FEATURE_PATH )
  
  
  print("------------------ PLOTTING HEATMAP & READ-RPC RELATION ------------------")
  
  ## I should add plotting heatmap for scAbsolute
  col = c("1" = "red", "2" = "blue", "A" = "green", "B" = "orange")
  
  panel_fun = function(index, levels) {
    grid.rect(gp = gpar(fill = col[levels[2]], col = "black"))
    txt = paste(levels, collapse = ",")
    txt = paste0(txt, "\n", length(index), " rows")
    grid.text(txt, 0.5, 0.5, rot = 0,
              gp = gpar(col = col[levels[1]]))
  }
  
  plotCopynumberHeatmap(object,file = HEATMAP_PATH,show_chromosome_names=TRUE)
  ggplot(object_df, aes(used.reads, rpc, color = alpha)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = samples[s],
         x = "Used Reads",
         y = "RPC",
         color="Alpha")
  ggsave(USEDREAD_RPC_PATH, width = 8, height = 6)
  
  
  ggplot(object_df, aes(used.reads, total.reads, color = alpha)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = samples[s],
         x = "Used Reads",
         y = "Total Reads",
         color="Alpha")
  ggsave(USED_TOTAL_READ_PATH, width = 8, height = 6)
  
  
  # Create the histogram
  histogram <- ggplot(data = object@phenoData@data, aes(x = used.reads / rpc)) +
    geom_histogram(binwidth = 1, fill = "skyblue", color = "black", aes(y=..count..)) +
    labs(x = "used.reads / rpc", y = "Frequency", title = "Histogram of used.reads / rpc") +
    theme_minimal()
  
  # Save the plot as a PDF
  pdf(USEDREAD_RPC_HIST_PATH)
  print(histogram)
  dev.off()
  
  
  print("------------------ PLOTTING USED READS, RPC and PLOIDY ------------------")
  
  ggplot(object_df, aes(used.reads, ploidy, color = rpc)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = samples[s],
         x = "Used Reads",
         y = "Ploidy",
         color="RPC")
  ggsave(paste0(OUTPUT,"/scAbsolute/rpc_usedreads_ploidy_2.jpg"), width = 8, height = 6)
  
  
  print("------------------ PLOTTING DROPOUT ------------------")
  
  dropout <- as.data.frame(table(cn_segmented[cn_segmented$segVal==0,]$chromosome) )
  pdf(DROPOUT_PATH_SEGMENT, width = 6, height = 3)
  barplot(dropout$Freq, names.arg = dropout$Var1, xlab = "", ylab = "Dropout Frequency (segments)", col = "skyblue",
          main = paste0("SLX-",samples[s]), cex.axis = 0.4, cex.lab = 0.4, cex.names = 0.4)
  dev.off()
  
  
  dropout <- as.data.frame(table(cn_binned[cn_binned$segVal==0,]$chromosome) )
  pdf(DROPOUT_PATH_BIN, width = 6, height = 3)
  barplot(dropout$Freq, names.arg = dropout$Var1, xlab = "", ylab = "Dropout Frequency (bins)", col = "skyblue",
          main = paste0("SLX-",samples[s]), cex.axis = 0.4, cex.lab = 0.4, cex.names = 0.4)
  dev.off()
  
  
  
  print("------------------ VISUALIZE CN FREQUENCY ------------------")
  pdf(CN_FREQ_PATH, width = 8, height = 5)
  cn_features_numeric=as.numeric(cn_features[["copynumber"]][,2])
  cn_freq <- as.data.frame(table(cn_features_numeric))
  
  barplot(cn_freq$Freq, names.arg = cn_freq$dat, xlab = "", ylab = "Copy Number Frequency", col = "skyblue",
          main = paste0("PEO1: SLX-",samples[s]), cex.axis = 0.7, cex.lab = 0.7)
  dev.off()
  
  print("------------------ VISUALIZE BP10MB FREQ FREQUENCY ------------------")
  pdf(BP10MB_FREQ_PATH, width = 8, height = 10)
  bp10MB_numeric=as.numeric(cn_features[["bp10MB"]][,2])
  bp10MB_freq <- as.data.frame(table(bp10MB_numeric))
  barplot(bp10MB_freq$Freq, names.arg = bp10MB_freq$dat, xlab = "", ylab = "BreakPoint Frequency", col = "skyblue",
          main = paste0("SLX-",samples[s]), cex.axis = 0.7, cex.lab = 0.7)
  dev.off()
  
  
  
  print("------------------ VISUALIZE osCN FREQ FREQUENCY ------------------")
  pdf(osCN_FREQ_PATH, width = 8, height = 10)
  ocCN_numeric=as.numeric(cn_features[["osCN"]][,2])
  ocCN_freq <- as.data.frame(table(ocCN_numeric))
  barplot(ocCN_freq$Freq, names.arg = ocCN_freq$dat, xlab = "", ylab = "osCN Frequency", col = "skyblue",
          main = paste0("SLX-",samples[s]), cex.axis = 0.7, cex.lab = 0.7)
  dev.off()
  
  
  
  print("------------------ VISUALIZE CHANGEPOINT FREQ FREQUENCY ------------------")
  pdf(CHANGEPOINT_PATH, width = 8, height = 10)
  changepoint_numeric=as.numeric(cn_features[["changepoint"]][,2])
  changepoint_freq <- as.data.frame(table(changepoint_numeric))
  barplot(changepoint_freq$Freq, names.arg = changepoint_freq$dat, xlab = "", ylab = "changepoint Frequency", col = "skyblue",
          main = paste0("SLX-",samples[s]), cex.axis = 0.7, cex.lab = 0.7)
  dev.off()
  
  
  print("------------------ VISUALIZE SEGSIZE FREQ FREQUENCY ------------------")
  pdf(SEGSIZE_PATH, width = 8, height = 10)
  segsize_numeric=as.numeric(cn_features[["segsize"]][,2])
  segsize_freq <- as.data.frame(table(segsize_numeric))
  
  barplot(segsize_freq$Freq, names.arg = segsize_freq$dat, xlab = "", ylab = "Segsize Frequency", col = "skyblue",
          main = paste0("SLX-",samples[s]), cex.axis = 0.7, cex.lab = 0.7)
  dev.off()
  
}

extractCopynumberFeatures( cn_segmented_list, cores = CORES, prePath=PREPATH, rmNorm = RMNORM )

print(paste0(PREPATH, "hg19.chrom.sizes.txt"))
chrlen<-read.table(paste0(PREPATH, "hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
gaps<-read.table(paste0(PREPATH, "gap_hg19.txt"),sep="\t",header=F,stringsAsFactors = F)
centromeres<-gaps[gaps[,8]=="centromere",]


abs_profiles <- cn_segmented_list


out<-c()
samps<-getSampNames(abs_profiles)
for(i in samps)
{
  if(class(abs_profiles)=="QDNAseqCopyNumbers")
  {
    segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
  }else
  {
    segTab<-abs_profiles[[i]]
    colnames(segTab)[4]<-"segVal"
  }
  chrs<-unique(segTab$chromosome)
  all_dists<-c()
  for(c in chrs)
  {
    if(nrow(segTab)>1)
    {
      starts<-as.numeric(as.character(segTab[segTab$chromosome==c,2]$start))[-1]
      segstart<-as.numeric(as.character(segTab[segTab$chromosome==c,2]$start))[1]
      ends<-as.numeric(as.character(segTab[segTab$chromosome==c,3]$end))
      segend<-ends[length(ends)]
      ends<-ends[-length(ends)]
      centstart<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
      centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
      chrend<-chrlen[substr(chrlen[,1],4,5)==c,2]
      ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
      ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
      ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
      ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
      ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
      ndist<-apply(ndist,1,min)
      
      all_dists<-rbind(all_dists,sum(ndist>0))
      all_dists<-rbind(all_dists,sum(ndist<=0))
    }
  }
  if(nrow(all_dists)>0)
  {
    # Make sure it's really numeric
    out<-rbind(out,cbind(ID=i,ct1=as.numeric(all_dists[,1])))
  }
}
rownames(out)<-NULL
data.frame(out,stringsAsFactors = F)



