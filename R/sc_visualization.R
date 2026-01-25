

#scUnique = "/home/schnei01/scUnique"
#BASEDIR="/home/schnei01/scAbsolute"
#source(file.path(BASEDIR, "R/load_dependencies.R"))
#load_dependencies(BASEDIR)
library(future, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(QDNAseq, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(ggbeeswarm, quietly = TRUE)
library("ggpubr", quietly = TRUE)
library(robustbase, quietly = TRUE)
library(readxl, quietly = TRUE)
library(dplyr, quietly = TRUE)

source(file.path("R/core.R"))
source(file.path( "R/visualization_helpers.R"))



pre_scUnique_filtering <- function(df,OUTPUT,mapd_cutoff=1.5,mapd_density_control= TRUE,gini_norm_cutoff= 1.5,alpha_cutoff=2,alpha_hard_cutoff =0.02,gini_density_control=TRUE,rpc_cutoff=15){
  filtered_rpc <- df %>% dplyr::filter(!replicating, rpc >= rpc_cutoff)
  
  filtered_mapd <- qc_mapd(filtered_rpc,cutoff_percentile=0.05, cutoff_value = mapd_cutoff, symmetric=FALSE,
                           densityControl=mapd_density_control, densityCutoff=0.10)
  filtered_mapd$dmapd.outlier <- filtered_mapd$dmapd.outlier[,"outlier"]
  filtered_gini <- qc_gini_norm(filtered_mapd, cutoff_value = gini_norm_cutoff, 
                                densityControl = gini_density_control)
  filtered_gini$dgini.outlier <- filtered_gini$dgini.outlier[,"outlier"]
  iq = qc_alpha(filtered_gini, cutoff_value=alpha_cutoff, hard_cutoff = alpha_hard_cutoff)
  iq$outlier = iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier
  iq$outlier <- iq$outlier[,"outlier"]
  print(iq$dmapd.outlier )
  fig_outliers = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=alpha, colour=outlier), 
                                                      size=0.8, alpha=1.0) + facet_wrap(~SLX, scales = "free_x")+
    theme_pubclean()
  fig_outliers
  ggsave(paste0(OUTPUT,"/cell_qc.jpg"), width = 3, height = 4)
  return(iq)
}

general_analysis_plots <- function(df,iq,OUTPUT){
  sample_name = df$UID[1]
  
  pdf(paste0(OUTPUT,"/normal_dist_alpha.pdf"), width = 3.5, height = 3, pointsize=7)
  cutoff_value <- 1.5
  x <- df$hmm.alpha
  hist(x, breaks = 60, col = "lightblue", main = "Filtering Based on Alpha Distribution",xlab=expression(alpha))
  abline(v = mean(x), col = "blue", lwd = 1)
  text(mean(x), par("usr")[4], labels = sprintf("mean(x).   %.2f", mean(x)),
       pos = 2, offset = 0.5, col = "blue", srt = 90, cex = 0.7)
  condition_region <- x [log(x) > mean(log(x)) + (cutoff_value * sd(log(x)))]
  abline(v = exp(mean(log(x)) + (cutoff_value * sd(log(x)))), col = "red", lwd = 1)
  text(exp(mean(log(x)) + (cutoff_value * sd(log(x)))), par("usr")[4], labels = sprintf("exp(mean(log(x)) + (cutoff_value * sd(log(x)))).   %.2f", exp(mean(log(x)) + (cutoff_value * sd(log(x))))),
       pos = 2, offset = 0.5, col = "red", srt = 90, cex = 0.7)
  abline(v = exp(mean(log(x)) + (sd(log(x)))), col = "green", lwd = 1)
  text(exp(mean(log(x)) + (sd(log(x)))), par("usr")[4], labels = sprintf("exp(mean(log(x)) + (sd(log(x)))).   %.2f", exp(mean(log(x)) + (sd(log(x))))),
       pos = 2, offset = 0.5, col = "green", srt = 90, cex = 0.7)
  hist(condition_region, breaks = 30, main = "Distribution with Filtering Condition", col = rgb(1, 0, 0, 0.5), add = TRUE,alpha=0.5)
  legend("topright", legend = c("Data", "Mean", "Filtering Region"),
         fill = c("lightblue", "blue", rgb(1, 0, 0, 0.5)), border = NA)
  dev.off()
  
  
  pdf(paste0(OUTPUT,"/total.reads_all.pdf"), width = 2, height = 3, pointsize=7)  # Adjust width and height as needed
  
  # Change the size of the font and figure
  par(
    cex.lab = 1,      # Set the size of axis labels
    cex.main = 1      # Set the size of the main title
  )
  
  # Create a histogram
  hist(df$total.reads, 
       main = "Distribution of Total Reads\n(all cells)",  # Title of the plot
       xlab = "Total Reads",                  # Label for x-axis
       ylab = "Probability",                    # Label for y-axis
       col = "skyblue",                       # Color of bars
       border = "black",                      # Border color of bars
       xlim = c(0, max(df$total.reads)),     # X-axis limits
       breaks = 40,                           # Number of bins
       prob = TRUE)                           # Density scale
  
  # Add a smoothed density curve
  lines(density(df$total.reads), col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Histogram", "Density"), fill = c("skyblue", "red"), bty = "n")
  
  dev.off()
  
  
  
  
  pdf(paste0(OUTPUT,"/total.reads_non_replicating.pdf"), width = 2, height = 3, pointsize=7)  # Adjust width and height as needed
  
  # Change the size of the font and figure
  par(
    cex.lab = 1,      # Set the size of axis labels
    cex.main = 1      # Set the size of the main title
  )
  
  # Create a histogram
  hist(iq$total.reads, 
       main = "Distribution of Total Reads\n(non-replicating)",  # Title of the plot
       xlab = "Total Reads",                  # Label for x-axis
       ylab = "Probability",                    # Label for y-axis
       col = "skyblue",                       # Color of bars
       border = "black",                      # Border color of bars
       xlim = c(0, max(iq$total.reads)),     # X-axis limits
       breaks = 40,                           # Number of bins
       prob = TRUE)                           # Density scale
  
  # Add a smoothed density curve
  lines(density(df$total.reads), col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Histogram", "Density"), fill = c("skyblue", "red"), bty = "n")
  
  dev.off()
  
  
  
  
  pdf(paste0(OUTPUT,"/total.reads_non_outliers.pdf"), width = 2, height = 3, pointsize=7)  # Adjust width and height as needed
  
  # Change the size of the font and figure
  par(
    cex.lab = 1,      # Set the size of axis labels
    cex.main = 1      # Set the size of the main title
  )
  
  # Create a histogram
  hist(iq$total.reads[iq$outlier==FALSE], 
       main = "Distribution of Total Reads\n(non-replicating)",  # Title of the plot
       xlab = "Total Reads",                  # Label for x-axis
       ylab = "Probability",                    # Label for y-axis
       col = "skyblue",                       # Color of bars
       border = "black",                      # Border color of bars
       xlim = c(0, max(iq$total.reads)),     # X-axis limits
       breaks = 40,                           # Number of bins
       prob = TRUE)                           # Density scale
  
  # Add a smoothed density curve
  lines(density(df$total.reads), col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Histogram", "Density"), fill = c("skyblue", "red"), bty = "n")
  
  dev.off()
  
  
  
  
  
  fig_cellcycle = ggplot(data = df) + geom_quasirandom(aes(x=SLX, y=cycling_activity,
                                                           color=replicating), size=0.8, alpha=1.0) +
    facet_wrap(~UID, scales = "free_x") + theme_pubclean()
  fig_cellcycle
  ggsave(paste0(OUTPUT,"/cellcycle.jpg"), width = 3, height = 4)
  
  
  print(iq)
  fig_outliers = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=alpha, colour=outlier), 
                                                      size=0.8, alpha=1.0) + facet_wrap(~SLX, scales = "free_x")+theme_pubclean()
  fig_outliers
  ggsave(paste0(OUTPUT,"/cell_qc.jpg"), width = 3, height = 4)
  
  
  
  
  fig_cellcycle = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=dmapd,color=dmapd.outlier), size=0.8, alpha=1.0) +
    facet_wrap(~UID, scales = "free_x") + theme_pubclean()
  fig_cellcycle
  ggsave(paste0(OUTPUT,"/dmapd_cutoff.jpg"), width = 3, height = 4)
  
  
  
  fig_cellcycle = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=dgini,
                                                           color=dgini.outlier), size=0.8, alpha=1.0) +
    facet_wrap(~UID, scales = "free_x") + theme_pubclean()
  fig_cellcycle
  ggsave(paste0(OUTPUT,"/dgini_cutoff.jpg"), width = 3, height = 4)
  
  
  
  

  ggplot(iq, aes(used.reads, rpc, color = ploidy)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Used Reads",
         y = "RPC",
         color="Ploidy")
  ggsave(paste0(OUTPUT,"/rpc_usedreads_ploidy_replicatings.jpg"), width = 8, height = 6)
  
  
  
  ggplot(iq, aes(used.reads, rpc, color = as.factor(outlier))) +
    geom_point(size=1) +
    
    # Discrete color scale
    scale_color_brewer(type = "qual", palette = "Set1") +
    
    theme_light() +
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Used Reads",
         y = "RPC",
         color = "Outlier")
  ggsave(paste0(OUTPUT,"/rpc_usedreads_outliers_replicatings.jpg"), width = 8, height = 6)
  
  
  
  ggplot(iq, aes(alpha, rpc, color = ploidy)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Alpha",
         y = "RPC",
         color="Ploidy")
  ggsave(paste0(OUTPUT,"/alpha_rpc_ploidy_replicatings.jpg"), width = 8, height = 6)
  
  
  ggplot(df, aes(rpc,alpha , color = as.factor(replicating))) +
    geom_point(size=1) +
    
    # Discrete color scale
    scale_color_brewer(type = "qual", palette = "Set1") +
    
    theme_light() +
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "RPC",
         y = "Alpha",
         color = "Replicating")
  ggsave(paste0(OUTPUT,"/alpha_rpc_replicating_all.jpg"), width = 8, height = 6)
  
  
  ggplot(df, aes(alpha,used.reads , color = as.factor(replicating))) +
    geom_point(size=1) +
    
    # Discrete color scale
    scale_color_brewer(type = "qual", palette = "Set1") +
    
    theme_light() +
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Alpha",
         y = "Used reads",
         color = "Replicating")
  ggsave(paste0(OUTPUT,"/alpha_used_reads.jpg"), width = 8, height = 6)
  
  

   
  
  
}


plot_heatmap <- function(object,OUTPUT,SAMPLENAME,post_name,col_cluster){
  cn = object@assayData$copynumber
  df_no_na <- na.omit(cn)
  
  png(paste0(OUTPUT,"/",SAMPLENAME,"_",post_name,".png"), width = 7000, height = 2000, res = 300)
  
  pheatmap(t(df_no_na), 
           cluster_rows=TRUE,
           cluster_cols=col_cluster,
           display_numbers = FALSE, 
           number_color="black",
           show_rownames = FALSE,
           labels_col = "",
           fontsize_row = 7, 
           fontsize_col = 7 )
  
  dev.off()
}

plot_heatmap_t <- function(object,OUTPUT,SAMPLENAME,post_name,row_cluster){
  cn = object@assayData$copynumber
  df_no_na <- na.omit(cn)
  df_no_na = t(df_no_na)
  
  png(paste0(OUTPUT,"/",SAMPLENAME,"_",post_name,".png"), width = 7000, height = 3500, res = 300)
  
  pheatmap(df_no_na, 
           cluster_rows=row_cluster,
           cluster_cols=TRUE,
           display_numbers = FALSE, 
           number_color="black",labels_row = "")
  dev.off()
}
