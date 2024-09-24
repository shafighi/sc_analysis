
library(gridExtra)
library(pheatmap) 
library(Biobase)
library(gridExtra)
library(ggplot2)
library(dplyr)
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/sc_visualization.R")

current_directory <- getwd()
print(paste("Current Working Directory:", current_directory))


get_summary_of_outliers <- function(object,OUTLIERS_CELLBASE_PATH,NON_OUTLIER_OBJ_PATH,OUTLIER_SUMMARY_PATH,SLX_value,UID_value){
 
  # Modify the data frame before passing it to the predict_replicating function
  modified_df <- Biobase::pData(object) %>%
    tidyr::separate(name, into = c("TECHNOLOGY", "cellid"), sep = "_(?=[^_]+$)", remove = FALSE) %>%
    mutate(SLX = SLX_value, UID = UID_value)
  
  # Pass the modified data frame to predict_replicating
  df <- predict_replicating(
    modified_df,
    batch = "technology",
    cutoff_value = 1,
    iqr_value = 1.5
  )
  
  if (length(unique(df$cellid))>1){
    all_cellids <- unique(df$cellid)
  }else if((length(unique(df$name))>1)){
    print("Problem!!!")
    all_cellids <- unique(df$name)
  }else{
    print("Problem: I could not find the column of names!")
  }
  
  na_alpha_cells <- df[!complete.cases(df$hmm.alpha), ]
  df <- df[complete.cases(df$hmm.alpha), ]
  rep_cells = df[df$replicating==TRUE,]$name
  condition_to_remove <- df$name %in% rep_cells
  non_rep_object <- object[, !condition_to_remove]
  #iq = pre_scUnique_filtering(df[df$replicating==FALSE,],OUTPUT,mapd_cutoff=2,mapd_density_control= TRUE,gini_norm_cutoff=2,alpha_cutoff=2,alpha_hard_cutoff =0.05,gini_density_control=TRUE,rpc_cutoff=15)
  mapd_cutoff=1.5
  mapd_density_control= TRUE
  gini_norm_cutoff= 1.5
  alpha_cutoff=2
  alpha_hard_cutoff =0.05
  gini_density_control=TRUE
  rpc_cutoff=15
  densityCutoff = 0.1
  cutoff_percentile=0.05
  filtered_rpc <- df[df$replicating==FALSE,] %>% dplyr::filter(!replicating, rpc >= rpc_cutoff)
  
  filtered_mapd <- qc_mapd(filtered_rpc,cutoff_percentile=cutoff_percentile, cutoff_value = mapd_cutoff, symmetric=FALSE,
                           densityControl=mapd_density_control, densityCutoff=densityCutoff)
  filtered_mapd$dmapd.outlier <- filtered_mapd$dmapd.outlier[,"outlier"]
  filtered_gini <- qc_gini_norm(filtered_mapd, cutoff_value = gini_norm_cutoff, 
                                densityControl = gini_density_control, densityCutoff=densityCutoff)
  filtered_gini$dgini.outlier <- filtered_gini$dgini.outlier[,"outlier"]
  iq = qc_alpha(filtered_gini, cutoff_percentile=cutoff_percentile, cutoff_value=alpha_cutoff, hard_cutoff = alpha_hard_cutoff)
  iq$outlier = iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier
  
  print(paste0("Out of ",nrow(df)," cells, ",sum((df$replicating)), " are replicating. Out of the ",sum((!df$replicating))," non-replicating cells, mapd outliers: ", sum(iq$dmapd.outlier)," gini outliers: ", sum(iq$dgini.outlier), ", and rpc outliers: " , sum(df[df$replicating==FALSE,]$rpc<15), ", and alpha outliers: ",sum(iq$alpha.outlier), " and cells that are outlier at least in one of these categories: ",sum(iq$outlier)))
  
  non_outlier_cells = iq[iq$outlier==FALSE,]$name
  condition_to_stay <- pData(non_rep_object)$name %in% non_outlier_cells
  
  non_outlier_object <- non_rep_object[, condition_to_stay]
  saveRDS(non_outlier_object,NON_OUTLIER_OBJ_PATH)
  
  
  print("------------------ NORMAL CELLS ------------------")
  # Apply the function to each column (excluding the first one if it is 'Position')
  # Function to split row names into chromosome, start, and end
  
  calculate_metrics <- function(col) {
    total_values <- length(col)
    num_2 <- sum(col == 2, na.rm = TRUE)
    num_na <- sum(is.na(col))
    num_other <- total_values - num_2 - num_na
    percentage_2 <- (num_2 / (total_values - num_na)) * 100
    
    return(c(num_2, num_na, num_other, percentage_2))
  }
  split_row_names <- function(df) {
    parts <- strsplit(rownames(df), '[:-]')
    data.frame(
      chromosome = sapply(parts, '[', 1),
      start = as.integer(sapply(parts, '[', 2)),
      end = as.integer(sapply(parts, '[', 3))
    )
  }
  
  cn_non_out <- non_outlier_object@assayData$copynumber
  cn_non_out <- cbind(split_row_names(cn_non_out), cn_non_out)
  metrics <- sapply(cn_non_out[-1], calculate_metrics)
  
  # Convert the result to a dataframe and set column names
  metrics_df <- as.data.frame(t(metrics))
  colnames(metrics_df) <- c("num_2", "num_na", "num_other", "percentage_2")
  
  # Print the results
  print(metrics_df)
  normals= metrics_df[metrics_df$percentage_2>90,]
  
  summary_df <- as.data.frame(t(c(nrow(df),sum((df$replicating)),sum(iq$dmapd.outlier),sum(iq$dgini.outlier),sum(df[df$replicating==FALSE,]$rpc<15),sum(iq$alpha.outlier),ncol(non_outlier_object@assayData$copynumber),nrow(normals))))
  colnames(summary_df) <- c("Processed Cells", "Replicating", "Mapd Outliers", "Gini Outliers","RPC Outliers","Alpha Outliers","Good Quality Cells","Normal Cells")
  
  not_rep = df[df$replicating==FALSE,]
  not_rep[not_rep$rpc<15,]$cellid
  
  # Plot the table
  table_plot <- tableGrob(summary_df)
  
  saveRDS(summary_df, OUTLIER_SUMMARY_PATH)
  
  
  
  
  # Create the summary dataframe
  if (length(unique(df$cellid))>1){
  summary_df_outliers <- data.frame(cellid = all_cellids) %>%
    mutate(
      replicating = cellid %in% df[df$replicating,]$cellid,
      dmap_outlier = cellid %in% iq[iq$dmapd.outlier,]$cellid,
      gini_outlier = cellid %in% iq[iq$gini.outlier,]$cellid,
      rpc_outlier = cellid %in% not_rep[not_rep$rpc<15,]$cellid,
      alpha_outlier = cellid %in% iq[iq$alpha.outlier,]$cellid,
      na_alpha_outliers = cellid %in% na_alpha_cells$cellid
    ) %>%
    mutate(across(-cellid, as.integer))
  }else if((length(unique(df$name))>1)){
    summary_df_outliers <- data.frame(cellid = all_cellids) %>%
      mutate(
        replicating = cellid %in% df[df$replicating,]$name,
        dmap_outlier = cellid %in% iq[iq$dmapd.outlier,]$name,
        gini_outlier = cellid %in% iq[iq$gini.outlier,]$name,
        rpc_outlier = cellid %in% not_rep[not_rep$rpc<15,]$name,
        alpha_outlier = cellid %in% iq[iq$alpha.outlier,]$name,
        na_alpha_outliers = cellid %in% na_alpha_cells$name
      ) %>%
      mutate(across(-cellid, as.integer))
  } else {
    print("There is a problem with the name of the cells")
  }
  # Get a summary of how many cells are in each category
  summary_counts <- summary_df_outliers %>%
    summarise(across(-cellid, sum))
  
  print(summary_counts)
  saveRDS(summary_df_outliers, OUTLIERS_CELLBASE_PATH)
  
}

ALLSAMPLES=c("23003","24078","24533","23359","23526","24173","24534","24174","24535","23527","24175","24536","23961","24441","24835","23965","24489","24007","24490","24076","24518")

ALLSAMPLES <- c("24173","24489","24007","23965","24490","24174")
ALLSAMPLES=c("24077","24130","24491","23303","24532","23528")
ALLSAMPLES=c("23355","25146","25148","24831","24832","24833")

ALLSAMPLES=c("25146","25147","25149","24831","24832","24833")
ALLSAMPLES=c("24518","24491")
categories <- c("PEO STOP","PEO MISSENSE")

ALLSAMPLES=c("24911","24912","24913")
categories <- c("FFPE","FFPE","FFPE")

bin_size <- "100"

for (i in 1:length(ALLSAMPLES)){
  SAMPLENAME = ALLSAMPLES[i]
  SAMPLEONJ = paste0("SLX-",SAMPLENAME,"_",bin_size)
  ALLCELLS=file.path("Documents/sc_analysis/scAboslute-obj",paste0(SAMPLEONJ,'.rds')) 
  OUTPUT= file.path("Documents/sc_analysis/post-scAbsolute", SAMPLEONJ)
  if (file.exists(OUTPUT)) {
    cat("Directory is already there!", OUTPUT, "\n")
  } else {
    dir.create(OUTPUT)
    cat("Directory is created:", OUTPUT, "\n")
  }
  OUTLIERS_CELLBASE_PATH = paste0(OUTPUT,"/cellbased_outliers.rds")#summary_outliers.rds
  NON_OUTLIER_OBJ_PATH = paste0(OUTPUT,"/",SAMPLENAME,"_non_outlier.rds")
  OUTLIER_SUMMARY_PATH = paste0(OUTPUT,"/outlier_summary.rds")#summary_df.rds
  object = readRDS(ALLCELLS)
  SLX_value <- paste0("SLX-",SAMPLENAME)
  UID_value <- categories[i]
  get_summary_of_outliers(object,OUTLIERS_CELLBASE_PATH,NON_OUTLIER_OBJ_PATH,OUTLIER_SUMMARY_PATH,SLX_value,UID_value)
}

