get_summary_of_outliers <- function(object,
                                    SLX = NULL,
                                    UID = NULL,
                                    # QC thresholds
                                    rpc_cutoff = 25,
                                    mapd_cutoff = 2,
                                    mapd_density_control = FALSE,
                                    gini_norm_cutoff = 2,
                                    alpha_cutoff = 1.5,
                                    alpha_hard_cutoff = 0.05,
                                    gini_density_control = TRUE,
                                    densityCutoff = 0.1,
                                    cutoff_percentile = 0.05,
                                    # Replicating cell detection
                                    replicating_cutoff_value = 1.5,
                                    replicating_iqr_value = 1.5,
                                    # Normal cell threshold
                                    normal_threshold = 95,
                                    # Cell cycle plot
                                    generate_cellcycle_plot = FALSE,
                                    cellcycle_path = NULL,
                                    # Output options
                                    save_results = FALSE,
                                    outlier_cellbase_path = NULL,
                                    non_outlier_obj_path = NULL,
                                    outlier_summary_path = NULL,
                                    normals_path = NULL) {

  # Build phenotype dataframe for QC functions
  modified_df <- Biobase::pData(object) %>%
    # Remove trailing underscore, if there is one
    dplyr::mutate(name = sub("_[0-9]*$", "", name)) %>%
    tidyr::separate(name, into = c("TECHNOLOGY", "cellid"), sep = "_(?=[^_]+$)", remove = FALSE)
  if (!is.null(SLX)) modified_df$SLX <- SLX
  if (!is.null(UID)) modified_df$UID <- UID

  # Pass the modified data frame to predict_replicating (uses batch grouping internally)
  df <- predict_replicating(
    modified_df,
    batch = "technology",
    cutoff_value = replicating_cutoff_value,
    iqr_value = replicating_iqr_value
  )

  # Optionally generate cell cycle plot
  if (isTRUE(generate_cellcycle_plot) && !is.null(cellcycle_path)) {
    fig_cellcycle <- ggplot2::ggplot(data = df) +
      ggbeeswarm::geom_quasirandom(ggplot2::aes(x = SLX, y = cycling_activity, color = replicating), size = 0.8, alpha = 1.0) +
      ggplot2::facet_wrap(~UID, scales = "free_x") +
      ggpubr::theme_pubclean()
    ggplot2::ggsave(cellcycle_path, fig_cellcycle, width = 3, height = 4)
  }
  
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
  # use parameters passed to the function (defaults above)
  
  #print(df)
  #print(sum(df$replicating==FALSE))
  filtered_rpc <- df %>% dplyr::filter(!replicating, rpc >= rpc_cutoff)
  
  #print(filtered_rpc)
  filtered_mapd <- qc_mapd(filtered_rpc, cutoff_percentile = cutoff_percentile, cutoff_value = mapd_cutoff, symmetric = FALSE,
                           densityControl = mapd_density_control, densityCutoff = densityCutoff)
  #print(filtered_mapd$dmapd.outlier)
  #filtered_mapd$dmapd.outlier <- filtered_mapd$dmapd.outlier[,"outlier"]
  
  filtered_gini <- qc_gini_norm(filtered_mapd, cutoff_value = gini_norm_cutoff,
                                densityControl = gini_density_control, densityCutoff = densityCutoff)
  # qc_gini_norm returns a data.frame with a dgini.outlier column (data.frame with outlier), normalize it
  if (is.data.frame(filtered_gini$dgini.outlier)) filtered_gini$dgini.outlier <- filtered_gini$dgini.outlier[, "outlier"]

  iq <- qc_alpha(filtered_gini, cutoff_percentile = cutoff_percentile, cutoff_value = alpha_cutoff, hard_cutoff = alpha_hard_cutoff)
  

  # create logical outlier column (flatten data.frame columns if necessary)
  if (is.data.frame(iq$dmapd.outlier)) dmap_flag <- iq$dmapd.outlier[, 1] else dmap_flag <- iq$dmapd.outlier
  if (is.data.frame(iq$dgini.outlier)) dgini_flag <- iq$dgini.outlier[, 1] else dgini_flag <- iq$dgini.outlier
  if (is.data.frame(iq$alpha.outlier)) alpha_flag <- iq$alpha.outlier[, 1] else alpha_flag <- iq$alpha.outlier
  iq$outlier <- dmap_flag | dgini_flag | alpha_flag
  
  message("Out of ", nrow(df), " cells, ", sum(df$replicating), " are replicating. Out of the ", sum(!df$replicating), " non-replicating cells: mapd outliers: ", sum(iq$dmapd.outlier), ", gini outliers: ", sum(iq$dgini.outlier), ", rpc outliers: ", sum(df[df$replicating == FALSE, ]$rpc < rpc_cutoff), ", alpha outliers: ", sum(iq$alpha.outlier), ", total outliers (at least one category): ", sum(iq$outlier))
  non_outlier_cells <- iq[!iq$outlier, ]$name
  condition_to_stay <- Biobase::pData(non_rep_object)$name %in% non_outlier_cells

  non_outlier_object <- non_rep_object[, condition_to_stay]
  if (isTRUE(save_results) && !is.null(non_outlier_obj_path)) saveRDS(non_outlier_object, non_outlier_obj_path)
  
  
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
  #print(metrics_df)
  normals= metrics_df[metrics_df$percentage_2>normal_threshold,]
  
  
  not_rep = df[df$replicating==FALSE,]
  rep_only = df[df$replicating==TRUE,]
  not_rep[not_rep$rpc<rpc_cutoff,]$cellid
  
  
  #summary_df <- as.data.frame(t(c(nrow(df),sum((df$replicating)),sum(iq$dmapd.outlier),sum(iq$dgini.outlier),sum(df[df$replicating==FALSE,]$rpc<15),sum(iq$alpha.outlier),ncol(non_outlier_object@assayData$copynumber),nrow(normals))))
  #colnames(summary_df) <- c("Processed Cells", "Replicating", "Mapd Outliers", "Gini Outliers","RPC Outliers","Alpha Outliers","Good Quality Cells","Normal Cells")
  
  
  summary_df <- as.data.frame(t(c(nrow(df),sum((df$replicating)),sum(rep_only$rpc<rpc_cutoff),sum(df[df$replicating==FALSE,]$rpc<rpc_cutoff),sum(iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier),ncol(non_outlier_object@assayData$copynumber),nrow(normals))))
  colnames(summary_df) <- c("post-scAbsolute", "Replicating", "Replicating(low-RPC)", "Outliers(RPC)", "Outliers(Alpha/Mapd/Gini,post-RPC)", "PassedQC(incl.Normal)", "Normal")
  
  # Plot the table
  #table_plot <- tableGrob(summary_df)
  
  if (isTRUE(save_results) && !is.null(normals_path)) saveRDS(normals, normals_path)
  if (isTRUE(save_results) && !is.null(outlier_summary_path)) saveRDS(summary_df, outlier_summary_path)
  
  # Create the summary dataframe
  if (length(unique(df$cellid))>1){
    summary_df_outliers <- data.frame(cellid = all_cellids) %>%
      mutate(
        replicating = cellid %in% df[df$replicating,]$cellid,
        dmap_outlier = cellid %in% iq[iq$dmapd.outlier,]$cellid,
        gini_outlier = cellid %in% iq[iq$dgini.outlier,]$cellid,
        rpc_outlier_non_replicating = cellid %in% not_rep[not_rep$rpc<rpc_cutoff,]$cellid,
        rpc_outlier_replicating = cellid %in% rep_only[rep_only$rpc<rpc_cutoff,]$cellid,
        alpha_outlier = cellid %in% iq[iq$alpha.outlier,]$cellid,
        na_alpha_outliers = cellid %in% na_alpha_cells$cellid,
        gini_alpha_mapd = cellid %in% iq[iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier,]$cellid
      ) %>%
      mutate(across(-cellid, as.integer))
  }else if((length(unique(df$name))>1)){
    summary_df_outliers <- data.frame(cellid = all_cellids) %>%
      mutate(
        replicating = cellid %in% df[df$replicating,]$name,
        dmap_outlier = cellid %in% iq[iq$dmapd.outlier,]$name,
        gini_outlier = cellid %in% iq[iq$dgini.outlier,]$name,
        rpc_outlier_non_replicating = cellid %in% not_rep[not_rep$rpc<rpc_cutoff,]$cellid,
        rpc_outlier_replicating = cellid %in% rep_only[rep_only$rpc<rpc_cutoff,]$cellid,
        alpha_outlier = cellid %in% iq[iq$alpha.outlier,]$name,
        na_alpha_outliers = cellid %in% na_alpha_cells$name,
        gini_alpha_mapd = cellid %in% iq[iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier,]$name
      ) %>%
      mutate(across(-cellid, as.integer))
  } else {
    print("There is a problem with the name of the cells")
  }
  # Get a summary of how many cells are in each category
  summary_counts <- summary_df_outliers %>%
    summarise(across(-cellid, sum))
  
  # optionally persist the cell-based outlier table
  if (isTRUE(save_results) && !is.null(outlier_cellbase_path)) saveRDS(summary_df_outliers, outlier_cellbase_path)

  # Store QC criteria used
  qc_criteria <- list(
    rpc_cutoff = rpc_cutoff,
    mapd_cutoff = mapd_cutoff,
    mapd_density_control = mapd_density_control,
    gini_norm_cutoff = gini_norm_cutoff,
    alpha_cutoff = alpha_cutoff,
    alpha_hard_cutoff = alpha_hard_cutoff,
    gini_density_control = gini_density_control,
    densityCutoff = densityCutoff,
    cutoff_percentile = cutoff_percentile,
    replicating_cutoff_value = replicating_cutoff_value,
    replicating_iqr_value = replicating_iqr_value,
    normal_threshold = normal_threshold
  )

  # return results as a list (pure function unless save_results requested)
  return(list(
    df = df,
    iq = iq,
    summary_df = summary_df,
    summary_df_outliers = summary_df_outliers,
    non_outlier_object = non_outlier_object,
    normals = normals,
    summary_counts = summary_counts,
    qc_criteria = qc_criteria
  ))
}