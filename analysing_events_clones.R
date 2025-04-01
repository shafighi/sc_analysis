


df_rcna_post <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.df_rcna_post.RDS")
df_rcna <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.df_rcna.RDS")
df_events <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.df_events.RDS")
df_rejected <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.df_rejected.RDS")
finalCN <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.finalCN.RDS")
n_events <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.n_events.RDS")
tree_profiles <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.tree_profiles.RDS")
uniqueEvents <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.uniqueEvents.RDS")
df_pass_post <- readRDS("../../Volumes/Fl/scUnique_results/PEO1-evol-15_100/PEO1-evol-15_100.df_pass_post.RDS")




# Initialize a list to store the count of events for each SINCEL
sincel_event_counts <- list()

# Iterate over each key in uniqueEvents
for (key in names(uniqueEvents)) {
  # Check if the key starts with 'SINCEL'
  if (startsWith(key, "SINCEL")) {
    # Extract the SINCEL name (e.g., 'SINCEL-283')
    sincel_name <- strsplit(key, "_")[[1]][1]
    # Count the number of events for this key
    event_count <- length(uniqueEvents[[key]]$events)
    # Add the event count to the corresponding SINCEL name
    if (is.null(sincel_event_counts[[sincel_name]])) {
      sincel_event_counts[[sincel_name]] <- 0
    }
    sincel_event_counts[[sincel_name]] <- sincel_event_counts[[sincel_name]] + event_count
  }
}

# Print the number of events for each SINCEL name
for (sincel_name in names(sincel_event_counts)) {
  cat(sincel_name, ": ", sincel_event_counts[[sincel_name]], " events\n", sep = "")
}


library(dplyr)
# Create entries in data frame for rows not in df_pass
df_allcells2 = df_pass_post %>%
  #dplyr::group_by(GROUP) %>%
  #dplyr::mutate(n_rcna_cells = length(unique(cellindex))) %>%
  #dplyr::group_by(UID, SLX, GROUP, TIMEPOINT, subcondition, cellindex, cellname) %>%
  dplyr::group_by(cellname)%>%
  dplyr::summarise(n=n(), has_rCNA=TRUE,# n_rcna_cells=n_rcna_cells,
                   total_size=sum(width), max_size=max(width), median_size=median(width),
                   dist_centromere_mean = mean(dist_centromere),
                   dist_telomere_max = max(dist_telomere),
                   dist_telomere_min = min(dist_telomere),
                   cn_peak = max(mean_reads / rpc),
                   cn_peak_count = sum((mean_reads / rpc)[as.numeric(h1_l) == 8] >= 10),
                   cn_peak_mean = mean((mean_reads / rpc)[as.numeric(h1_l) == 8]),
                   start_min=min(start), end_min=min(end),
                   start_max=max(start), end_max=max(end),
                   cn_jump = median(as.numeric(h1_l)-as.numeric(h0_l)),
                   cn_jump_max = max(as.numeric(h1_l)-as.numeric(h0_l)),
                   max_amplitude = max(abs(as.numeric(h1_l)-as.numeric(h0_l))),
                   high_amplitude_count = sum(as.numeric(h1_l) == 8),
                   median_amplitude=median(abs(as.numeric(h1_l)-as.numeric(h0_l))),
                   mean_amplitude=mean(abs(as.numeric(h1_l)-as.numeric(h0_l))))

similar = df_pass_post %>%
  #dplyr::group_by(GROUP) %>%
  #dplyr::mutate(n_rcna_cells = length(unique(cellindex))) %>%
  #dplyr::group_by(UID, SLX, GROUP, TIMEPOINT, subcondition, cellindex, cellname) %>%
  dplyr::group_by(start,end,chroms_l) %>%
  dplyr::summarise(n=n(),cells = cellname)
  
  
  
  
  