
sample = "SLX-24130-50-qc_100"


UE <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".uniqueEvents.RDS"))
OUTPUT <- paste0("../../Volumes/Fl/scUnique_results/",sample,"/")
#dataframae for chr specific event frequency in cells
chromosomes <- c(paste0("chr", 1:24), "chrX", "chrY")
cell_counts <- data.frame(node = character(0), matrix(NA, nrow = 0, ncol = length(chromosomes) + 1))
colnames(cell_counts)[-1] <- chromosomes


#dataframae for specifications of each event in cells
events_df <- data.frame(node = character(0),matrix(NA, nrow = 0, ncol = 6))
colnames(events_df)[-1] <- c("start","end","length","h0","h1","h0toh1","chr")

ue_object<- UE
nodes <- names(ue_object)

for (j in 1:length(nodes)){
  ue <- ue_object[[nodes[j]]]$events
  ue_chr <- c()
  if(length(ue)>0){
    for (k in 1:length(ue)){
      ue_chr <- c(ue_chr, ue[[k]]$chr)
      
      events_df[nrow(events_df) + 1, ] <- NA  # Add a new row
      events_df[nrow(events_df), "node"] <- nodes[j]  # Assign the row name
      events_df[nrow(events_df), "start"] <- ue[[k]]$start
      events_df[nrow(events_df), "end"] <- ue[[k]]$end
      events_df[nrow(events_df), "length"] <- ue[[k]]$end-ue[[k]]$start
      events_df[nrow(events_df), "h0"] <- ue[[k]]$h0
      events_df[nrow(events_df), "h1"] <- ue[[k]]$h1
      events_df[nrow(events_df), "h0toh1"] <- paste0(ue[[k]]$h0," to ",ue[[k]]$h1)
      events_df[nrow(events_df), "chr"] <- ue[[k]]$chr
      
    }       
  }
  ue_chr_freq <- as.data.frame(table(ue_chr))
  
  cell_counts[nrow(cell_counts) + 1, ] <- NA  # Add a new row
  cell_counts[nrow(cell_counts), "node"] <- nodes[j]  # Assign the row name
  
  for (i in 1:nrow(ue_chr_freq)) {
    chromosome <- ue_chr_freq[i, "ue_chr"]
    frequency <- ue_chr_freq[i, "Freq"]
    column_index <- match(chromosome, colnames(cell_counts))
    cell_counts[nrow(cell_counts), column_index] <- frequency
  }
  cell_counts[is.na(cell_counts)] <- 0
  rownames(cell_counts)[nrow(cell_counts)] <- nodes[j]
}
cell_counts <- cell_counts[rowSums(is.na(cell_counts)) != ncol(cell_counts), ]

write.csv(events_df, paste0(OUTPUT,"ue_df.csv"), row.names = FALSE)

