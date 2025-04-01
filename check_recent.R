

PEO_50 <- paste0("../../Volumes/Fl/scUnique_results/PEO50_qc_100/")
PEO1_50 <- paste0("../../Volumes/Fl/scUnique_results/SLX-23965-50-qc_100/")
PEO4_50 <- paste0("../../Volumes/Fl/scUnique_results/SLX-24077-50-qc_100/")
PEO6_50 <- paste0("../../Volumes/Fl/scUnique_results/SLX-24130-50-qc_100/")

# Read the CSV file into a dataframe
df_peo <- read.csv(paste0(PEO_50,"ue_df.csv"))
df_peo1 <- read.csv(paste0(PEO1_50,"ue_df.csv"))
df_peo4 <- read.csv(paste0(PEO4_50,"ue_df.csv"))
df_peo6 <- read.csv(paste0(PEO6_50,"ue_df.csv"))

common_events <- function(df1,df2){
  # Remove the 'node' column before comparing
  df1_filtered <- select(df1, -node)
  df2_filtered <- select(df2, -node)
  
  # Find common rows
  common_rows <- inner_join(df1_filtered, df2_filtered, by = colnames(df1_filtered))
  
  # Count the number of matching rows
  num_common_rows <- nrow(common_rows)
  print(num_common_rows)
}

exact_common_events <- function(df1,df2){

  # Find common rows
  common_rows <- inner_join(df1, df2, by = colnames(df1))
  
  # Count the number of matching rows
  num_common_rows <- nrow(common_rows)
  print(num_common_rows)
}

nrow(df_peo)
nrow(df_peo1)
nrow(df_peo4)
nrow(df_peo6)
common_events(df_peo1,df_peo)
exact_common_events(df_peo1,df_peo)
exact_common_events(df_peo4,df_peo)
exact_common_events(df_peo6,df_peo)





# Load required libraries
library(dplyr)

# Sample DataFrames (Replace with your actual data)
df1 <- df_peo1
df2 <- df_peo4
df3 <- df_peo6
df_merged <- df_peo

# Function to extract sample name from node
extract_sample <- function(node) {
  strsplit(node, "_")[[1]][1]  # Extract first part (e.g., "SINCEL-211")
}

# Add Sample Column
df1$sample <- sapply(df1$node, extract_sample)
df2$sample <- sapply(df2$node, extract_sample)
df3$sample <- sapply(df3$node, extract_sample)
df_merged$sample <- sapply(df_merged$node, extract_sample)

# Merge all sample dataframes
df_all_samples <- bind_rows(df1, df2, df3)

# Create sets for event comparison (using chr, start, and end)
samples_events <- unique(df_all_samples[, c("chr", "start", "end")])
merged_events <- unique(df_merged[, c("chr", "start", "end")])

# Find overlaps and new events
common_events <- intersect(merge(samples_events, merged_events, by = c("chr", "start", "end")), samples_events)
new_events <- setdiff(merged_events, samples_events)
missing_events <- setdiff(samples_events, merged_events)

# Find shared events across multiple samples
df_all_samples$event <- paste(df_all_samples$chr, df_all_samples$start, df_all_samples$end, sep = "_")
event_counts <- df_all_samples %>%
  group_by(event) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  filter(n_samples > 1)

# Calculate how many and what percentage of unique events in the merged appear in each df1, df2, df3
merged_event_count <- nrow(merged_events)

df1_events <- unique(df1[, c("chr", "start", "end")])
df2_events <- unique(df2[, c("chr", "start", "end")])
df3_events <- unique(df3[, c("chr", "start", "end")])

df1_common <- nrow(intersect(df1_events, merged_events))
df2_common <- nrow(intersect(df2_events, merged_events))
df3_common <- nrow(intersect(df3_events, merged_events))

df1_percentage_merged <- (df1_common / merged_event_count) * 100
df2_percentage_merged <- (df2_common / merged_event_count) * 100
df3_percentage_merged <- (df3_common / merged_event_count) * 100

# Calculate percentage of events in df1, df2, and df3 that are found in merged
df1_percentage_sample <- (df1_common / nrow(df1_events)) * 100
df2_percentage_sample <- (df2_common / nrow(df2_events)) * 100
df3_percentage_sample <- (df3_common / nrow(df3_events)) * 100

# Print results
cat("Total events in merged:", merged_event_count, "\n")
cat("Events in merged that were in individual samples:", nrow(common_events), "\n")
cat("New events found in merged:", nrow(new_events), "\n")
cat("Events from individual samples that disappeared in merged:", nrow(missing_events), "\n")
cat("Events shared across multiple samples:", nrow(event_counts), "\n")

# Print event details
# cat("\nNew events in merged:\n")
# print(new_events)

# cat("\nEvents shared across samples:\n")
# print(event_counts)

# Print the number and percentage of unique events in merged appearing in each sample
cat("\nEvents from merged appearing in df1:", df1_common, "(", df1_percentage_merged, "%)\n")
cat("Events from merged appearing in df2:", df2_common, "(", df2_percentage_merged, "%)\n")
cat("Events from merged appearing in df3:", df3_common, "(", df3_percentage_merged, "%)\n")

# Print the percentage of events in df1, df2, df3 that appeared in merged
cat("\nPercentage of events in df1 appearing in merged:", df1_percentage_sample, "%\n")
cat("Percentage of events in df2 appearing in merged:", df2_percentage_sample, "%\n")
cat("Percentage of events in df3 appearing in merged:", df3_percentage_sample, "%\n")


# Define the percentages
percentages <- c(df1_percentage_sample, df2_percentage_sample, df3_percentage_sample)
names <- c("df1", "df2", "df3")

pdf("../../Volumes/Fl/sc_analysis/all_samples_23july2024/percentage_barplot.pdf", width = 2, height = 5)  # Adjust size if needed
barplot(percentages, names.arg = names, col = c("#FFB347", "#77DD77", "#779ECB"),
        main = "Percentage of Events in Merged", ylab = "Percentage", ylim = c(0, 100), cex.names = 0.8)
text(x = seq_along(percentages), y = percentages, labels = paste0(round(percentages, 1), "%"), pos = 3, cex = 0.8)
dev.off() 

