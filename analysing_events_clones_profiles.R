library(dplyr)
source("R/core.R")
source("R/visualization_helpers.R")


sample <- "SLX-25393-15-qc_100"
sample <- "SLX-25394_B4-15-qc_100"

tree <- readRDS(paste0("/Volumes/Fl/scUnique_results/",sample,"/",sample,".tree.RDS"))
finalCNfreq <- readRDS(paste0("/Volumes/Fl/scUnique_results/",sample,"/",sample,".finalCN.RDS"))
background_freq <- readRDS(paste0("/Volumes/Fl/scUnique_results/",sample,"/",sample,".profiles_background_freq.RDS"))
dt <- as.dendrogram(tree)
plotUniqueEvents(finalCNfreq, background_freq, dt, difference_only = TRUE, fontsize_chr=22, fontsize_leg_title=22, fontsize_leg_label=18, row_title="")

finalCNfreq[, labels(as.dendrogram(tree))]

uniqueEventEncoding = matrix(-1, nrow=dim(finalCNfreq)[[1]], ncol=dim(finalCNfreq)[[2]])
uniqueEventIndices = finalCNfreq@assayData$copynumber != background_freq & !is.na(finalCNfreq@assayData$copynumber)
uniqueEventEncoding[uniqueEventIndices] = finalCNfreq@assayData$copynumber[uniqueEventIndices]
table(uniqueEventIndices)



df_pass_post <- readRDS("/Volumes/Fl/scUnique_results/SLX-25393-15-qc_100/SLX-25393-15-qc_100.df_pass_post.RDS")
df_pass_post_B4 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_B4-15-qc_100/SLX-25394_B4-15-qc_100.df_pass_post.RDS")
df_pass_post_C9 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_C9-15-qc_100/SLX-25394_C9-15-qc_100.df_pass_post.RDS")
df_pass_post_D5 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D5-15-qc_100/SLX-25394_D5-15-qc_100.df_pass_post.RDS")
df_pass_post_D7 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D7-15-qc_100/SLX-25394_D7-15-qc_100.df_pass_post.RDS")
df_pass_post_D11 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D11-15-qc_100/SLX-25394_D11-15-qc_100.df_pass_post.RDS")
df_pass_post_G3 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_G3-15-qc_100/SLX-25394_G3-15-qc_100.df_pass_post.RDS")


parent_C9 <- inner_join(df_pass_post, df_pass_post_C9, by = c("start", "end" ,"chroms_l")) 
parent_B4 <- inner_join(df_pass_post, df_pass_post_B4, by = c("start", "end", "chroms_l"))
parent_D5 <- inner_join(df_pass_post, df_pass_post_D5, by = c("start", "end", "chroms_l"))
parent_D7 <- inner_join(df_pass_post, df_pass_post_D7, by = c("start", "end", "chroms_l"))
parent_D11 <- inner_join(df_pass_post, df_pass_post_D11, by = c("start", "end", "chroms_l"))
parent_G3 <- inner_join(df_pass_post, df_pass_post_G3, by = c("start", "end", "chroms_l"))


exact_parent_C9 <- inner_join(df_pass_post, df_pass_post_C9, by = c("start", "end" ,"chroms_l","h0_l",  "h1_l")) 
exact_parent_B4 <- inner_join(df_pass_post, df_pass_post_B4, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_parent_D5 <- inner_join(df_pass_post, df_pass_post_D5, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_parent_D7 <- inner_join(df_pass_post, df_pass_post_D7, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_parent_D11 <- inner_join(df_pass_post, df_pass_post_D11, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_parent_G3 <- inner_join(df_pass_post, df_pass_post_G3, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))




exact_D5_G3 <- inner_join(df_pass_post_D5, df_pass_post_G3, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D5_D7 <- inner_join(df_pass_post_D5, df_pass_post_D7, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D5_C9 <- inner_join(df_pass_post_D5, df_pass_post_C9, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D5_D11 <- inner_join(df_pass_post_D5, df_pass_post_D11, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D5_B4 <- inner_join(df_pass_post_D5, df_pass_post_B4, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))

exact_D7_G3 <- inner_join(df_pass_post_D7, df_pass_post_G3, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D7_D5 <- inner_join(df_pass_post_D7, df_pass_post_D5, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D7_C9 <- inner_join(df_pass_post_D7, df_pass_post_C9, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D7_D11 <- inner_join(df_pass_post_D7, df_pass_post_D11, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D7_B4 <- inner_join(df_pass_post_D7, df_pass_post_B4, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))

exact_D11_G3 <- inner_join(df_pass_post_D11, df_pass_post_G3, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D11_D5 <- inner_join(df_pass_post_D11, df_pass_post_D5, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D11_C9 <- inner_join(df_pass_post_D11, df_pass_post_C9, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D11_D7 <- inner_join(df_pass_post_D11, df_pass_post_D7, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_D11_B4 <- inner_join(df_pass_post_D11, df_pass_post_B4, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))

exact_B4_G3 <- inner_join(df_pass_post_B4, df_pass_post_G3, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_B4_D5 <- inner_join(df_pass_post_B4, df_pass_post_D5, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_B4_C9 <- inner_join(df_pass_post_B4, df_pass_post_C9, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_B4_D7 <- inner_join(df_pass_post_B4, df_pass_post_D7, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_B4_D11 <- inner_join(df_pass_post_B4, df_pass_post_D11, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))

exact_C9_G3 <- inner_join(df_pass_post_C9, df_pass_post_G3, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_C9_D5 <- inner_join(df_pass_post_C9, df_pass_post_D5, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_C9_D11 <- inner_join(df_pass_post_C9, df_pass_post_D11, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_C9_D7 <- inner_join(df_pass_post_C9, df_pass_post_D7, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_C9_B4 <- inner_join(df_pass_post_C9, df_pass_post_B4, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))

exact_G3_C9 <- inner_join(df_pass_post_G3, df_pass_post_C9, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_G3_D5 <- inner_join(df_pass_post_G3, df_pass_post_D5, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_G3_D11 <- inner_join(df_pass_post_G3, df_pass_post_D11, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_G3_D7 <- inner_join(df_pass_post_G3, df_pass_post_D7, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))
exact_G3_B4 <- inner_join(df_pass_post_G3, df_pass_post_B4, by = c("start", "end", "chroms_l","h0_l",  "h1_l"))



D5_G3 <- inner_join(df_pass_post_D5, df_pass_post_G3, by = c("start", "end", "chroms_l"))
D5_D7 <- inner_join(df_pass_post_D5, df_pass_post_D7, by = c("start", "end", "chroms_l"))
D5_C9 <- inner_join(df_pass_post_D5, df_pass_post_C9, by = c("start", "end", "chroms_l"))
D5_D11 <- inner_join(df_pass_post_D5, df_pass_post_D11, by = c("start", "end", "chroms_l"))
D5_B4 <- inner_join(df_pass_post_D5, df_pass_post_B4, by = c("start", "end", "chroms_l"))

D7_G3 <- inner_join(df_pass_post_D7, df_pass_post_G3, by = c("start", "end", "chroms_l"))
D7_D5 <- inner_join(df_pass_post_D7, df_pass_post_D5, by = c("start", "end", "chroms_l"))
D7_C9 <- inner_join(df_pass_post_D7, df_pass_post_C9, by = c("start", "end", "chroms_l"))
D7_D11 <- inner_join(df_pass_post_D7, df_pass_post_D11, by = c("start", "end", "chroms_l"))
D7_B4 <- inner_join(df_pass_post_D7, df_pass_post_B4, by = c("start", "end", "chroms_l"))

D11_G3 <- inner_join(df_pass_post_D11, df_pass_post_G3, by = c("start", "end", "chroms_l"))
D11_D5 <- inner_join(df_pass_post_D11, df_pass_post_D5, by = c("start", "end", "chroms_l"))
D11_C9 <- inner_join(df_pass_post_D11, df_pass_post_C9, by = c("start", "end", "chroms_l"))
D11_D7 <- inner_join(df_pass_post_D11, df_pass_post_D7, by = c("start", "end", "chroms_l"))
D11_B4 <- inner_join(df_pass_post_D11, df_pass_post_B4, by = c("start", "end", "chroms_l"))

B4_G3 <- inner_join(df_pass_post_B4, df_pass_post_G3, by = c("start", "end", "chroms_l"))
B4_D5 <- inner_join(df_pass_post_B4, df_pass_post_D5, by = c("start", "end", "chroms_l"))
B4_C9 <- inner_join(df_pass_post_B4, df_pass_post_C9, by = c("start", "end", "chroms_l"))
B4_D7 <- inner_join(df_pass_post_B4, df_pass_post_D7, by = c("start", "end", "chroms_l"))
B4_D11 <- inner_join(df_pass_post_B4, df_pass_post_D11, by = c("start", "end", "chroms_l"))

C9_G3 <- inner_join(df_pass_post_C9, df_pass_post_G3, by = c("start", "end", "chroms_l"))
C9_D5 <- inner_join(df_pass_post_C9, df_pass_post_D5, by = c("start", "end", "chroms_l"))
C9_D11 <- inner_join(df_pass_post_C9, df_pass_post_D11, by = c("start", "end", "chroms_l"))
C9_D7 <- inner_join(df_pass_post_C9, df_pass_post_D7, by = c("start", "end", "chroms_l"))
C9_B4 <- inner_join(df_pass_post_C9, df_pass_post_B4, by = c("start", "end", "chroms_l"))

G3_C9 <- inner_join(df_pass_post_G3, df_pass_post_C9, by = c("start", "end", "chroms_l"))
G3_D5 <- inner_join(df_pass_post_G3, df_pass_post_D5, by = c("start", "end", "chroms_l"))
G3_D11 <- inner_join(df_pass_post_G3, df_pass_post_D11, by = c("start", "end", "chroms_l"))
G3_D7 <- inner_join(df_pass_post_G3, df_pass_post_D7, by = c("start", "end", "chroms_l"))
G3_B4 <- inner_join(df_pass_post_G3, df_pass_post_B4, by = c("start", "end", "chroms_l"))


# Group by three columns and summarize
grouped_df <- G3_C9 %>%
  group_by(start, end, chroms_l) %>%
  summarise( count = n(), .groups = "drop")
grouped_df

# Combine all dataframes into a list
df_list <- list(df_pass_post, df_pass_post_B4, df_pass_post_C9, df_pass_post_D5, df_pass_post_D7, df_pass_post_D11, df_pass_post_G3)

# Function to find matching rows between two dataframes
find_matching_rows <- function(df1, df2) {
  inner_join(df1, df2, by = c("start", "end", "chroms_l"))
}

# Initialize the first dataframe as the base
matching_rows <- df_list[[2]]

# Iterate over the remaining dataframes and find matching rows
for (i in 2:length(df_list)) {
  matching_rows <- find_matching_rows(matching_rows, df_list[[i]])
}

# Count the number of matching rows
num_matching_rows <- nrow(matching_rows)

# Print the matching rows and their count
print(matching_rows)
print(paste("Number of matching rows:", num_matching_rows))




tree_profiles_B4 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_B4-15-qc_100/SLX-25394_B4-15-qc_100.tree_profiles.RDS")

medicc_tree_B4 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_B4-15-qc_100/SLX-25394_B4-15-qc_100.medicc_tree.RDS")

plot(medicc_tree_B4)


tree_profiles <- readRDS("/Volumes/Fl/scUnique_results/SLX-25393-15-qc_100/SLX-25393-15-qc_100.tree_profiles.RDS")
tree_profiles_B4 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_B4-15-qc_100/SLX-25394_B4-15-qc_100.tree_profiles.RDS")
tree_profiles_C9 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_C9-15-qc_100/SLX-25394_C9-15-qc_100.tree_profiles.RDS")
tree_profiles_D5 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D5-15-qc_100/SLX-25394_D5-15-qc_100.tree_profiles.RDS")
tree_profiles_D7 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D7-15-qc_100/SLX-25394_D7-15-qc_100.tree_profiles.RDS")
tree_profiles_D11 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D11-15-qc_100/SLX-25394_D11-15-qc_100.tree_profiles.RDS")
tree_profiles_G3 <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_G3-15-qc_100/SLX-25394_G3-15-qc_100.tree_profiles.RDS")

tree_profiles[[30]]

for (i in 15:29){
  print(sum(!(tree_profiles_D7[[1]][[4]] == tree_profiles[[i]][[4]])))
  #print(sum(abs(tree_profiles_D7[[1]][[5]] - tree_profiles[[i]][[5]])))
}

finalCN <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D7-15-qc_100/SLX-25394_D7-15-qc_100.finalCN.RDS")

pData(finalCN)["SINCEL-283_SC_GEMNEG_Plate_430_D6",]

finalCN@assayData$copynumber[,"SINCEL-283_SC_GEMNEG_Plate_430_D6"]

# finalCN@assayData$copynumber[,"SINCEL-274_Plate_414_P22"]  # Sample not found in this dataset


