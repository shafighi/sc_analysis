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

#finalCN <- readRDS("/Volumes/Fl/scUnique_results/SLX-25394_D7-15-qc_100/SLX-25394_D7-15-qc_100.finalCN.RDS")

#pData(finalCN)["SINCEL-283_SC_GEMNEG_Plate_430_D6",]

#finalCN@assayData$copynumber[,"SINCEL-283_SC_GEMNEG_Plate_430_D6"]

# finalCN@assayData$copynumber[,"SINCEL-274_Plate_414_P22"]  # Sample not found in this dataset


