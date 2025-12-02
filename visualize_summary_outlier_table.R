#install.packages("kableExtra")
#install.packages("webshot")
#install.packages("htmltools")
#install.packages("rmarkdown")
library(kableExtra)
library(webshot)
library(htmltools)

library(dplyr)
library(knitr)
library(gt)
library(webshot2)
library(kableExtra)
library(knitr)
library(kableExtra)

library(flextable)
library(officer)

# Define the samples and categories
ALLSAMPLES <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24174","24173","24489","24007","23965","24490","24491","24518","24835","24532")
categories <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","PEO1\nMISSENSE","PEO1\nSTOP","PEO1","PEO1")



# Define the samples and categories
ALLSAMPLES <- c("23355","25146","24831","24832","24833")
ALLSAMPLES=c("25146","25147","25149","24831","24832","24833")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
ALLSAMPLES=c("24518","24491")
categories <- c("PEO STOP","PEO MISSENSE")

ALLSAMPLES=c("24911","24912","24913")
categories <- c("FFPE","FFPE","FFPE")


ALLSAMPLES <- c("23962","23963","23964","24831","24832","24833","25146","25147","25148","25149","24911","24912","24913","24914","24915","24806","24807","24808","24809")
categories <- c("FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")

ALLSAMPLES=c("24757")
categories <- c("Organoid")

ALLSAMPLES=c("24757")
categories <- c("Organoid")


ALLSAMPLES=c("25393","25394","25394_B4","25394_C9","25394_D5","25394_D7","25394_D11","25394_G3")
categories <- c("PEO1 parent","PEO1 children","PEO1 B4","PEO1 C9","PEO1 D5","PEO1 B7","PEO1 D11","PEO1 G3")


bin_size <- "100"

# Initialize an empty data frame to store the combined outlier summaries
combined_outlier_summary <- data.frame()

for (i in 1:length(ALLSAMPLES)) {
  print(i)
  SAMPLENAME <- ALLSAMPLES[i]
  print(SAMPLENAME)
  SAMPLEONJ <- paste0("SLX-",SAMPLENAME,"_",bin_size)
  OUTPUT <- file.path("../../Volumes/Fl/sc_analysis/post-scAbsolute", SAMPLEONJ)
  
  # Read the outlier_summary for the current sample
  outlier_summary <- readRDS(paste0(OUTPUT, "/outlier_summary.rds"))
  print(outlier_summary)
  # Add columns for the sample name and category
  outlier_summary$Sample <- SAMPLENAME
  outlier_summary$Category <- categories[i]
  
  # Combine the current outlier summary with the accumulated data
  combined_outlier_summary <- rbind(combined_outlier_summary, outlier_summary)
}

# Optionally, set the row names to the sample names if desired
# rownames(combined_outlier_summary) <- combined_outlier_summary$Sample

# Define the custom order for categories
#categories <- c("PEO1", "PEO4", "PEO6", "PEO1\nMISSENSE", "PEO1\nSTOP", "PEO14", "PEO23", "CIOV3", "CIOV6", "HCT116", "HCT116\nBRCA2 -/-", "UWB1.289", "UWB1.289\nBRCA")

# Convert the Category column to a factor with specified levels
####xf$Category <- factor(combined_outlier_summary$Category, levels = categories)

# Order the dataframe by Category
sorted_data <- combined_outlier_summary %>%
  arrange(Category) %>%
  select(Sample, Category, everything())

sorted_data <- combined_outlier_summary
# Calculate the sum for numeric columns
sum_row <- sorted_data %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(Sample = "Total", Category = "")

# Add the sum row to the dataframe
sorted_data <- bind_rows(sorted_data, sum_row)


pdf("../../Volumes/Fl/sc_analysis/all_samples_23july2024/PEO1_evol_95.pdf", height = 6, width = 8)
print(kable(sorted_data, format = "latex", booktabs = TRUE) %>%
        kable_styling(latex_options = c("striped", "hold_position")))
dev.off()

# Create and save the table as PNG
kable(sorted_data, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  add_header_above(c(" " = 2, "Datasets" = ncol(sorted_data) - 2)) %>%
  row_spec(nrow(sorted_data), bold = TRUE, background = "#EFEFEF") %>%
  save_kable(file = "../../Volumes/Fl/sc_analysis/all_samples_23july2024/PEO1_evol_95.png", zoom = 2)

write.csv(sorted_data, "../../Volumes/Fl/sc_analysis/all_samples_23july2024/PEO1_evol.csv", row.names = FALSE)


ft <- flextable(sorted_data)

# Save as PDF
save_as_image(ft, path = "../../Volumes/Fl/sc_analysis/all_samples_23july2024/PEO1_evol.png")
pdf("../../Volumes/Fl/sc_analysis/all_samples_23july2024/PEO1_evol.pdf", height = 6, width = 8)
grid::grid.raster(png::readPNG("PEO1_evol.png"))
dev.off()
