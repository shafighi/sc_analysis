library(dplyr)
library(forcats)
library(caret)
library(randomForest)
library(e1071)
library(class)
library(ggplot2)
library("flexmix")
library(data.table)
library(tidyr)
fitComponentGM<-function(dat,dist="norm",seed=77777,model_selection="BIC",min_prior=0.001,niter=1000,nrep=1,min_comp=2,max_comp=10)
{
  control<-new("FLXcontrol")
  control@minprior<-min_prior
  control@iter.max<-niter
  
  set.seed(seed)
  if(dist=="norm")
  {
    if(min_comp==max_comp)
    {
      fit<-flexmix::flexmix(dat ~ 1,model=flexmix::FLXMCnorm1(),k=min_comp,control=control)
    }else{
      fit<-flexmix::stepFlexmix(dat ~ 1,model = flexmix::FLXMCnorm1(),k=min_comp:max_comp,nrep=nrep,control=control)
      fit<-flexmix::getModel(fit,which=model_selection)
    }
    
  }else if(dist=="pois")
  {
    if(min_comp==max_comp)
    {
      fit<-flexmix::flexmix(dat ~ 1,model=flexmix::FLXMCmvpois(),k=min_comp,control=control)
    }else{
      fit<-flexmix::stepFlexmix(dat ~ 1,model = flexmix::FLXMCmvpois(),k=min_comp:max_comp,nrep=nrep,control=control)
      fit<-flexmix::getModel(fit,which=model_selection)
    }
  }
  fit
}



######################################################.     PARAMETERS.      ######################################################

seed=77777
min_comp=2
max_comp=20
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA")  # Adjust this to your preferred order

samples <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24174","24173","24489","24007","23965","24490","24831","24832","24833","25146","25147","25148","25149")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO23","PEO23","CIOV3","PEO4","CIOV6","PEO6","HCT116","CIOV6","PEO1","HCT116\nBRCA2 -/-","PEO14","UWB1.289\nBRCA","PEO1","PEO1","UWB1.289","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE","FFPE")
category_order <- c("PEO1", "PEO4","PEO6","CIOV3","CIOV6","PEO14","PEO23","HCT116","HCT116\nBRCA2 -/-","UWB1.289","UWB1.289\nBRCA","FFPE")  # Adjust this to your preferred order


samples <- c("23003","23303","23359","23526","24077","24130","24532","24007","23965")
category <- c("PEO1","PEO1","PEO1","PEO4","PEO4","PEO6","PEO1","PEO1","PEO1")
category_order <- c("PEO1", "PEO4","PEO6")  # Adjust this to your preferred order



all_samples_ue <- list()

# Loop to create sublists and add them to the main list
for (i in 1:length(samples)) {
  OUTPUT <- paste0('Documents/sc_analysis/scUnique-obj/SLX-',samples[i],'/')
  chromosomes <- c(paste0("chr", 1:22))
  UE_cells_PATH  <- paste0(OUTPUT,'unique_events_cells_compact.rds')
  ue_df <- readRDS(UE_cells_PATH)
  #ue_df_sum <- rowSums(ue_df[chromosomes])
  all_samples_ue[[i]] <- ue_df$length
}

# Assuming all_samples_ue is a list, convert it to a data frame
df <- data.frame(
  Sample = rep(samples, sapply(all_samples_ue, length)),
  Value = unlist(all_samples_ue)
)

# Create a mapping of Sample to Category
sample_categories <- data.frame(
  Sample = unique(df$Sample),
  Category = category
)

# Join the category information to the main dataframe
df <- df %>%
  left_join(sample_categories, by = "Sample")

df <- df %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)),
         group = interaction(Category, Sample, drop = TRUE, sep = "_"))


custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
                   "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB")
# Create the plot



dat<-as.numeric(df$Value)
ue_length_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                               min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

listComp = unlist(ue_length_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps['mean',], SD = allComps['sd',], Weight = prior(ue_length_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]


# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of Unique Events Lengths",
     xlab = "Length",
     ylab = "Frequency")

# Save the plot as a PNG file
png(paste0(OUTPUT,"changepoint_hist.png"), width = 500, height = 400, units = "px", res = 300)
dev.off()





# Step 2: Extract posterior probabilities
post_probs <- posterior(ue_length_mm)

# Step 3: Find the component with the highest posterior probability for each data point
assigned_components <- apply(post_probs, 1, which.max)
max_posteriors <- apply(posterior_probs, 1, max)

# View the assigned components
head(assigned_components)



# Assuming you have two dataframes:
# labels_df: contains the sample labels, with a column 'Label' and sample identifiers as row names
# posterior_df: contains posterior probabilities with samples as rows and components as columns

# Step 1: Combine the two dataframes by binding them together (assuming they have the same row order)
combined_df <- cbind(df, post_probs)

# Step 2: Sort the dataframe by the label column
combined_df_sorted <- combined_df[order(combined_df$group), ]

# Step 3: Remove the labels for plotting purposes (heatmap is just for posterior probabilities)
posteriors_sorted <- combined_df_sorted[, -1]  # Assuming the 'Label' is the first column

# Step 4: Plot the heatmap
library(pheatmap)

# Use pheatmap for creating the heatmap, displaying labels on the side
pheatmap(as.matrix(posteriors_sorted[c("1","2","3","4","5","6","7","8","9")]),
         cluster_rows = FALSE,  # Do not cluster rows, they are sorted by label
         cluster_cols = TRUE,   # You can cluster components (columns) if needed
         annotation_row = combined_df_sorted[, "group", drop = FALSE],  # Annotate rows by label
         show_rownames = TRUE,   # Show row names (samples)
         show_colnames = TRUE    # Show column names (components)
)

combined_df <- cbind(combined_df, assigned_components)

count_table <- table(combined_df$group, combined_df$assigned_components)
# Assuming you have the table result from the previous step:


# Step 1: Convert the table to a dataframe
count_df <- as.data.frame.matrix(count_table)

# Step 2: Normalize the counts (each row sums to 1)
normalized_df <- count_df %>%
  mutate(across(everything(), ~ . / rowSums(count_df)))

# Step 3: Plot the heatmap
library(pheatmap)

# Plot the heatmap, showing normalized values inside
pheatmap(as.matrix(normalized_df),
         cluster_rows = FALSE,   # Keep rows (samples) in the current order
         cluster_cols = TRUE,   # Disable clustering of columns
         display_numbers = TRUE, # Show the normalized values inside each cell
         number_format = "%.2f", # Format the numbers to 2 decimal places
         fontsize_number = 10,   # Font size for numbers inside cells
         color = colorRampPalette(c("white", "blue"))(50))  # Heatmap color gradient

