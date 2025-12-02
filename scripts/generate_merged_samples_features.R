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
#SINCEL-194 SINCEL_210 SINCEL-245 SINCEL-246 SINCEL-219 SINCEL-211
#PEO1 PEO4. BRCA NULL PEO4 PEO1 
#sample_names <- c("24007","23526","24489","24490","24077","23965")
sample_names <- c("23003","23303","23359","23526","23527","23528","23961","24077","24078","24130","24175","24441","24532","24173","24489","24007","23965","24490","24174")

chromosomes <- c(paste0("chr", 1:22))#, "chrX", "chrY")
INPUTFILES <- 'Documents/sc_analysis/scUnique-obj/'
OUTPUT <- paste0('Documents/sc_analysis/all_samples_23july2024/')
if (file.exists(OUTPUT)) {
  cat("Directory is already there!", OUTPUT, "\n")
} else {
  dir.create(OUTPUT, recursive = TRUE)
  cat("Directory is created:", OUTPUT, "\n")
}  

seed=77777
min_comp=2
max_comp=10
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000



######################################################     LOADING      ######################################################

# Initialize an empty list to store the sample data
samples <- list()

# Define bin_size
bin_size <- 100  # Replace with your actual bin size

# Loop through each sample name to read the data and store in the list
for (sample_name in sample_names) {
  print(paste0("SAMPLE NAME: ", sample_name))
  output1 <- paste0(INPUTFILES,'SLX-', sample_name, '/')
  feature_path <- paste0(output1, "cn_", sample_name, '_', bin_size, "_features.rds")
  
  sample_data <- readRDS(feature_path)
  sample_data$node <- paste0(sample_name, "_", sample_data$node)
  
  UE_cells_PATH  <- paste0(output1,'unique_events_cells.rds')
  sample_ue <- readRDS(UE_cells_PATH)
  sample_data$uechr <- pivot_longer(sample_ue, 
                                    cols = starts_with("chr"), 
                                    names_to = "chr", 
                                    values_to = "cnt")
  sample_data$uechr <- sample_data$uechr[, c("node", "cnt", "chr")]
  
  samples[[sample_name]] <- sample_data
}



# Initialize an empty list to store combined features
CN_features <- list()

# List of features to combine
features <- c("bp10MB", "osCN", "bpchrarm", "segsize", "changepoint", "copynumber","uechr")

# Loop through each feature and combine the corresponding data from all samples
for (feature in features) {
  # Use lapply to extract the feature from each sample and combine them with rbind
  CN_features[[feature]] <- do.call(rbind, lapply(samples, `[[`, feature))
}
CN_features$uechr <- as.data.frame(CN_features$uechr[grepl("SINCEL", CN_features$uechr$node), ])



thisOut = list()
dat<-as.numeric(CN_features[["segsize"]][,2])
segsize_mm<-fitComponentGM(dat,seed=seed,model_selection=model_selection,
                           min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

listComp = unlist(segsize_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps['mean',], SD = allComps['sd',], Weight = prior(segsize_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "segsize" ]] = dtSolution

# Save the plot as a PNG file


png(paste0(OUTPUT,"segsize_hist.png"),  width = 800, height = 1200,units = "px", res = 300)

hist(dat, col = "darkgreen",
     main = "Histogram of segsize",
     xlab = "segsize",
     ylab = "Frequency")
dev.off()


dat<-as.numeric(CN_features[["bp10MB"]][,2])
bp10MB_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

listComp = unlist(bp10MB_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps, Weight = prior(bp10MB_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "bp10MB" ]] = dtSolution

# Save the plot as a PNG file
png(paste0(OUTPUT,"bp10MB_hist.png"), width = 800, height = 1200, units = "px", res = 300)
# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of bp10MB",
     xlab = "bp10MB",
     ylab = "Frequency")
dev.off()




dat<-as.numeric(CN_features[["osCN"]][,2])
osCN_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

listComp = unlist(osCN_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps, Weight = prior(osCN_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "osCN" ]] = dtSolution

# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of osCN",
     xlab = "osCN",
     ylab = "Frequency")

# Save the plot as a PNG file
png(paste0(OUTPUT,"osCN_hist.png"), width = 500, height = 400, units = "px", res = 300)
dev.off()



dat<-as.numeric(CN_features[["bpchrarm"]][,2])
bpchrarm_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                            min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
listComp = unlist(bpchrarm_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps, Weight = prior(bpchrarm_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "bpchrarm" ]] = dtSolution

# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of bpchrarm",
     xlab = "bpchrarm",
     ylab = "Frequency")

# Save the plot as a PNG file
png(paste0(OUTPUT,"bpchrarm_hist.png"), width = 500, height = 400, units = "px", res = 300)
dev.off()



dat<-as.numeric(CN_features[["changepoint"]][,2])
changepoint_mm<-fitComponentGM(dat,seed=seed,model_selection=model_selection,
                               min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

listComp = unlist(changepoint_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps['mean',], SD = allComps['sd',], Weight = prior(changepoint_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "changepoint" ]] = dtSolution

# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of changepoint",
     xlab = "changepoint",
     ylab = "Frequency")

# Save the plot as a PNG file
png(paste0(OUTPUT,"changepoint_hist.png"), width = 500, height = 400, units = "px", res = 300)
dev.off()



dat<-as.numeric(CN_features[["copynumber"]][,2])
copynumber_mm<-fitComponentGM(dat,seed=seed,model_selection=model_selection,
                              nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000)
listComp = unlist(copynumber_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps['mean',], SD = allComps['sd',], Weight = prior(copynumber_mm))
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "copynumber" ]] = dtSolution

# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of copynumber",
     xlab = "copynumber",
     ylab = "Frequency")

# Save the plot as a PNG file
png(paste0(OUTPUT,"copynumber_hist.png"), width = 500, height = 400, units = "px", res = 300)
dev.off()







dat<-as.numeric(CN_features[["uechr"]]$cnt)
uechr_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                               nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000)
listComp = unlist(uechr_mm@components)
allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
dtSolution <- data.table(Mean = allComps, Weight = prior(uechr_mm))#, SD = allComps['sd',]
dtSolution = dtSolution[ order(dtSolution$Mean) ]
thisOut[[ "uechr" ]] = dtSolution

png(paste0(OUTPUT,"uechr_hist.png"), width = 800, height = 1200, units = "px", res = 300)
# Create the histogram
hist(dat, col = "darkgreen",
     main = "Histogram of unique events",
     xlab = "unique events",
     ylab = "Frequency")

# Save the plot as a PNG file

dev.off()


models <- list(uechr=uechr_mm,segsize=segsize_mm,bp10MB=bp10MB_mm,osCN=osCN_mm,bpchrarm=bpchrarm_mm,changepoint=changepoint_mm,copynumber=copynumber_mm)



###########################################################


# Calculate the sum-of-posterior matrix for the extracted copy number features
allEcnf = CN_features
allModels = thisOut
UNINFPRIOR=FALSE

allFeatures = names(allModels)
lMats = lapply(allFeatures, function(thisFeature) {
  
  print(thisFeature)
  thisEcnf = allEcnf[[ thisFeature ]]
  thisModel = allModels[[ thisFeature ]]
  
  dat = as.numeric( thisEcnf[,2] )
  # We want a posterior, hence likelihood (density) times prior (weight)
  if( ncol(thisModel) == 2 ) {
    # Poisson model
    print("Poisson-based posterior")
    if(UNINFPRIOR){
      postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
    } else {
      postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )    
    }
    
  } else {
    # Gaussian model
    print("Gauss-based posterior")
    if(UNINFPRIOR){
      postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
    } else {
      postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
    }
  }
  
  # Normalise densities to probabilities
  postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
  postDatScaled$Sample = thisEcnf[,1]
  matSxC = aggregate(. ~ Sample, postDatScaled, sum)
  rownames(matSxC) = matSxC$Sample
  matSxC$Sample = NULL
  matSxC = as.matrix(matSxC)
  
  # Should be sorted but just to be sure
  #matSxC = matSxC[ , order(thisModel[,"Mean"]) ]
  colnames(matSxC) = paste0( thisFeature, 1:ncol(matSxC) )
  
  return(matSxC)
  
} )

allMats = do.call(cbind, lMats)

saveRDS(allMats,paste0(OUTPUT,"sum_of_posterior.rds"))
