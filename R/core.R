## Distance for the scUnique package

# as the name says
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

## convert QDNAseq objects to object with invalid position removed
coord_to_valid <- function(index, object){
  valid = binsToUseInternal(object)
  rn = rownames(object)[index]
  new_index = which(rownames(object[valid, ]) %in% rn)
  return(new_index)
}

## convert from NA free (only valid) coordinates to QDNAseq object coordinates
valid_to_coord <- function(index, object){
  valid = binsToUseInternal(object)
  valid_index = which(cumsum(valid) == index)[[1]]
  return(valid_index)
}


# predict gini outliers cells in each dataset
qc_gini <- function(dft, factor_iqr=3, cutoff_value=2.0, cutoff_percentile=0.03){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("gini" %in% colnames(dft))
  
  # analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      
      # DEBUG
      #x = df_ploidy %>% dplyr::filter(UID %in% c("UID-NNA-TN3"))
      #x = df_ploidy %>% dplyr::filter(UID %in% c("UID-DLP-SA1090"))
      
      # compute mode of distribution
      Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }
      #hist(x$cycling_activity)
      x$gini_mode = Mode(round(x$gini, digits = 3))
      
      # sd of distribution (assume this is half-normal) below mode:
      sd_est = sqrt((1 / (1 - (2/pi))) * (sd(x$gini[x$gini < x$gini_mode], na.rm=TRUE)^2))
      
      x$gini.outlier_sd = x$gini > max(x$gini_mode + cutoff_value *sd_est, quantile(x$gini, (1-cutoff_percentile)))
      
      # gini boxplot outlier
      Q_3 = quantile(x$gini, 0.75)
      Q_1 = quantile(x$gini, 0.25)
      IQR = Q_3 - Q_1
      upper_whisker = min(max(x$gini), Q_3 + factor_iqr * IQR)
      
      x$gini.outlier = x$gini > upper_whisker
      
      #ggplot(data=x) + geom_quasirandom(aes(x="NA", y=gini, color=gini.outlier))
      
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}

# compute per dataset alpha value cutoff (we use hmm.alpha by default, can be overriden by specifying alpha.select)
qc_alpha <- function(dft, factor_iqr=3, cutoff_percentile=0.05, cutoff_value=2.0, hard_cutoff=0.05){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("hmm.alpha" %in% colnames(dft))
  
  # MAPD analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      
      if(!("alpha.select" %in% colnames(x))){
        x$alpha.select = x$hmm.alpha
      }
      # DEBUG
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A73044A"))
      # x = df2 %>% dplyr::filter(SLX %in% c("SLX-A90553C"))
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A96139A"))
      Q_3 = quantile(x$alpha.select, 0.75, na.rm=TRUE)
      Q_1 = quantile(x$alpha.select, 0.25, na.rm=TRUE)
      IQR = Q_3 - Q_1
      upper_whisker = min(max(x$alpha.select), Q_3 + factor_iqr * IQR)
      lower_whisker = max(min(x$alpha.select), Q_1 - factor_iqr * IQR) 
     
      x$alpha.outlier_percentile = x$alpha.select > quantile(x$alpha.select, cutoff_percentile) | x$alpha.select >  hard_cutoff
      # log normal distribution
      x$alpha.outlier = log(x$alpha.select) > mean(log(x$alpha.select)) + (cutoff_value * sd(log(x$alpha.select)))  | x$alpha.select >  hard_cutoff
      x$alpha.outlier_iqr = x$alpha.select > upper_whisker  | x$alpha.select >  hard_cutoff
      
      #ggplot(data=x) + geom_histogram(aes(x=l2, color=l2.outlier), bins=100)
      #ggplot(data=x) + geom_quasirandom(aes(x="NA", y=alpha, color=alpha.outlier))
      
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}


# compute per dataset rpc cutoff
qc_error <- function(dft, factor_iqr=3, cutoff_percentile=0.05, cutoff_value=2.0, symmetric=FALSE){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("l2" %in% colnames(dft))
  
  # MAPD analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      
      # DEBUG
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A73044A"))
      # x = df2 %>% dplyr::filter(SLX %in% c("SLX-A90553C"))
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A96139A"))
      Q_3 = quantile(x$l2, 0.75)
      Q_1 = quantile(x$l2, 0.25)
      IQR = Q_3 - Q_1
      upper_whisker = min(max(x$l2), Q_3 + factor_iqr * IQR)
      lower_whisker = max(min(x$l2), Q_1 - factor_iqr * IQR) 
     
      if(symmetric){
        x$error.outlier_percentile = x$l2 > quantile(x$l2 , 1-cutoff_percentile) | x$l2 < quantile(x$l2, cutoff_percentile)
        x$error.outlier_sd = abs(x$l2 - mean(x$l2)) > cutoff_value * sd(x$l2)
        x$error.outlier = x$l2 > upper_whisker | x$l2 < lower_whisker
      }else{
        x$error.outlier_percentile = x$l2 < quantile(x$l2, cutoff_percentile)
        x$error.outlier_sd = x$l2 < mean(x$l2) - (cutoff_value * sd(x$l2))
        x$error.outlier = x$l2 > upper_whisker
      }
      
      #ggplot(data=x) + geom_histogram(aes(x=l2, color=l2.outlier), bins=100)
      #ggplot(data=x) + geom_quasirandom(aes(x="NA", y=l2, color=l2.outlier))
      
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}


# compute per dataset rpc cutoff
qc_rpc <- function(dft, cutoff_percentile=0.05, cutoff_value=2.0, symmetric=TRUE){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("rpc" %in% colnames(dft))
  
  # MAPD analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      
      # DEBUG
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A73044A"))
      # x = df2 %>% dplyr::filter(SLX %in% c("SLX-A90553C"))
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A96139A"))
     
      if(symmetric){
        x$rpc.outlier_percentile = x$rpc > quantile(x$rpc, 1-cutoff_percentile) | x$rpc < quantile(x$rpc, cutoff_percentile)
        x$rpc.outlier = abs(x$rpc - mean(x$rpc)) > cutoff_value * sd(x$rpc)
      }else{
        x$rpc.outlier_percentile = x$rpc < quantile(x$rpc, cutoff_percentile)
        x$rpc.outlier = x$rpc < mean(x$rpc) - (cutoff_value * sd(x$rpc))
      }
      
      #ggplot(data=x) + geom_histogram(aes(x=rpc, color=rpc.outlier), bins=100)
      #ggplot(data=x) + geom_quasirandom(aes(x="NA", y=rpc, color=rpc.outlier))
      
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}

# compute per dataset prediction quality cutoff
qc_predict <- function(dft, cutoff_percentile=0.05, cutoff_value=2.0){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("prediction.quality" %in% colnames(dft))
  
  # MAPD analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      
      # DEBUG
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A73044A"))
      # x = df2 %>% dplyr::filter(SLX %in% c("SLX-A90553C"))
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A96139A"))
      
      x$predict.outlier_percentile = x$prediction.quality < quantile(x$prediction.quality, cutoff_percentile, na.rm=TRUE)
      x$predict.outlier = x$prediction.quality < mean(x$prediction.quality, na.rm=TRUE) - (cutoff_value * sd(x$prediction.quality, na.rm=TRUE))
      
      #ggplot(data=x) + geom_histogram(aes(x=prediction.quality, color=predict.outlier), bins=100)
      #ggplot(data=x) + geom_quasirandom(aes(x="NA", y=rpc, color=rpc.outlier))
      
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}


# generic kernel density based filter function
kde2d_filter <- function(x, y, q = 0.05, n = 200, plt = FALSE) {
  stopifnot(length(x) == length(y))
  z <- MASS::kde2d(x, y, n = n)
  cuts <- data.frame(x = cut(x, seq(z$x[1], z$x[n], length.out = n + 1)), 
                     y = cut(y, seq(z$y[1], z$y[n], length.out = n + 1)), 
                     point_number = seq_along(x))
  dens <- expand.grid(x = levels(cuts$x), y = levels(cuts$y))
  dens$z <- as.vector(z$z)
  dens <- merge(cuts, dens, sort = FALSE, all.x = TRUE)
  dens$z[is.na(dens$z)] <- 0
  pts <- dens$point_number[dens$z >= quantile(dens$z, prob = q)]
  
  if (plt) {
    plot(x, y)
    points(x[pts], y[pts], pch = 20, col = 'red')
    contour(z, drawlabels = FALSE, add = TRUE, lwd = 2, nlevels = 10, 
            col = rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
  }
  
  #invisible(pts)
  res=rep(TRUE, length(x))
  res[pts] = FALSE
  invisible(res)
}

# compute per dataset mapd cutoff (using dmapd value, i.e. mapd corrected for rpc)
qc_mapd <- function(dft, cutoff_percentile=0.05, cutoff_value=2.0, symmetric=FALSE, densityControl=FALSE, densityCutoff=0.05){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("mapd" %in% colnames(dft))
  stopifnot("rpc" %in% colnames(dft))
  
  # MAPD analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      #model = stats::loess(mapd ~ rpc, data=x, span=0.3)
      #pred = predict(model, newdata=dplyr::tibble(rpc = x$rpc))
      #x$dmapd = x$mapd - pred
      
      # DEBUG
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A73044A"))
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A90553C"))
      #x = df2 %>% dplyr::filter(SLX %in% c("SLX-A96139A"))
      
      outlier_mapd = kde2d_filter(x$rpc, x$mapd, q=densityCutoff, plt=FALSE)
      model = robustbase::lmrob(mapd ~ I(1/rpc), data=x)
      #pred = predict(model, newdata=dplyr::tibble(rpc=x$rpc))
      x$dmapd.rweights = model$rweights
      x$dmapd.residuals = residuals(model)
      x$dmapd = model$residuals / model$scale
      
      #ggplot(data=x) + geom_point(aes(x=rpc, y=mapd)) + 
      #  geom_point(aes(x=rpc, y=pred), color="red")
      #
      #ggplot(data=x) + geom_point(aes(x=rpc, y=dmapd, color=mapd > 0.3)) + 
      #  geom_point(aes(x=rpc, y=pred), color="red")
      
      if(symmetric){
        x$dmapd.outlier_percentile = abs(x$dmapd) > quantile(abs(x$dmapd), 1-cutoff_percentile)
        x$dmapd.outlier = abs(x$dmapd) > cutoff_value
      }else{
        x$dmapd.outlier_percentile = x$dmapd > quantile(x$dmapd, 1-cutoff_percentile)
        x$dmapd.outlier = x$dmapd > cutoff_value
      }
      
      if(densityControl){
        x$dmapd.outlier = x$dmapd.outlier | outlier_mapd
      }
      
      #ggplot(data=x) + geom_point(aes(x=rpc, y=dmapd, color=dmapd.outlier)) + 
      #  geom_point(aes(x=rpc, y=pred), color="red")

      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}


# compute per dataset gini cutoff (using dgini value, i.e. gini corrected for rpc)
qc_gini_norm <- function(dft, cutoff_percentile=0.05, cutoff_value=2.0, symmetric=FALSE, densityControl=FALSE, densityCutoff=0.05){
  
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("gini_normalized" %in% colnames(dft))
  stopifnot("rpc" %in% colnames(dft))
  
  # Gini coefficient analysis within each dataset
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>%
    dplyr::group_modify(function(x, y){

      # DEBUG
#      x = df3 %>% dplyr::filter(UID %in% c("UID-DLP-SA1044", "UID-DLP-SA928"), cellcycle %in% c("G1"))
#      x = df3 %>% dplyr::filter(UID %in% c("UID-DLP-SA1044"), cellcycle %in% c("G1"))
#      ggplot(data=x) + geom_point(aes(x=sqrt(rpc), y=gini_normalized, color=outlier_gini)) +
#        geom_smooth(aes(x=sqrt(rpc), y=gini_normalized), method=robustbase::lmrob)
      outlier_gini = kde2d_filter(x$rpc, x$gini_normalized, q=densityCutoff, plt=FALSE)
      
      model = robustbase::lmrob(gini_normalized ~ I(rpc), data=x)
      #summary(model)
      pred = predict(model, newdata=dplyr::tibble(rpc=x$rpc))
      x$dgini.rweights = model$rweights
      x$dgini.residuals = model$residuals
      x$dgini= model$residuals / model$scale
      
      if(symmetric){
        x$dgini.outlier_percentile = abs(x$dgini) > quantile(abs(x$dgini), 1-cutoff_percentile)
        x$dgini.outlier = abs(x$dgini) > cutoff_value
      }else{
        x$dgini.outlier_percentile = x$dgini > quantile(x$dgini, 1-cutoff_percentile)
        x$dgini.outlier = x$dgini > cutoff_value
      }
      
      if(densityControl){
        x$dgini.outlier = x$dgini.outlier | outlier_gini
      }
      
      #ggplot(data=x) + geom_point(aes(x=rpc, y=dgini, color=dgini.outlier)) + 
      #  geom_point(aes(x=rpc, y=pred), color="red")

      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}

# simple script to identify cells in dataset with branch length longer than threshold
get_branch_length <- function(tree, threshold=30){
  require(dendextend)
  #require(phylogram)
  
  dend = as.dendrogram(tree)
  
  subtrees <- dendextend::partition_leaves(dend)
  leaves <- subtrees[[1]]   # assume top node is used by all sub trees
  
  pathRoutes <- function(leaf) {
    which(sapply(subtrees, function(x) leaf %in% x))
  }
  
  paths <- lapply(leaves, pathRoutes)
  heights <- get_nodes_attr(dend, "height")
  leaf_height = unlist(lapply(paths, function(x){
    height = min(setdiff(heights[x], 0))
    return(height)
  }))
  names(leaf_height) = labels(dend)
  
  return(leaf_height)
}

#' binsToUseInternal
#'
#' \code{binsToUseInternal} modified method from QDNAseq package to select bins to use
binsToUseInternal <- function(object){
  require(Biobase, warn.conflicts = FALSE)
  if ("use" %in% colnames(Biobase::fData(object))) {
    return(Biobase::fData(object)$use)
  }else{
    warning("No filter available for this object.")
  }
  
  return(rep(TRUE, times=nrow(object)))
}

#' copyNumberDistance
#'
#' \code{copyNumberDistance} wrapper around MEDICC copy number distance as implemented by 
#' Zeira, Ron, and Benjamin J. Raphael. "Copy number evolution with weighted aberrations in cancer."
# Bioinformatics 36.Supplement_1 (2020): i344-i352.
copyNumberDistance <- function(input, weighted_distance=FALSE, BASEDIR=BASEDIR){
  
  require(QDNAseq, quietly=TRUE)
  require(reticulate, quietly=TRUE)
  reticulate::source_python(file.path(BASEDIR, "scUnique/R/WeightedCopyNumberDistanceFunctions.py"), convert=TRUE)
  # reticulate::source_python("~/scUnique/R/WeightedCopyNumberDistanceFunctions.py", convert=TRUE)
  
  if("QDNAseqCopyNumbers" %in% class(input)){
    stopifnot(!is.null(assayDataElement(input, "copynumber")))
    
    mat = Biobase::assayDataElement(input, "copynumber")
    is_valid = apply(mat, 1, function(x) sum(is.na(x))) == 0
    mat = t(unname(round(mat[is_valid, ], digits = 0)))
    n_cells = dim(mat)[[1]]
  }else{
    stopifnot(all(round(input, digits = 0) == input | is.na(input)))
    stopifnot(dim(input)[[1]] < dim(input)[[2]])
    is_valid = apply(input, 2, function(x) sum(is.na(x))) == 0
    mat = input[, is_valid]
    n_cells = dim(mat)[[1]]
  }
  rownames(mat) = NULL
  colnames(mat) = NULL
  
  # NOTE, run_CopyNumberDistance is python interface function
  d = run_CopyNumberDistance(mat=mat, weighted=FALSE)
  cn_dist = as.dist(d)
  stopifnot(dim(cn_dist)[[1]] == n_cells)
  
  return(cn_dist)
}


## Create a segmentation from unique events (in copynumber slot) ====
clusterCopynumberProfile <- function(object, cluster_profiles, membership, unique_events){
  
  n_cells = dim(object)[[2]]
  copynumber_slot = matrix(NA, nrow=dim(object)[[1]], ncol=dim(object)[[2]])
  segmentation_slot = matrix(NA, nrow=dim(object)[[1]], ncol=dim(object)[[2]])
  # cluster_profiles
  # membership
  # unique_events
  
  for(cellindex in 1:n_cells){
    cellname = colnames(object)[[cellindex]]
    
    kic = which(unlist(lapply(membership, function(x) cellindex %in% x)))
    copynumber_slot[, cellindex] = cluster_profiles[, kic]
    segmentation_slot[,cellindex] = cluster_profiles[, kic]
    gr = unique_events[unique_events$name == cellname]
    n_len = length(gr)
    if(n_len == 0){
      next;
    }
    for(ni in 1:n_len){
      g = gr[ni]
      values = copynumber_slot[g$start.index:g$end.index, cellindex]
      values[!is.na(values)] = g$h1
      copynumber_slot[g$start.index:g$end.index, cellindex] = values
      # segmentation_slot[g$start.index:g$end.index, cellindex] = values
    }
  }
  
  ## Fully segmented object
  clusterCN = new("QDNAseqCopyNumbers",
                bins=Biobase::featureData(object),
                copynumber=copynumber_slot,
                phenodata=Biobase::pData(object))
  Biobase::assayDataElement(clusterCN, "calls") <- Biobase::assayDataElement(object, "calls")
  Biobase::assayDataElement(clusterCN, "probdloss") <- Biobase::assayDataElement(object, "probdloss")
  Biobase::assayDataElement(clusterCN, "segmented") <- segmentation_slot
  
  return(clusterCN)
}

# get a medicc compatible data frame from QDNAseq object (copynumber is based on copynumber slot)
getMediccTable <- function(object, mediccSizeCutoff=3, mediccDiffAmplitude=1, valueField="copynumber"){
  require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE)
  
  CN = Biobase::assayDataElement(object, valueField)
  valid = binsToUseInternal(object)
  
  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)
  
  # remove small and small amplitude segments
  CN_clean = apply(CN[valid, ], 2, function(x){
    rl2 = Rle(x)
    candidates = rl2@lengths <= mediccSizeCutoff
    candidates_amplitude = c(FALSE, abs(diff(rl2@values)) <= mediccDiffAmplitude)
    candidates_replacement = c(0, rl2@values)
    rl2@values = ifelse(candidates & candidates_amplitude, candidates_replacement, rl2@values)
    x2 = as.vector(rl2)
    
    return(x2)
    })
  CN[valid, ] = CN_clean
  
  ids = paste0("chr",chromosome,":", apply(CN, 1, function(x){
    v=paste(x, sep="-", collapse="-"); return(v)}))
  rl = Rle(ids[valid])
  
  export_single_sample <- function(sample_index){
    cn = CN[valid, sample_index]
    sample = colnames(object)[[sample_index]]
    
    cn_values = cn[start(rl)]
    chromosomeValues = chromosome[valid][start(rl)]
    
    # MEDICC2 follows the BED convention for segment coordinates, i.e. segment start is at 0 and the segment end is non-inclusive.
    gr = GenomicRanges::GRanges(seqnames = chromosomeValues, IRanges(start(rl),
                                                      width=runLength(rl)+1, sample=sample, copynumber = cn_values))
    
    df = data.frame(gr) %>% dplyr::rename("chrom"="seqnames", "sample_id"="sample", "cn_a"=copynumber) %>% dplyr::select(-strand, -width) %>%
      dplyr::relocate(sample_id, chrom, start, end, cn_a)
    
    return(df)
  }
  
  dfs = dplyr::bind_rows(lapply(1:ncol(object), export_single_sample))
  dfs$chrom = ifelse(dfs$chrom == "X", "23", dfs$chrom)
  dfs$chrom = ifelse(dfs$chrom == "Y", "24", dfs$chrom)
  
  return(dfs)
}


# get a SITKA compatible data frame from QDNAseq object (copynumber is based on copynumber slot)
getSITKATable <- function(object, value="copynumber"){
  require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE)
  
  CN = Biobase::assayDataElement(object, value)
  valid = binsToUseInternal(object)
  
  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)
  
  export_single_sample_sitka <- function(sample_index){
    
    cn = CN[valid, sample_index]
    sample = colnames(object)[[sample_index]]
    
    ids = paste0("chr",chromosome[valid],":", sapply(cn, function(x){
      v=paste(x, sep="-", collapse="-"); return(v)}))
    rl = Rle(ids)
    
    l = length(rl@values)
    cn_values = cn[start(rl)]
    chromosomeValues = chromosome[valid][start(rl)][2:l]
    start_values = start(rl)[2:l]
    
    gr = GenomicRanges::GRanges(seqnames = chromosomeValues, IRanges(start_values-1,
                                                                     width=2, sample=sample,
                                                                     copynumber_l = cn_values[1:(l-1)],
                                                                     copynumber_r = cn_values[2:l]))
    
    df = data.frame(gr) %>% dplyr::rename("chrom"="seqnames", "sample_id"="sample", "cn_l"=copynumber_l, "cn_r"=copynumber_r) %>% dplyr::select(-strand, -width) %>%
      dplyr::relocate(sample_id, chrom, start, end, cn_l, cn_r)
    
    return(df)
  }
  
  dfs = dplyr::bind_rows(lapply(1:ncol(object), export_single_sample_sitka))
  
  return(dfs %>% dplyr::mutate(loci=paste0(chrom, "_", start, "_", end), tipInclusionProbabilities=1.0) %>% 
           dplyr::rename("cells"="sample_id") %>% dplyr::select(cells, loci, tipInclusionProbabilities))
}

# only apply confirmed events on profiles unique:
finalizeCopynumberProfile = function(object, profiles_background, grl_unique){
  
  profiles_unique = profiles_background
  finalUnique = matrix(-1, nrow=nrow(profiles_unique), ncol=ncol(profiles_unique))
  events = unlist(grl_unique)
  valid = binsToUseInternal(object)
  
  # finalUnique is a matrix displaying ONLY unique events (baseline is -1)
  finalUnique[!valid, ] = NA
  rownames(finalUnique) = rownames(profiles_background)
  colnames(finalUnique) = colnames(profiles_background)
  for (ei in 1:length(events)){
    event = events[ei,]
    profiles_unique[,event$name][event$start.index:event$end.index] = event$h1
    finalUnique[,event$name][event$start.index:event$end.index] = event$h1
  }
  profiles_unique[!valid,] = NA
  ## update objects with final copy number profile
  finalCN = Biobase::assayDataElementReplace(object, "copynumber", profiles_unique)
  finalCN = Biobase::assayDataElementReplace(finalCN, "segmented", profiles_background)
  return(list("finalCN"=finalCN, "profiles_unique"=profiles_unique, "finalUnique"=finalUnique))
}


## algorithm to filter highly segmented cells from sample
identify_oversegmented = function(object){
  
  require(S4Vectors)
  valid = binsToUseInternal(object)
  n_cells = dim(object)[[2]]
  n_bins = dim(object)[[1]]
  cn = object@assayData$copynumber[valid,]
  breakpoints = apply(cn, 2, function(x) diff(x) != 0)
  n_breakpoints = apply(breakpoints, 2, function(x) sum(x, na.rm=TRUE))
  
  get_individual_breakpoints = function(index){
    
    x = breakpoints[, index, drop=FALSE]
    y = breakpoints[, -index]
    
    y_margin = apply(y, 1, sum)
    individual_breakpoints = sum(x & y_margin == 0)
    
    return(individual_breakpoints)
  }
  
  n_individual_breakpoints = sapply(seq(1, n_cells), get_individual_breakpoints)
  names(n_individual_breakpoints) = colnames(object)
  
  # distance between breakpoints
  breakpoint_distance = apply(cn, 2, function(x){
    rl=Rle(x);
    return(list(mean(rl@lengths), sd(rl@lengths)));
    })
  mean_breakpoint_distance = unlist(lapply(breakpoint_distance, "[[", 1))
  sd_breakpoint_distance = unlist(lapply(breakpoint_distance, "[[", 2))
  
  # distance between breakpoints, by chromosomes
  breakpoint_distance_chromosomes = apply(cn, 2, function(x){
    na=names(x)
    chr=unlist(lapply(str_split(na, ":"), "[[", 1))
    rl=Rle(x);
    bps=rep(0, length(x))
    bps[cumsum(rl@lengths)] = 1
    bps_per_chrom = lapply(split(bps, chr), sum)
    return(list(mean(unlist(bps_per_chrom)), sd(unlist(bps_per_chrom))));
  })
  mean_bps_per_chromosome = unlist(lapply(breakpoint_distance_chromosomes, "[[", 1))
  sd_bps_per_chromosome = unlist(lapply(breakpoint_distance_chromosomes, "[[", 2))
  

  return(list(n_individual_breakpoints=n_individual_breakpoints, n_breakpoints=n_breakpoints,
              mean_breakpoint_distance=mean_breakpoint_distance, sd_breakpoint_distance=sd_breakpoint_distance,
              mean_bps_per_chromosome=mean_bps_per_chromosome, sd_bps_per_chromosome=sd_bps_per_chromosome))
}

## assign cells (QDNAseq object) to clones (list of copy number profiles), based on MEDICC distance
assign_CellsToClones_MEDICC1 <- function(cells, clones, ncores=8, BASEDIR=BASEDIR){
  require(parallel)
  require(doParallel)
  require(foreach)
  require(nnet)
  valid = binsToUseInternal(cells)
  cnames = colnames(cells)
  #mat = matrix(data=NA, nrow=dim(cells)[[2]], ncol=length(clones))
  
  clone_names = names(clones)
  n_cells = ncol(cells)
  pb <- txtProgressBar(min = 1, max = n_cells, initial = 0, width=100, style = 3, file=stdout())
  cl <- parallel::makeCluster(ncores, outfile="/dev/null")
  doParallel::registerDoParallel(cl)
  mat = foreach::foreach(index1 = 1:n_cells, .combine='rbind', .inorder=TRUE, .multicombine=TRUE,
                   .packages=c("reticulate", "tensorflow", "QDNAseq", "Biobase")) %dopar% {
                     
                     
    if(index1 %% 100 == 1){
      setTxtProgressBar(pb, index1)
    }
    e = base::list2env(list("BASEDIR"=BASEDIR, "cells"=cells, "valid"=valid))
    source(file.path(e$BASEDIR, "scUnique/R/core.R"))
    ma = matrix(data=NA, nrow=1, ncol=length(clones))
    colnames(ma) = clone_names
    
      #print(index1)
      repCellProfile = Biobase::assayDataElement(cells[e$valid, index1], "copynumber")
    for(index2 in 1:length(clones)){
      #print(index2)
      cloneProfile = clones[[index2]] 
      stopifnot(length(repCellProfile[,1]) == length(cloneProfile))
      mt = cbind(repCellProfile, cloneProfile)
      dst = copyNumberDistance(t(mt), weighted_distance = TRUE, BASEDIR=e$BASEDIR)
      ma[1, index2] = dst
    }
  return(ma)
  }
  close(pb)
  parallel::stopCluster(cl)
  
  # note ties in min are broken at random (nnet package - which.is.max)
  exclude_cols = colnames(mat) %in% c("0", "-1")
  dft = dplyr::tibble(name = colnames(cells), clone=clone_names[!exclude_cols][apply(mat[, !exclude_cols], 1, function(x) nnet::which.is.max(-x))],
                      clone_neutral=clone_names[apply(mat, 1, function(x) nnet::which.is.max(-x))], dist=apply(mat, 1, min))
  
  return(list("assignment"=dft, "distances"=mat))
}


convert_cProfiles_to_cnProfiles <- function(cProfiles, object){
  cpr = do.call("rbind", cProfiles)
  valid = binsToUseInternal(object)
  newCopynumber = matrix(data=NA, nrow=dim(object)[[1]], ncol=length(cProfiles))
  newCopynumber[valid, ] = t(cpr)
  
  df <- data.frame(name=names(cProfiles))
  metaData <- data.frame(labelDescription=c("name"))
  
  newCN = new("QDNAseqCopyNumbers",
              bins=featureData(object),
              copynumber=newCopynumber,
              phenodata=AnnotatedDataFrame(data=df, varMetadata=metaData))
  
  Biobase::assayDataElement(newCN, "segmented") <- Biobase::assayDataElement(newCN, "copynumber")
  colnames(newCN) = names(cProfiles)

  return(newCN)
}


assign_CellsToClones_MEDICC2 <- function(clones, cells=NULL, ncores=8){
  
  require(tensorflow)
  require(S4Vectors)
  require(zoo)
  require(stringr)
  require(proxy)
  require(RColorBrewer)
  
  require(foreach)
  require(parallel)
  require(doParallel)
  
  stopifnot("QDNAseqCopyNumbers" %in% class(clones))
  n_cells = ifelse(is.null(cells), 1, dim(cells)[[2]])
  if(!is.null(cells)){
    description = Biobase::pData(cells) %>% dplyr::select(name)
    Biobase::pData(cells) = description
    pb <- txtProgressBar(min = 1, max = n_cells, initial = 0, width=100, style = 3, file=stdout())
  }
  
  require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
  reticulate::py_discover_config()
  e = new.env()
  reticulate::source_python(file.path("~/scUnique/interface_medicc2.py"), convert=TRUE, envir=e)
  
  
  start_time_medicc <- Sys.time()
  dist_mat = dplyr::tibble()
  for(cellindex in 1:n_cells){
    
    if(!is.null(cells)){
      object = combineQDNASets(list(clones, cells[, cellindex]))
      cellname = colnames(cells)[[cellindex]]
    }else{
      object = clones
    }

  df_medicc = getMediccTable(object) # NOTE, by default accessing copynumber slot!
  start_time_medicc <- Sys.time()
  out = reticulate::py_capture_output({
  results = e$interface_medicc(df_medicc %>% dplyr::select(sample_id, chrom, start, end, cn_a), as.integer(ceiling(ncores*(2/3))))
  })
  
  pdms = results[[5]];
  
  if(is.null(cells)){
    return(list("distances"=pdms))
  }
  
  if(cellindex %% 25 == 1){
    setTxtProgressBar(pb, cellindex)
  }
  
  
  vec = pdms[[cellname]]
  names(vec) = paste0("C", rownames(pdms))
  df0 = as_tibble(t(vec))
  df0[[paste0("C", cellname)]] = NULL
  dist_mat = dplyr::bind_rows(dist_mat,
         dplyr::bind_cols(dplyr::tibble(name=cellname), df0))
  }
  end_time_medicc <- Sys.time()
  print(paste0("Run MEDICC algorithm ",  difftime(end_time_medicc,start_time_medicc,units="mins")))
  
  assignment_long = dist_mat %>% tidyr::pivot_longer(-name, names_to="clone", values_to="distance")
  assignment = assignment_long %>% dplyr::mutate(ties = sample(n())) %>% 
    dplyr::group_by(name) %>%
    dplyr::arrange(distance, ties, .by_group = TRUE) %>%
    dplyr::slice(1)
    
  return(list("assignment"=assignment, "distances"=dist_mat, "assignment_long"=assignment_long))
}


## obtain clonal profiles from ctree (cluster assignment, cutree), from dendrogram (dt), and from medicc profiles
get_ClonalProfiles <- function(ctree, dendt, cprofiles){
  result = list()
  cnodes = unlist(lapply(cprofiles, `[[`, 1))
  for(ct in sort(unique(ctree))){
    if(ct == 0){
      next
    }
    leaves = labels(dendt)[ctree == ct]
    phy = phylogram::as.phylo.dendrogram(dendt)
    
    l = c()
    for(i in 1:1000){
      le = sample(leaves, size=0.5*length(leaves), replace=FALSE)
      mrca_node = ape::getMRCA(phy, le)
      mrca_nodelabel = phy$node.label[mrca_node-Ntip(phy)]
      l = c(l, mrca_nodelabel)
    }
    tbl = table(l)
    if(length(tbl) != 1){
      warning("Non-overlapping clusters")
      print(tbl)
    }
    mrca_consensus = names(tbl[tbl == min(tbl)])
    entry = cprofiles[[which(mrca_consensus == cnodes)]]
    
    #mrca_node = ape::getMRCA(phy, leaves)
    #mrca_nodelabel = phy$node.label[mrca_node-Ntip(phy)]
    #entry = cprofiles[[which(mrca_nodelabel == cnodes)]]
    
    result[[ct]] = entry[[4]]
  }
  return(result)
}

# get most recent common ancestor from dendrogram and medicc profiles based on selection of leaves
get_MRCA <- function(dendt, cprofiles, leaves){
  
  cnodes = unlist(lapply(cprofiles, `[[`, 1))
  #leaves = labels(dendt) 
  
  phy = phylogram::as.phylo.dendrogram(dendt)
  mrca_node = ape::getMRCA(phy, leaves)
  mrca_nodelabel = phy$node.label[mrca_node-Ntip(phy)]
  
  entry = cprofiles[[which(mrca_nodelabel == cnodes)]]
    
  return(entry[[4]])
}

## core methods for copynumber analysis
# Copyright 2022, Michael Schneider, All rights reserved.

#' binsToUseInternal
#'
#' \code{binsToUseInternal} modified method from QDNAseq package to select bins to use
binsToUseInternal <- function(object){
  require(Biobase, warn.conflicts = FALSE)
  if ("use" %in% colnames(Biobase::fData(object))) {
    return(Biobase::fData(object)$use)
  }else{
    warning("No filter available for this object.")
  }

  return(rep(TRUE, times=nrow(object)))
}

#' getGenomeInformation
#'
#' \code{getGenomeInformation} access QDNAseq object genome information
getGenomeInformation <- function(object){
  return(Biobase::experimentData(object)@other)
}

#' createGR 
#'
#' \code{createGR} create a GRanges representation 
#'
#' @param object QDNAseq object
#' @param value character add QNDAseq assayData element of object to the GRanges object as a metadata column [copynumber, calls, segmented]
#' @param valid index include only elements at valid positions
#' @param asRleList convert GRanges output to RleList
#' 
#' @return GRanges object (or RleList, if asRleList is set)
createGR <- function(object, value=NULL, valid=NULL, asRleList=FALSE, variable=NULL){
  require(GenomicRanges, warn.conflicts = FALSE)
  stopifnot(dim(object)[[2]] == 1)

  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)

  tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
  if(asRleList){
    chunkStart = do.call("c", lapply(rle(seqName)$lengths, function(x) 1:x))
    chunkEnd = do.call("c", lapply(rle(seqName)$lengths, function(x) 1:x))
  }else{
    chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
    chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
  }
  if(!is.null(variable)){
    gr = GRanges(chromosome, IRanges(chunkStart, chunkEnd), annotation=variable)
  }else{
    gr = GRanges(chromosome, IRanges(chunkStart, chunkEnd))
  }
  

  if(!is.null(value)){
    val = Biobase::assayDataElement(object, value)
    if(!is.null(valid)){
      val[!valid] = NA
    }

    df = data.frame(val=as.numeric(val))
    colnames(df) = value
    mcols(gr) = df

    if(asRleList){
      gr = GenomicRanges::mcolAsRleList(gr, value)
    }
  }else{
    if(asRleList){
      stop("asRleList not supported if value is not set")
    }
    if(!is.null(valid)){
      gr = gr[valid]
    }
  }

  return(gr)
}

#' retrieveAssayData 
#'
#' \code{retrieveAssayData} extract assayData from QDNAseq object
#'
#' @param object QDNAseq object
#' @param index sample index into QDNAseq object
#' @param value character QNDAseq assayData element of object to the GRanges object as a metadata column [copynumber, calls, segmented]
#' @param valid index include only elements at valid positions
#' @param asRleList convert GRanges output to RleList
#' 
#' @return RleList with value as metadata
#' 
#' @export
retrieveAssayData <- function(object, index, value, valid=NULL){
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  stopifnot(length(index) == 1)
  return(createGR(object[,index,drop=FALSE], value=value, valid=valid, asRleList=TRUE))
}

#' retrieveSegmentation 
#'
#' \code{retrieveSegmentation} create a GRanges representation 
#'
#' @param object QDNAseq object
#' @param index sample index into QDNAseq object
#' @param valid index include only elements at valid positions
#' @param formatIR return IRanges instead of GRanges
#' 
#' @return return GRanges object (or IRanges if flag is set) containing segmentation
#' 
#' NOTE
#' segment encoding is processed as follows:
#   1. Remove NaNs that are covered by two equal segments (30 NA 30) -> 30
#   2. Keep NaNs that border two unequal segments (30 NA 31) -> 30 NA 31
#   3. Remove NaNs at start or end (NA 30) -> 30, (30 NA) -> 30
#   4. Do this on a per chromosome level (chromosome encoding in chromosome variable)
#   5. Keep NaN elements in segmented
retrieveSegmentation <- function(object, index, valid=NULL, formatIR=FALSE, slot="segmented"){

  if(class(object)[[1]] != "QDNAseqCopyNumbers"){
    segments = object
    stopifnot(!is.null(rownames(object)))
    stopifnot(is.null(valid))
  }else{
    stopifnot(!is.null(Biobase::assayDataElement(object, slot)))
    if(is.null(valid)){
      segments = Biobase::assayDataElement(object, slot)[,index,drop=FALSE]
    }else{
      segments = Biobase::assayDataElement(object, slot)[valid,index,drop=FALSE]
    }
  }
  
  replacena <- function(rl) {
    lv = rl@values
    ll = rl@lengths
    stopifnot(is.numeric(lv))
    indx <- is.na(lv)
    mat = matrix(c(lv, ll), nrow=2, byrow = TRUE)
    
    if(all(is.na(lv))){
      return(list(lv=lv, ll=ll))
    }
    
    for(x in 1:length(lv)){
      if(is.na(lv[x])){
        
        # left border
        if(x == 1){
          mat[1,x+1] = mat[1,2]
          mat[2,x+1] = mat[2,1] + mat[2,2]
          # mat[1,x] = NA
          mat[2,x] = 0
          next
        }
        
        # right border
        if(x == length(lv)){
          mat[1,x-1] = mat[1,x-1]
          mat[2,x-1] = mat[2,x-1] + mat[2,x]
          # mat[1,x] = NA
          mat[2,x] = 0
          next
        }
        
        ## normal case
        left = x-1
        right = x+1
        
        if(mat[1,left] == mat[1,right]){
          mat[1,x] = mat[1,left]
          mat[2,x] = mat[2,left] + mat[2,x] + mat[2,right]
          # mat[1,left] = NA
          mat[2,left] = 0
          # mat[1,right] = NA
          mat[2,right] = 0
        }
        
      }
    }
    
    lv = mat[1,]
    ll = mat[2,]
    return(list(lv=lv, ll=ll))
  }
  
  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)
  
  segments.cleaned = do.call("RleList", lapply(split(segments, chromosome), function(x){
    tmp = replacena(Rle(x))
    lv = tmp$lv
    ll = tmp$ll
    remove_idx = ll == 0
    lv = lv[!remove_idx]
    ll = ll[!remove_idx]
    vec = as.vector(Rle(values=lv, lengths = ll))
    stopifnot(length(vec) == length(x))
    ## NOTE: dummy encode NA as -1, to avoid filtering segments in bindAsGRanges -> replace values again after transform
    rl = Rle(as.numeric(vec))
    rl@values[is.na(rl@values)] = -1
    return(rl)
  }))
  
  # transform back to GRanges representation
  grs_new = bindAsGRanges(segmentation = segments.cleaned)
  grs_new$segmentation[grs_new$segmentation == -1] = NA
  # see note in function documentation
  
  stopifnot(sum(grs_new@ranges@width) == dim(object)[[1]])
  
  if(formatIR){
    n_bins = dim(object)[[1]]
    ir = grs_new@ranges
    ir_start = c(0, cumsum(ir@width))+1
    ir_end = cumsum(ir@width)
    ir = IRanges(start=ir_start[1:(length(ir_start)-1)], end=ir_end)
    return(ir)    
  }
  
  return(grs_new)
}

#' selectChromosomes 
#'
#' \code{selectChromosomes} select subregion of QDNAseq based on chromosome names either to be included or excluded
#'
#' @param object QDNAseq object
#' @param include character chromosome names (either in format "1" or "chr1") to include
#' @param exclude character chromosome names (either in format "1" or "chr1") to exclude
#' 
#' @return object containing only selected chromosomes
selectChromosomes <- function(object, include=NULL, exclude=NULL){
  require(S4Vectors, warn.conflicts = FALSE)
  
  if(!is.null(include) && !is.null(exclude)){
    stop("Choose either include or exclude regions.")
  }
  stopifnot(is.null(include) || all(is.character(include)))
  stopifnot(is.null(exclude) || all(is.character(exclude)))
  
  
  ploidyRegion=c(as.character(seq(1,22)), "X", "Y")
  if(!is.null(include)){
    if(all(base::startsWith(include, "chr"))){
      include = base::sub("^chr", "", include)
    }
    ploidyRegion = include
  }
  
  
  if(!is.null(exclude)){
    if(all(base::startsWith(exclude, "chr"))){
      exclude = base::sub("^chr", "", exclude)
    }
    ploidyRegion = ploidyRegion[!(ploidyRegion %in% exclude)]
  }
  
  gr = createGR(object[,1,drop=FALSE])
  indices = S4Vectors::decode(seqnames(gr) %in% ploidyRegion)
  
  return(object[indices, ,drop=FALSE])
}

#' smoothOperator
#' 
#' \code{smoothOperator} apply the summary function func to the data z. The summary function is computed 
#' for each segmented defined in the GRanges or IRanges object gr. Trim length reads are removed from the start and end of each segment
#'
#' @param z RleList or numeric vector (this will usually be the raw read count vector)
#' @param gr GRanges object with segmentation (z is RleList) or alternatively IRanges object (z is vector)
#' @param func function returning summary statistics for raw data (for a given segment)
#' @param minLength Numeric minimum number of bins to form a segment
#' @param trimLength Numeric number of bins to trim from start and end of segments
#' 
#' @return fit object containing best fit among the iterations from the SVI algorithm
smoothOperator <- function(z, gr, func, minLength=30, trimLength=0){
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
  require(IRanges, quietly=TRUE, warn.conflicts = FALSE)
  
  if(class(gr)[[1]] == "IRanges"){
    index = gr@width >= minLength
  }else{
    if(class(gr)[[1]] == "GRanges"){
      index = gr@ranges@width >= minLength
    }else{
      stop("RangeObject not supported.")
    }
  }
  
  vw = viewApply(Views(z, gr[index]), function(x) {
    
    # necessary to avoid XDouble errors
    x = as.numeric(x)
    
    y = x[(1+trimLength):(length(x)-trimLength)]
    y = func(y)
    
    return(y)
  }, simplify = TRUE)
  
  if(class(gr)[[1]] == "GRanges"){
    vw = do.call("c", lapply(names(vw), function(x) return(unlist(vw[[x]]))))
  }
  
  return(vw)
}

#' applyScale
#'
#' \code{applyScale} helper function to apply a given scaling factor to the counts in the form of a linear transformation
#'
#' @param CN QDNAseq object
#' @param scale Numeric Linear scaling factor, i.e. a in y = a x + b
#' @param shift Numeric Offest factor in linear transformation, i.e. b in y = a x + b
#'
#' @return scaled copy number object
applyScale <- function(CN, scale=1, shift=0){
  
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  result = new("QDNAseqCopyNumbers",
               bins=Biobase::featureData(CN),
               copynumber=Biobase::assayDataElement(CN, "copynumber"), #t((diag(as.numeric(unname(scale)), nrow=dim(CN)[[2]]) %*% t(Biobase::assayDataElement(CN, "copynumber"))) + shift),
               phenodata=Biobase::pData(CN))
  Biobase::assayDataElement(result, "segmented") <- t(diag(as.numeric(unname(scale)), nrow=dim(CN)[[2]]) %*% t(Biobase::assayDataElement(CN, "segmented"))) + shift
  if(!is.null(Biobase::assayDataElement(CN, "calls"))){
    Biobase::assayDataElement(result, "calls") <- Biobase::assayDataElement(CN, "calls")
  }
  if(!is.null(Biobase::assayDataElement(CN, "probdloss"))){
    Biobase::assayDataElement(result, "probdloss") <- Biobase::assayDataElement(CN, "probdloss")
  }
  Biobase::protocolData(result) <- Biobase::protocolData(CN)
  return(result)
}

#' computeRPC
#'
#' \code{computeRPC} compute average reads per copy per bin (rpc) value for a given QDNAseq object (based on calls and segments).
#' segments need to be a proxy for copy number state, i.e. a scaled data
#' 
#' @param object a QDNAseq object
#' @param rpc_value 
#' @param ploidy_value 
#' 
#' @return average reads per copy per bin
computeRPC <- function(object, rpc_value="rpc", ploidy_value="ploidy"){
  
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  
  valid = binsToUseInternal(object)
  rpc_raw = Biobase::assayDataElement(object, "calls")[valid,,drop=FALSE] / Biobase::assayDataElement(object, "segmented")[valid,,drop=FALSE]
  rpc_raw[is.infinite(rpc_raw)] <- NA
  rpc = apply(rpc_raw, 2, sum, na.rm=TRUE) / sum(valid, na.rm=TRUE)
  ploidy = mean(Biobase::assayDataElement(object, "segmented")[valid,,drop=FALSE], na.rm=TRUE)
  
  valid = binsToUseInternal(object) & !(base::startsWith(rownames(object), "X:") | base::startsWith(rownames(object), "Y:"))
  rpc_raw = Biobase::assayDataElement(object, "calls")[valid,,drop=FALSE] / Biobase::assayDataElement(object, "copynumber")[valid,,drop=FALSE]
  rpc_raw[is.infinite(rpc_raw)] <- NA
  rpc_final = apply(rpc_raw, 2, sum, na.rm=TRUE) / apply(rpc_raw, 2, function(x) sum(!is.na(x)))
  
  pd = Biobase::pData(object)
  pd[[rpc_value]] = rpc
  pd[[ploidy_value]] = ploidy
  if(!is.null(pd[["rpc"]])){
    pd[["rpc.old"]] = pd[["rpc"]]  
  }
  pd[["rpc"]] = rpc_final
  Biobase::pData(object) = pd
  
  return(object)
}

#' computeGini
#'
#' \code{computeGini} compute Gini coefficient for a given QDNAseq object (based on calls)
#' 
#' @param scaledCN a QDNAseq object
#' 
#' @return Gini coefficient
computeGini <- function(scaledCN){
  
  stopifnot(class(scaledCN)[[1]] == "QDNAseqCopyNumbers")
  Y = Biobase::assayDataElement(scaledCN, "calls")
  Y = Y[which(!is.na(Y[, 1])), ,drop=FALSE]
  
  igini <- function(x, weights=rep(1,length=length(x))){
    # source: https://rdrr.io/cran/reldist/src/R/gini.r
    ox <- order(x)
    x <- x[ox]
    weights <- weights[ox]/sum(weights)
    p <- cumsum(weights)
    nu <- cumsum(weights*x)
    n <- length(nu)
    nu <- nu / nu[n]
    sum(nu[-1]*p[-n]) - sum(nu[-n]*p[-1])
  }
  
  ngini <- function(i){
    obj = scaledCN[, i]
    stopifnot(dim(obj)[[2]] == 1)
    valid = binsToUseInternal(obj)
    Y = obj@assayData$calls[valid, ]
    C = obj@assayData$copynumber[valid, ]

    info = rownames(obj)
    seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)

    dft = dplyr::tibble(copy=C, count=Y, chromosome = chromosome[valid])
    dftp = dft %>% dplyr::group_by(chromosome, copy) %>% dplyr::summarise(n = n(), normalized_gini = igini(count), .groups="keep") %>%
      dplyr::ungroup() %>% dplyr::filter(copy > 0)
    
    va = dftp$n > 10
    normalized_gini = mean(rep(dftp$normalized_gini[va], times=dftp$n[va]), na.rm=TRUE)
    return(normalized_gini)
  }
  
  gini = apply(Y, 2, igini)
  gini_normalized = unlist(lapply(1:ncol(scaledCN), ngini))
  # scaledCN@phenoData@data$gini = gini
  pd = Biobase::pData(scaledCN)
  pd$gini = gini
  pd$gini_normalized = gini_normalized
  Biobase::pData(scaledCN) = pd
  return(scaledCN)
}

#' computeMAPD
#'
#' \code{computeMAPD} compute MAPD coefficient for a given QDNAseq object (based on calls)
#' 
#' @param scaledCN a QDNAseq object
#' 
#' @return MAPD coefficient
computeMAPD <- function(scaledCN){
  
  stopifnot(class(scaledCN)[[1]] == "QDNAseqCopyNumbers")
  Y = Biobase::assayDataElement(scaledCN, "calls")
  Y = Y[which(!is.na(Y[, 1])), ,drop=FALSE]
  
  mapd <- function(x){
    dd = (-1 * diff(x)) / mean(x)
    return(median( abs(dd - median(dd)) ))
  }
  
  mapd = apply(Y, 2, mapd)
  pd = Biobase::pData(scaledCN)
  pd$mapd = mapd
  Biobase::pData(scaledCN) = pd
  return(scaledCN)
}

#' computeL2
#'
#' \code{computeL2} compute l2 error for a given QDNAseq object (based on calls and scaling)
#' 
#' @param scaledCN a QDNAseq object
#' 
#' @return l2 error (normalized)
computeL2 <- function(scaledCN){
  
  stopifnot(class(scaledCN)[[1]] == "QDNAseqCopyNumbers")
  stopifnot("rpc" %in% colnames(Biobase::pData(scaledCN)))
  valid = binsToUseInternal(scaledCN)
  cnv = Biobase::assayDataElement(scaledCN, "copynumber")[valid, ,drop=FALSE]
  rpcs = Biobase::pData(scaledCN)[["rpc"]]
  error = (t(t(Biobase::assayDataElement(scaledCN, "calls")[valid, ,drop=FALSE]) / rpcs) - cnv)
  l2 = apply(error, 2, function(x) sqrt(sum(x^2)))
  
  pd = Biobase::pData(scaledCN)
  pd$l2 = l2 
  Biobase::pData(scaledCN) = pd
  return(scaledCN)
}

# Compute information theoretic measures for scaled and segmented QDNAseq object
computeInfotheo <- function(scaledCN){
  
  require(infotheo)
  if(is.null(featureData(scaledCN)[["replicationTiming"]])){
    warning("Replication timing not available for binsize - aborting")
    return(scaledCN)
  }
  
  computeInfotheoSingle <- function(i){
    object = scaledCN[, i]
    valid = binsToUseInternal(object)
    rT = featureData(object)[["replicationTiming"]][valid]
    gc = featureData(object)[["gc"]][valid]
    X = Biobase::assayDataElement(object, "copynumber")[valid]
    Y = Biobase::assayDataElement(object, "calls")[valid]
    covariate_info = estimate_gc_correction(object)[valid]
    
    # replication time
    cellcycle.cmi_xy = condinformation(X, Y, S=infotheo::discretize(rT))
    cellcycle.cmi_yrT = condinformation(Y, infotheo::discretize(rT), S=X)
    
    # cell quality
    mi.cmi_yGC = condinformation(Y, infotheo::discretize(gc), S=X)
    mi.yGC = mutinformation(Y, infotheo::discretize(gc))
    
    mi.cmi_xy = condinformation(X, Y, infotheo::discretize(covariate_info))
    mi.xy = mutinformation(X, Y)
    mi.xm = mutinformation(X, infotheo::discretize(covariate_info))
    
    return(c(cellcycle.cmi_xy, cellcycle.cmi_yrT, mi.cmi_yGC, mi.yGC, mi.cmi_xy, mi.xy, mi.xm))
  }
  
  ve = do.call("rbind", lapply(1:ncol(scaledCN), computeInfotheoSingle))
  colnames(ve) = c("cellcycle.cmi_xy", "cellcycle.cmi_yrT", "mi.cmi_yGC", "mi.yGC", "mi.cmi_xy", "mi.xy", "mi.xm")
  
  pd = Biobase::pData(scaledCN)
  ve = as.data.frame(ve)
  rownames(ve) = rownames(pd)
  pd2 = dplyr::bind_cols(pd, ve)
  Biobase::pData(scaledCN) = pd2
  return(scaledCN)
}

# Skewness and kurtosis and their standard errors as implement by SPSS
#
# Reference: pp 451-452 of
# http://support.spss.com/ProductsExt/SPSS/Documentation/Manuals/16.0/SPSS 16.0 Algorithms.pdf
# 
# See also: Suggestion for Using Powerful and Informative Tests of Normality,
# Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.,
# The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321
spssSkew=function(x){
  x = x[!is.na(x)]
  if(length(x) < 10){
    return( NA)
  }
  w=length(x)
  m1=mean(x)
  # m2=sum((x-m1)^2)
  m3=sum((x-m1)^3)
  # m4=sum((x-m1)^4)
  s1=sd(x)
  skew=w*m3/(w-1)/(w-2)/s1^3
  # sdskew=sqrt( 6*w*(w-1) / ((w-2)*(w+1)*(w+3)) )
  # kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)
  # sdkurtosis=sqrt( 4*(w^2-1) * sdskew^2 / ((w-3)*(w+5)) )
  # se only valid for normality assumption
  # mat=matrix(c(skew,kurtosis, sdskew,sdkurtosis), 2,
  #            dimnames=list(c("skew","kurtosis"), c("estimate","se")))
  # return(mat)
  return(skew)
}

# Skewness and kurtosis and their standard errors as implement by SPSS
#
# Reference: pp 451-452 of
# http://support.spss.com/ProductsExt/SPSS/Documentation/Manuals/16.0/SPSS 16.0 Algorithms.pdf
# 
# See also: Suggestion for Using Powerful and Informative Tests of Normality,
# Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.,
# The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321
spssKurtosis=function(x){
  x = x[!is.na(x)]
  if(length(x) < 10){
    return(NA)
  }
  w=length(x)
  m1=mean(x)
  m2=sum((x-m1)^2)
  m4=sum((x-m1)^4)
  s1=sd(x)
  kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)
  return(kurtosis)
}


#' computeHigherOrderMoments
#'
#' \code{computeHigherOrderMoments} compute median skewness across segments for a given QDNAseq object (based on calls)
#' 
#' @param object a (segmented) QDNAseq object
#' 
#' @return higher order moments estimate
computeHigherOrderMoments <- function(object, minLength=20, trimLength=1){
  
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  stopifnot(all(c("calls", "segmented") %in% Biobase::assayDataElementNames(object)))
  require(future.apply, warn.conflicts = FALSE)
  require(matrixStats, warn.conflicts = FALSE)

  processSkewness <- function(obj){
    stopifnot(dim(obj)[[2]] == 1)
    
    valid_bins=binsToUseInternal(obj)
    counts = retrieveAssayData(obj, 1, value="calls", valid=valid_bins)
    gr = retrieveSegmentation(obj, 1)
  
    l_skewness = smoothOperator(counts, gr, function(x){spssSkew(x)}, trimLength = trimLength, minLength=minLength)
    l_kurtosis = smoothOperator(counts, gr, function(x){spssKurtosis(x)}, trimLength = trimLength, minLength=minLength)
    l_weight = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = trimLength, minLength=minLength)
    
    return(list(skewness=matrixStats::weightedMedian(l_skewness, w=l_weight, na.rm=TRUE), kurtosis=matrixStats::weightedMedian(l_kurtosis, w=l_weight, na.rm=TRUE)))
  }
  
  results = future.apply::future_lapply(lapply(colnames(object), function(x) object[, x]), processSkewness)
  skewness = sapply(results, "[[", 1)
  kurtosis = sapply(results, "[[", 2)
  Biobase::pData(object)[["skewness"]] = skewness
  Biobase::pData(object)[["kurtosis"]] = kurtosis
  return(object)
}

#' addObservedVariance
#'
#' \code{addObservedVariance} compute observed variance (as defined by QDNAseq and displayed on the QDNAseq plots as standard deviation). 
#' See sdDiffTrim computation from QDNAseq.
#' 
#' @param object QDNAseq object
#' @param trim 
#' @param scale
#' 
#' @return object with additional slot for "observed.variance"
addObservedVariance <- function(object, trim=0.001, scale=TRUE) {
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  require(Biobase, warn.conflicts = FALSE)
  require(matrixStats, warn.conflicts = FALSE)

  x = Biobase::assayDataElement(object, "copynumber")
  if (scale)
    x <- t(t(x) / base::colMeans(x, na.rm=TRUE))

  vars = apply(x, 2L, function(y) matrixStats::varDiff(y, trim=trim, na.rm=TRUE))
  Biobase::pData(object)[["observed.variance"]] = vars

  return(object)
}

#' readData
#' 
#' Read bam files into a QDNAseq object.
#'
#' \code{readData} returns a QDNAseq copynumber object containing all the cells specified in the bamfiles argument.
#' This is a wrapper function around QDNAseq functionality to read in bam files.
#'
#' @param bamfiles A character vector of full (BAM) file paths.
#' @param binSize Numeric See the QDNAseq package for valid values for binSize. 
#' @param species Character We support Human and Mouse.
#' @param filterChromosomes A vector of chromosome names. By default, the mitochondrial DNA are filtered. 
#' @param extendedBlacklisting filter additional bins based on low read counts (centromer regions)
#'
#' @return A QDNAseq object containing the reads from cells specified through the bamfiles argument
#'
#' @export
readData <- function(bamfiles, binSize, species = "Human", filterChromosomes=c("MT"), genome="hg19",
                     extendedBlacklisting=TRUE,...) {
  
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  require(methods, quietly=TRUE, warn.conflicts = FALSE) # Rscript bugfix for QDNAseq (https://github.com/ccagc/QDNAseq/issues/49)
  require(GenomicRanges, quietly=TRUE, warn.conflicts = FALSE)
  require(Rsamtools, quietly=TRUE, warn.conflicts = FALSE)
  require(dplyr, quietly=TRUE, warn.conflicts = FALSE)
  require(readr, quietly=TRUE, warn.conflicts = FALSE)
  
  transformFun = "none" # non-linear transformations are not supported
  # QDNAseq only supports these bin sizes out of the box
  bins = c(1, 5, 10, 15, 30, 50, 100, 500, 1000)
  stopifnot(binSize %in% bins || (binSize %in% c(200, 2000, 5000) && genome %in% c("hg19", "GRCh37", "GRCm38")))
  

  # required for smaller binsizes
  options(future.globals.maxSize= 4096*1024^2)
  
  if(species == "Human"){
    if(!(genome %in% c("GRCh37", "GRCh38", "hg19", "hg38"))){
      stop("Genome not supported")
    }
    
    if(genome == "hg19"){
      genome = "GRCh37"
    }
    
    if(genome == "hg38"){
      genome = "GRCh38"
    }
    
    if(genome == "GRCh37"){
      require(QDNAseq.hg19, quietly=TRUE, warn.conflicts = FALSE)
      if(binSize %in% c(200, 2000, 5000)){
        binannotations = paste0("hg19.", binSize, "kbp.SR50")
        bins = readRDS(file.path(BASEDIR, "data/customBins", paste0(binannotations, ".RDS")))
      }else{
        binannotations = paste0("hg19.", binSize, "kbp.SR50")
        data(list = c(binannotations), package="QDNAseq.hg19")
        assign("bins", get(binannotations))
      }
    }else{
      if(genome == "GRCh38"){
        require(QDNAseq.hg38, quietly=TRUE, warn.conflicts = FALSE)
        binannotations = paste0("hg38.", binSize, "kbp.SR50")
        data(list = c(binannotations), package="QDNAseq.hg38")
        assign("bins", get(binannotations))
      }
    }
    # bins <- QDNAseq::getBinAnnotations(binSize=binSize)
  } else {
    if(species == "Mouse"){
      require(QDNAseq.mm10, quietly=TRUE, warn.conflicts = FALSE)
      binannotations = paste0("mm10.", binSize, "kbp.SR50")
      data(list = c(binannotations), package="QDNAseq.mm10")
      assign("bins", get(binannotations))
      # bins <- QDNAseq::getBinAnnotations(binSize=binSize, genome = "mm10")
    }else {
      stop("Species is not supported")
    }
  }
  
  isPairedReads = sapply(Rsamtools::BamFileList(bamfiles), Rsamtools::testPairedEndBam)
  if(!(all(isPairedReads) || all(!isPairedReads))){
    stop("Mixing single-end and paired-end bam files is not supported!")
  }
  isPaired = all(isPairedReads)
  

  
  if (all(endsWith(bamfiles, ".bam"))){
    readCounts = QDNAseq::binReadCounts(bins, bamfiles=bamfiles, isDuplicate = FALSE,
                                        isPaired = isPaired, isProperPair=NA,
                                        isFirstMateRead=ifelse(isPaired, TRUE, NA),
                                        isSecondMateRead=ifelse(isPaired, FALSE, NA))
    rawCounts = Biobase::assayDataElement(readCounts, "counts")
    Biobase::pData(readCounts)["bampath"] = dirname(bamfiles)
    Biobase::pData(readCounts)["bamfile"] = basename(bamfiles)
    Biobase::pData(readCounts)["isPaired"] = isPaired
  } else {
    stop("Check input file names, only bam files supported")
  }
  Biobase::pData(readCounts)["binSize"] = binSize
  Biobase::pData(readCounts)["species"] = species
  
  # Compute coverage (as an alternative to reads per cell metric)
  compute_read_size <- function(bamFile){
    yieldSize(bamFile) <- 10000
    # open(bamFile)
    reads = scanBam(bamFile)[[1]]$seq
    
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    return(getmode(reads@ranges@width))
  }
  read_size = sapply(Rsamtools::BamFileList(bamfiles), compute_read_size)
  Biobase::pData(readCounts)["read_size"] = as.numeric(read_size)
  
  if(species == "Human"){
    Biobase::pData(readCounts)["genome"] = genome
    coverage = (readCounts@phenoData@data$used.reads * (as.numeric(isPaired) + 1) * read_size) / 3.2e9
    Biobase::pData(readCounts)["coverage"] = coverage
  }
  
  if(species == "Human" & extendedBlacklisting){
    if(binSize < 10) stop("blacklisting doesn't support binsizes smaller than 10kb")
    # this is because our blacklisted regions have a resolution of 10kb
    if(binSize %in% c(200, 2000, 5000) && genome == "GRCh37"){
      # custom bins
      extendedFilter = Biobase::fData(readCounts)[["use"]]
    }
    readCounts = readFilter(readCounts, genome=genome)
    extendedFilter = Biobase::fData(readCounts)[["use"]]
  }else{
    extendedFilter = Biobase::fData(readCounts)[["use"]] | base::startsWith(rownames(readCounts), "X:") | base::startsWith(rownames(readCounts), "Y:")
  }
  
  # create QNDAseq object ####
  readCounts <- suppressWarnings(QDNAseq::estimateCorrection(readCounts, verbose=FALSE))
  readCountsFiltered <- QDNAseq::applyFilters(readCounts, residual=FALSE, blacklist=TRUE, mappability = NA, chromosomes = filterChromosomes, verbose=FALSE)
  Biobase::fData(readCountsFiltered)[["use"]] = extendedFilter & !is.na(Biobase::fData(readCounts)[["gc"]])
  fdat = Biobase::fData(readCountsFiltered)
  if (binSize %in% c(30, 50, 100, 500, 1000) & species == "Human"){
    repTime = readRDS(file.path(BASEDIR, "data/replicationTiming/replicationTiming_per_binSize.RDS"))[[as.character(binSize)]]
    fdat$replicationTiming = repTime$replicationTime
  } else if (binSize %in% c(30, 50, 100, 500, 1000) & species == "Mouse"){
    repTime = readRDS(file.path(BASEDIR, "data/replicationTiming/replicationTiming_per_binSize_mouse.RDS"))[[as.character(binSize)]]
    fdat$replicationTiming = repTime$replicationTime
  }

  copyNumbers <- new("QDNAseqCopyNumbers", bins=fdat,
                     copynumber=Biobase::assayDataElement(readCountsFiltered, "counts"),
                     phenodata=Biobase::pData(readCountsFiltered))
  
  cn = Biobase::assayDataElement(copyNumbers, "copynumber")
  keep_idx = binsToUseInternal(copyNumbers)
  cn[!keep_idx] = NA

  copyNumbers = Biobase::assayDataElementReplace(copyNumbers, "copynumber", cn)

  # add value for observed variance (sqrt of value in plots)
  copyNumbers = addObservedVariance(copyNumbers)

  # remove dummy variables - we do not apply GC correction at this stage
  Biobase::pData(copyNumbers)[["loess.span"]] = NULL
  Biobase::pData(copyNumbers)[["loess.family"]] = NULL

  copyNumbers = Biobase::assayDataElementReplace(copyNumbers, "calls", assayDataElement(copyNumbers, "copynumber"))
  copyNumbers = Biobase::assayDataElementReplace(copyNumbers, "probdloss", rawCounts)
  
  if(species == "Human"){
    Biobase::experimentData(readCounts)@other = list(species="Human", genome=genome)
  }else{
    Biobase::experimentData(readCounts)@other = list(species=species)
  }
  
  return(copyNumbers)
}

#' \code{readFilter} apply an updated bin filter to exclude low quality and unreliable bins from the estimation
#'
#' @param readCounts QDNAseq object
#' 
#' @return QDNAseq object with updated bin usability
readFilter <- function(readCounts, genome="GRCh37"){

  stopifnot(genome %in% c("GRCh37", "GRCh38"))
  if(genome == "GRCh37"){
    gr = createGR(readCounts[,1,drop=FALSE])
    excluded_regions_bed = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/final_exclude_regions_hg19.bed"),
                                           col_names = c("chromosome", "start", "end", "A" ,"B", "C"), col_types="ciiccc")
    excluded_regions = GRanges(seqnames=excluded_regions_bed$chromosome,
                               ranges = IRanges(start=excluded_regions_bed$start+1, end=excluded_regions_bed$end),
                               seqinfo = seqinfo(gr))
  }
  if(genome == "GRCh38"){
    warning("blacklist has been optimized for GRCh37 only")
    stop("Not supported")
  }

  hits <- findOverlaps(excluded_regions, gr)

  # Compute percent overlap and filter the hits:
  overlaps <- pintersect(excluded_regions[queryHits(hits)], gr[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(gr[subjectHits(hits)])
  hits <- hits[percentOverlap > 0.25]

  pass_qc = Biobase::fData(readCounts)[["use"]] | base::startsWith(rownames(readCounts), "X:") | base::startsWith(rownames(readCounts), "Y:")
  pass_qc[subjectHits(hits)] = FALSE
   
  Biobase::fData(readCounts)[["use"]] = pass_qc
  return(readCounts)
}


#' readFlagstat
#' 
#' Read flagstat file to get additional information on read quality
#'
#' @param bamfiles A character vector of (BAM) file names.
#'
#' @return a data frame with the flagstat information.
#'
#' @export
readFlagstat <- function(bamfiles){
  
  read_flagstat_sub <- function(fh){
    # Read flagstat results
    require(stringr, quietly=TRUE, warn.conflicts = FALSE)
    f <- readLines(fh)
    a <- data.frame( as.character(sub("\\.flagstat$", "", base::basename(fh))),
                     as.numeric(stringr::word(f[1]) )
                     , as.numeric(stringr::word(f[2]) )
                     , as.numeric(stringr::word(f[4])  )
                     , as.numeric(stringr::word(f[5])  )
                     , as.numeric(stringr::word(f[6])  )
                     , as.numeric(stringr::word(f[9])  )
                     , as.numeric(stringr::word(f[10]) )
                     , as.numeric(stringr::word(f[11]) )
                     , as.numeric(stringr::word(f[12]) )
                     , as.numeric(stringr::word(f[13]) )
    )
    names(a) <- NULL
    colnames(a) <- c("name", "qcTotal", "qcSecondary", "qcDuplicates", "qcMapped", "qcPairedSequencing", "qcProperlyPaired", "qcPairedWithItself",
                     "qcSingletons", "qcMateMapped", "qcMateMappedmapQ5")
    return(a)
  }
  
  fh = bamfiles
  fh_list = c()
  for (i in 1:length(fh)){
    if(endsWith(fh[i], ".bam")){
      fh_list <- c(fh_list, fh[i])
    }else{
      if(dir.exists(fh[i])){
        files <- list.files(path = fh[i], pattern = "\\.bam$")
        files = paste0(fh[i], files)
        fh_list <- c(fh_list, files)
      }
    }
  }
  
  flagstat_description = list()
  counter = 1
  flaglist = sub("\\.bam$", ".flagstat", fh_list)
  for(i in 1:length(flaglist)){
    if(file.exists(flaglist[i]) && !dir.exists(flaglist[i])){
      flagstat_description[[counter]] = read_flagstat_sub(flaglist[i])
      counter = counter + 1
    }
  }
  
  result <- do.call("rbind", flagstat_description)
  return(result)
}




readPosition <- function(object, filePaths, selectRegion=c(as.character(c(1:22)), "X"), debug=FALSE){
  
  require(Biobase)
  require(QDNAseq)
  require(MASS)
  require(GenomicRanges)
  require(IRanges)
  require(GenomicAlignments)
  require(S4Vectors)
  stopifnot("name" %in% colnames(Biobase::pData(object)))
  
  readPositionSingle <- function(obj, fP){
    valid = binsToUseInternal(obj)
    high_quality_regions = createGR(obj[,1])
    high_quality_regions = high_quality_regions[valid]
    chromosome = decode(seqnames(high_quality_regions))
    copynumber = Biobase::assayDataElement(obj, "copynumber")[valid]
    name = Biobase::pData(obj)[["name"]]

    read_position_table = readr::read_tsv(paste0(sub('\\.bam$', '', fP), ".position.tsv"), col_types = "cdd")
    
    read_table = GRanges(read_position_table$chromosome,
                     IRanges(read_position_table$start, read_position_table$start+abs(read_position_table$length-1)),
                     length=abs(read_position_table$length),
                     dist=c(diff(read_position_table$start), Inf))
    
    read_table = GenomicRanges::sort(read_table, ignore.strand=TRUE)
    # identify reads in high quality regions
    hits = findOverlaps(read_table, high_quality_regions, type="within")
    hq = read_table[queryHits(hits)]
    
    # peaks 
    chromosome = decode(seqnames(high_quality_regions))
    copynumber = Biobase::assayDataElement(obj, "copynumber")[valid]
    start = high_quality_regions@ranges@start 
    end = high_quality_regions@ranges@start + high_quality_regions@ranges@width -1
    
    cov = GenomicRanges::coverage(hq)
    peaks = IRanges::slice(cov, lower=10)
    peak_maxima = IRanges::viewMaxs(peaks)
    peak_means = IRanges::viewMeans(peaks)
    gr = reduce(GRanges(peaks))
    gr$maxima = unlist(peak_maxima)
    gr$means = unlist(peak_means)
    #gr$name = name
    #gr$status = "extrema"
    hits = findOverlaps(gr, high_quality_regions, type="within")
    gr$overlap_start = gr@ranges@start
    gr$overlap_end =  gr@ranges@start + gr@ranges@width - 1
    gr$bin_start = start[subjectHits(hits)]
    gr$bin_end =  end[subjectHits(hits)]
    gr$copy = copynumber[subjectHits(hits)]
    gr$chrom = chromosome[subjectHits(hits)]
      
    # overlap regions (average) 
    peaks_overlap = cov
    for(ls in names(peaks_overlap)){
      y = peaks_overlap[[ls]] 
      runValue(y)[runValue(y) == 1] = 0;
      peaks_overlap[[ls]] = y
    }
    no_overlap = cov
    for(ls in names(no_overlap)){
      y = no_overlap[[ls]] 
      runValue(y)[runValue(y) != 0] = -1;
      runValue(y)[runValue(y) == 0] = 1;
      runValue(y)[runValue(y) == -1] = 0;
      no_overlap[[ls]] = y
    }
    
    common_levels = base::intersect(names(peaks_overlap), unique(seqnames(high_quality_regions)))
    peaks_overlap = peaks_overlap[common_levels]
    cov = cov[common_levels]
    no_overlap = no_overlap[common_levels]
    high_quality_regions = keepSeqlevels(high_quality_regions, common_levels, pruning.mode="coarse")
    
    gr_overlap_2 = binnedAverage(high_quality_regions, peaks_overlap, "overlap")
    gr_overlap_all = binnedAverage(high_quality_regions, cov, "overlap")
    gr_no_overlap = binnedAverage(high_quality_regions, no_overlap, "no_overlap")
    
    dens_average = dplyr::tibble(chrom = chromosome, copy=copynumber,
                         overlap=gr_overlap_2$overlap,
                         overlap_all=gr_overlap_all$overlap,
                         no_overlap=gr_no_overlap$no_overlap, status="average")
    
    dens_extrema = as.data.frame(gr)
    if(nrow(dens_extrema) > 0){
      dens_extrema$status = "extrema"
    }
    dens = dplyr::bind_rows(dens_average, dens_extrema)
  
    dens$name = Biobase::pData(obj)[["name"]]
    dens$rpc = Biobase::pData(obj)[["rpc"]]
    
    return(dens)
  } 
  
  densTable = lapply(1:ncol(object), function(index) readPositionSingle(object[, index], filePaths[index]))
  densTable = do.call("rbind", densTable) %>% dplyr::ungroup()
  
  meta = base::data.frame(labelDescription = rep(NA, times=dim(densTable)[[2]]))
  rownames(meta) = colnames(densTable)
  outProtData <- new("AnnotatedDataFrame",
                      data=densTable, 
                      varMetadata = meta)
  Biobase::protocolData(object) = outProtData

  return(object)
}


#' initial_gc_correction
#'
#' \code{initial_gc_correction} corrects segmentedCounts object for gc and map variation
#'
#' @param segmentedCounts QDNAseq object with segmented slot
#' @param debug Boolean flag to add debug information to output
#' 
#'
#' @return A QDNAseq object with raw reads from call slot corrected and added to copynumber slot
#'
#' @export
initial_gc_correction <- function(segmentedCounts, debug=FALSE){
  suppressPackageStartupMessages(require(mgcv, quietly=TRUE, warn.conflicts = FALSE))
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  
  iterate_initial_gc_correction <- function(object){
    
    stopifnot(dim(object)[[2]] == 1)
    n_bins = dim(object)[[1]]
    valid = binsToUseInternal(object)

    reads = Biobase::assayDataElement(object[valid, 1], "calls")
    segs = Biobase::assayDataElement(object[valid, 1], "segmented")
    gc = Biobase::fData(object[valid, 1])[["gc"]]
    map = Biobase::fData(object[valid, 1])[["mappability"]]
    stopifnot(length(gc) == length(reads))

    # log transform seg level for better stability
    dat = data.frame(reads=as.vector(reads), segs=log(as.vector(segs) + 1),
                     reads_normalized=as.vector(reads)/as.vector(segs),
                     gc=gc, map=map)
    
    ## Build models
    if(all(map == 100)){
      formula.gam = formula(reads ~ segs + s(gc))
      newdata=data.frame(segs = dat$segs, gc=dat$gc)
    }else{
      formula.gam = formula(reads ~ segs + s(gc, map))
      newdata=data.frame(segs = dat$segs, gc=dat$gc, map=dat$map)
    }
    dat_model = dat %>% dplyr::filter(segs > 0)
    model <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    terms = predict(model, newdata=newdata, type="terms")
    resp = predict(model, newdata=newdata, type="response")
    # quantile(resp)
    # max_reads = sort(resp, decreasing = TRUE)[1:10]
    # dat[resp %in% max_reads,]
    # sum(dat$reads) - sum(round(resp))
    
    
    # baseline without gc vs gc-map model
    # formula.gam.baseline = formula(reads ~ segs + s(gc))
    formula.gam.baseline = formula(reads ~ segs)
    model.baseline <- mgcv::gam(formula.gam.baseline, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    
    # formula.gam.both = formula(reads ~ segs + s(gc) + s(map))
    # formula.gam.both = formula(reads ~ segs + s(gc, map))
    model.baseline.both <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    
    av = anova(model.baseline, model.baseline.both, test="Chisq")
    ## end baseline anova
    
    ## Remove covariate effect
    
    # response = exp(terms[,1] + terms[,2] + coef(model)[1])
    # stopifnot(all(abs(response - resp) < 1e-6))
    # mean(terms[,2])
    response_gc_and_map = exp(terms[,2])
    residuum = dat$reads * (1/response_gc_and_map)
    
    dat$pred_seg = terms[,1]
    dat$pred_cov = terms[,2]
    dat$residuum = residuum

    stopifnot(length(residuum) == length(reads))
    corrected = round(residuum, digits=0)
    dat$corrected_reads = corrected
    negative_elements = sum(corrected < 0, na.rm=TRUE)
    corrected[corrected < 0] = 0
    stopifnot(all(!is.na(corrected)))
    
    corrected_reads = Biobase::assayDataElement(object, "calls")
    stopifnot(length(corrected_reads[valid]) == length(corrected))
    corrected_reads[valid,1] = corrected
    
    corrected_reads_residuum = Biobase::assayDataElement(object, "calls")
    corrected_reads_residuum[valid,1] = residuum
    
    # sum(dat$reads) - sum(dat$corrected)
    # table(dat$corrected != dat$reads)
    # pl1 = ggplot(data=dat) + geom_point(aes(x=gc, y=reads_normalized)) +
    #   # geom_point(aes(x=gc, y=reads_normalized, color=as.factor(state)), size=0.1) +
    #   geom_point(aes(x=gc, y=reads_normalized), size=0.1) +
    #   geom_smooth(aes(x=gc, y=reads_normalized), formula="y ~ x", method="loess", span=0.7, color="red", se=FALSE) +
    #   geom_vline(xintercept = c(37, 43), color="red") + theme_cowplot()
    # pl1
    
    # pl2 = ggplot(data=dat1) + geom_point(aes(x=gc, y=reads)) +
    #   geom_point(aes(x=gc, y=reads, color=as.factor(state)), size=0.1) +
    #   geom_smooth(aes(x=gc, y=reads), formula="y ~ x", method="loess", span=0.5, color="red", se=FALSE) +
    #   geom_vline(xintercept = c(37, 43), color="red") + theme_cowplot()
    # pl2

    # ggplot(data=dat1) + geom_histogram(aes(x=ratio1)) + theme_cowplot() + coord_cartesian(xlim=c(-0.3, 0.3))
    # dat1$corrected_reads1 = dat1$reads - round(dat1$ratio1 * dat1$reads, digits=0)
    # dat1$corrected_reads2 = dat1$reads - round(dat1$ratio2 * dat1$reads, digits=0)
    # ggplot(data = dat1 %>%
    #     dplyr::select(reads, corrected_reads2) %>% tidyr::gather(condition, r, c(corrected_reads2, reads))) +
    #     geom_histogram(aes(x=r, fill=condition), binwidth = 1, position="identity", alpha=0.75) +
    #     theme_cowplot()
    # 
    # pl3 = ggplot(data=dat1) + geom_point(aes(x=gc, y=corrected_reads2)) +
    #   geom_smooth(aes(x=gc, y=corrected_reads2), formula="y ~ x", method="loess", span=0.75, color="red", se=FALSE) +
    #   geom_smooth(aes(x=gc, y=reads), formula="y ~ x", method="loess", span=0.75, color="blue", se=FALSE) +
    #   theme_cowplot()
    # pl3
    
    # # #Visualize mgcv model
    # options(rgl.useNULL=TRUE)
    # library(mgcViz) # conda dependencies: rgl, qgam, make surn KernSmooth is correct version (conda-forge)
    # 
    # model.viz <- getViz(model)
    # o <- plot( sm(model.viz, 1) )
    # o + theme_pubr()
    
    ## Read correction
    # dat1$corrected_reads = dat1$reads - round(ratio * dat1$reads, digits=0)
    
    # # smoothing per cn state
    # pl3 = ggplot(data=dat1) + geom_point(aes(x=gc, y=corrected_reads, color=as.factor(state))) + 
    #   geom_smooth(aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="red", se=FALSE) +
    #   geom_smooth(data=dat1 %>% dplyr::filter(state==1), aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="blue", se=FALSE) +
    #   geom_smooth(data=dat1 %>% dplyr::filter(state==2), aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="blue", se=FALSE) +
    #   geom_smooth(data=dat1 %>% dplyr::filter(state==3), aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="blue", se=FALSE) +
    #   geom_vline(xintercept = c(37, 43), color="red") + theme_cowplot() #+ coord_cartesian(ylim=c(200, 300)
    # pl3

    # ggplot(data = dat %>%
    #     dplyr::select(reads, corrected_reads) %>% tidyr::gather(condition, r, c(corrected_reads, reads))) +
    #     geom_histogram(aes(x=r, fill=condition), binwidth = 1, position="identity", alpha=0.75) +
    #     coord_cartesian(xlim=c(200, 500)) +
    #     theme_cowplot()
    
    object = Biobase::assayDataElementReplace(object, "copynumber", matrix(data=corrected_reads, ncol=1))
    desc = Biobase::pData(object)
    
    desc$gc_correction = TRUE
    desc$gc_correction.negative_elements = negative_elements
    desc$gc_correction.av.deviance = av[["Deviance"]][2]
    desc$gc_correction.av.p_value = av[["Pr(>Chi)"]][[2]]
    desc$gc_correction.dev.expl = summary(model)[["dev.expl"]]
    desc$gc_correction.alpha = 1/model$family$getTheta(TRUE)
    desc$gc_correction.delta_n_reads = sum(Biobase::assayDataElement(object, "calls"), na.rm = TRUE) - sum(Biobase::assayDataElement(object, "copynumber"), na.rm=TRUE)
    Biobase::pData(object) = desc
    return(list("obj"=object, list("name"=colnames(object), "dat"=dat, "model"=model)))
  }
  
  results = future.apply::future_lapply(lapply(colnames(segmentedCounts), function(x) segmentedCounts[, x]), iterate_initial_gc_correction)
  new_object = combineQDNASets(sapply(results, "[[", 1))
  debugInformation = sapply(results, "[[", 2)
  
  if(debug){
    return(list("new_object"=new_object, "debugInformation"=debugInformation) )
  }else{
    return(new_object)  
  }
}


#' estimate_gc_correction
#'
#' \code{estimate_gc_correction} estimate gc correction per cell
#'
#' @param object a QDNAseq object with segmented slot
#'
#' @return a n_bins x n_cells matrix with offset for gc/map correction
#'
#' @export
estimate_gc_correction <- function(segmentedCounts){
  
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)  
  n_cells = dim(segmentedCounts)[[2]]
  valid = binsToUseInternal(segmentedCounts)
  
  iterate_processing <- function(x){
    suppressPackageStartupMessages(require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE))
    obj = segmentedCounts[, x]
    
    stopifnot(dim(obj)[[2]] == 1)
    n_bins = dim(obj)[[1]]
    valid = binsToUseInternal(obj) & !is.na(Biobase::fData(obj)[["gc"]])
    
    reads = Biobase::assayDataElement(obj, "calls")[valid, 1, drop=TRUE]
    segs = Biobase::assayDataElement(obj, "segmented")[valid, 1, drop=TRUE]
    gc = Biobase::fData(obj[valid, 1])[["gc"]]
    map = Biobase::fData(obj[valid, 1])[["mappability"]]
    stopifnot(length(gc) == length(reads))
    
    # log transform seg level for better stability
    dat = data.frame(reads=as.vector(reads), segs=log(as.vector(segs) + 1),
                     reads_normalized=as.vector(reads)/as.vector(segs),
                     gc=gc, map=map)
    
    ## Build models
    if(all(map == 100)){
      formula.gam = formula(reads ~ segs + s(gc))
      newdata=data.frame(segs = dat$segs, gc=dat$gc)
    }else{
      formula.gam = formula(reads ~ segs + s(gc, map))  
      newdata=data.frame(segs = dat$segs, gc=dat$gc, map=dat$map)
    }
    dat_model = dat %>% dplyr::filter(segs > 0)
    model <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    terms = predict(model, newdata=newdata, type="terms")
    
    gc_terms = terms[,2]
    cor = rep(0, n_bins)
    cor[valid] = exp(gc_terms)
    return(cor)
  }
  
  xs <- seq(1, dim(segmentedCounts)[[2]])
  results = future.apply::future_lapply(as.list(xs), iterate_processing)
  mat = t(do.call("rbind", results))
  return(mat)
}


#' estimate_overdispersion
#'
#' \code{estimate_overdispersion} estimate alpha value per cell
#'
#' @param object a QDNAseq object with copynumber slot
#'
#' @return vector of alpha values
#'
#' @export
estimate_overdispersion <- function(object, initial_estimate_alpha=0.1, robust=NULL){
  # robust can be used to exclude copy number states, e.g. c(0, 1, 8)
  valid = binsToUseInternal(object) 
  
  for(cellindex in 1:ncol(object)){
    
    alpha_cell <- tryCatch({
    cn_states = Biobase::assayDataElement(object, "copynumber")[valid, cellindex]
    # robust can be max copy number state (estimate for high values uncertain)
    if(!is.null(robust)){
      cn_levels = setdiff(unique(cn_states), c(0, robust))
    }else{
      cn_levels = setdiff(unique(cn_states), 0)
    }
    
    alpha = c()
    alpha_estimate = c()
    alpha_estimate_size = c()
    rpc = Biobase::phenoData(object[, cellindex])[["rpc"]]
    
    for(cn_l in cn_levels){
      
      tmp <- MASS::fitdistr(Biobase::assayDataElement(object, "calls")[valid, cellindex][cn_states == cn_l],
                            dnbinom, list(mu=rpc * cn_l, size=1/initial_estimate_alpha),
                            method = "L-BFGS-B",lower = c(0.01,1e-3),upper = c(10000,1e9))
      alpha_estimate_cnl = 1.0/tmp$estimate[["size"]]
      alpha_estimate = c(alpha_estimate, alpha_estimate_cnl)
      alpha_estimate_size = c(alpha_estimate_size, sum(cn_levels == cn_l))
    }
    alpha_cell = stats::weighted.mean(alpha_estimate, alpha_estimate_size)
    return(alpha_cell)
   },
    error=function(cond) {
      # Choose a return value in case of error
      return(NA)
    })    
    
    alpha = c(alpha, alpha_cell)
  }
  return(alpha)
}


#' combineQDNASets
#'
#' \code{combineQDNASets} combines multiple QDNAseq objects to a single object, useful for parallelisation across cells.
#'
#' @param ... list of QDNAseq objects
#' @param dropProtocol Boolean size of protocolData can be enormous if all metadata is included, can be droped for convenience (recommend extracting info separately if required)
#'
#' @return A QDNAseq object containing all the input objects
#'
#' @export
combineQDNASets <- function(..., reduceMetadata=FALSE, dropProtocol=FALSE){
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  require(dplyr, quietly=TRUE, warn.conflicts = FALSE)


  x = ...[[1]]

  if(length(...) == 1){
    # idempotence
    return(x)
  }

  # keep information from x, append assayData and other data from other elements
  if("segmented" %in% names(Biobase::assayData(x))){
    outSegmented <- do.call(cbind, lapply(..., Biobase::assayDataElement, "segmented"))
  }
  outCopynumber <- do.call(cbind, lapply(..., Biobase::assayDataElement, "copynumber"))
  if("calls" %in% names(Biobase::assayData(x))){
    outCalls <- do.call(cbind, lapply(..., Biobase::assayDataElement, "calls"))
  }
  if("probdloss" %in% names(Biobase::assayData(x))){
    outProbdloss <- do.call(cbind, lapply(..., Biobase::assayDataElement, "probdloss"))
  }

  ## check pheno data is compatible and indicate source of error
  test = lapply(..., pData)
  n_names = unlist(lapply(..., BiocGenerics::colnames))
  n_items = unlist(lapply(test, length))
  n_max = max(n_items)
  if(any(n_items < n_max)){
    print(n_names[n_items == n_max])
    print("---")
    print(n_names[n_items < n_max])
    warning(paste0("INCONSISTENT METADATA -> CHECK ", sum(n_items < n_max), " OFFENDING SAMPLES"))
    if(reduceMetadata){
      commonNames = colnames(test[[1]])
      for(item in 1:length(test)){
        commonNames = intersect(commonNames, colnames(test[[item]]))
        if(length(setdiff(colnames(test[[item]]), commonNames)) > 0){
          print(paste0("Items in question: ", setdiff(colnames(test[[item]]), commonNames)))
        }
      }
    }else{
      stop("stop to avoid inconsistent metadata")
    }
  }else{
    commonNames = colnames(test[[1]])
  }

  meta = Biobase::varMetadata(x)
  meta = meta[rownames(meta) %in% commonNames,,drop=FALSE]
  outPhenoData <- new("AnnotatedDataFrame",
                      data=base::do.call(base::rbind, lapply(..., function(x){y = pData(x); return(y[, commonNames, drop=FALSE])})),
                      varMetadata = meta)
  
  outFeatureData <- new("AnnotatedDataFrame",
                        data=fData(x),
                        varMetadata = Biobase::fvarMetadata(x))

  outProtocolData = protocolData(x)
  if(dim(outProtocolData)[[2]] > 0){
    if(!dropProtocol){
      protDat = data.frame(dplyr::bind_rows(lapply(..., function(x) Biobase::protocolData(x)@data)))
      rownames(protDat) = as.character(seq(1, dim(protDat)[[1]]))
    }else{
      protDat = data.frame(dplyr::bind_rows(lapply(..., function(x){return(Biobase::protocolData(x)@data %>% dplyr::filter(status == "extrema"))})))
      rownames(protDat) = as.character(seq(1, dim(protDat)[[1]]))
    }
    outProtocolData@data = protDat
  }

  result = new("QDNAseqCopyNumbers",
               bins=outFeatureData,
               copynumber=outCopynumber,
               phenodata=outPhenoData)
  if("segmented" %in% names(Biobase::assayData(x))){
    Biobase::assayDataElement(result, "segmented") <- outSegmented
  }
  if("calls" %in% names(Biobase::assayData(x))){
    Biobase::assayDataElement(result, "calls") <- outCalls
  }
  if("probdloss" %in% names(Biobase::assayData(x))){
    Biobase::assayDataElement(result, "probdloss") <- outProbdloss
  }
  Biobase::protocolData(result) <- outProtocolData
  
  return(result)
}


#' predict_replicating
#'
#' \code{predict_replicating} update metadata with (predicted) value of cycling status
#   predict replicating cells in each dataset
#  (requires dft to have UID, SLX, and optionally TECHNOLOGY covariates, depending on batch variable)
#   cycling activity measure is cellcycle.cmi_yrT
#'
#' @param dft phenoData object of QDNAseq object after running through scAbsolute pipeline
#' @param cutoff_value Numeric number of standard deviations to use for outlier exclusion (cellcycle.cmi_yrT, fitted with normal dist)
#' @param iqr_value Numeric number of iqr ranges to exclude outlier (cellcycle.kendall.repTime.weighted.median.cor.corrected)
#' @param batch Character batch to use to group data on, before applying normal approximation
#' @param hard_threshold Boolean remove data points higher than hard_threshold cycling activity for data sets with
#' imbalanced cases (i.e. where cycling mode is higher than hard_threshold)
#'
#' @return updated data frame with cycling cell information
#'
#' @export
predict_replicating <- function(dft, cutoff_value=1.50, iqr_value=1.50, batch="sample", hard_threshold=NULL){
  
  stopifnot(batch %in% c("all", "sample", "technology"))
  require(stats)
  if(batch == "all"){
    stopifnot("UID" %in% colnames(dft))
    stopifnot("SLX" %in% colnames(dft))
    stopifnot("TECHNOLOGY" %in% colnames(dft))
    dftin = dft %>% dplyr::group_by(UID, SLX, TECHNOLOGY)
  }
  if(batch == "sample"){
    stopifnot("UID" %in% colnames(dft))
    stopifnot("SLX" %in% colnames(dft))
    dftin = dft %>% dplyr::group_by(UID, SLX)
  }
  if(batch == "technology"){
    stopifnot("TECHNOLOGY" %in% colnames(dft))
    dftin = dft %>% dplyr::group_by(TECHNOLOGY)
  }
  stopifnot(("cycling_activity_" %in% colnames(dft)) || ("cellcycle.kendall.repTime.weighted.median.cor.corrected" %in% colnames(dft)))
  
  # analysis within each dataset
  dft1 = dftin %>% 
    dplyr::group_modify(function(x, y){
      
      # DEBUG
      #x = df_ploidy %>% dplyr::filter(UID %in% c("UID-NNA-TN3"))
      #x = df_ploidy %>% dplyr::filter(UID %in% c("UID-DLP-SA1090"))
      if(!("cycling_activity" %in% colnames(dft))){
        x$cycling_activity = x$cellcycle.kendall.repTime.weighted.median.cor.corrected
      }
      
      # compute mode of distribution
      x$cycling_median = median(round(x$cycling_activity, digits = 3))
     
      # gini boxplot outlier - for cycling activity
      Q_3 = quantile(x$cycling_activity, 0.75)
      Q_1 = quantile(x$cycling_activity, 0.25)
      IQR = Q_3 - Q_1
      upper_whisker = min(max(x$cycling_activity), Q_3 + iqr_value * IQR)
      
      # sd of distribution (assume this is half-normal) below mode:
      lower_whisker = max(min(x$cycling_activity), Q_1 - (IQR*1.5)) #lower whisker for boxplot
      sd_est = sqrt((1 / (1 - (2/pi))) * (sd(x$cycling_activity[(x$cycling_activity < x$cycling_median) & (x$cycling_activity > lower_whisker)], na.rm=TRUE)^2))
      
      # gini boxplot outlier - for kendall test (only extrema)
      if("cellcycle.cmi_yrT" %in% colnames(x) && !is.null(iqr_value)){
        Q_3 = quantile(x$cellcycle.cmi_yrT, 0.75)
        Q_1 = quantile(x$cellcycle.cmi_yrT, 0.25)
        IQR = Q_3 - Q_1
        upper_whisker = min(max(x$cellcycle.cmi_yrT), Q_3 + (iqr_value * IQR))
        x$cycling_cutoff_iqr = upper_whisker
      }
      
      x$cycling_sd_est = sd_est
      x$cycling_cutoff_sd = x$cycling_median + (cutoff_value * x$cycling_sd_est)
      if("cellcycle.cmi_yrT" %in% colnames(x) && !is.null(iqr_value)){
        x$replicating = (x$cycling_activity > x$cycling_cutoff_sd) | (x$cellcycle.cmi_yrT > x$cycling_cutoff_iqr)
      }else{
        x$replicating = x$cycling_activity > x$cycling_cutoff_sd
      }

      if(!is.null(hard_threshold)){
        x$replicating = x$replicating | (x$cycling_activity > hard_threshold)
      }
      
      
      
      #ggplot(data=x) + geom_quasirandom(aes(x="NA", y=cycling_activity, color=replicating))
      
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}

# flag ploidy outlier in phenoData based on SLX, UID grouping - helper function
qc_ploidy <- function(dft, range=0.5){
  require(stats)
  stopifnot("UID" %in% colnames(dft))
  stopifnot("SLX" %in% colnames(dft))
  stopifnot("ploidy" %in% colnames(dft))
  
  dft1 = dft %>% dplyr::group_by(UID, SLX) %>% 
    dplyr::group_modify(function(x, y){
      
      #x = df2 %>% dplyr::filter(SLX == "SLX-A73044A")
      # NOTE: median is important here to get unbiased estimate of mode
      mp = median(x$ploidy, na.rm=TRUE)
      x$ploidy.outlier = (x$ploidy < (mp - range)) | (x$ploidy > (mp + range))
    
      return(x)
    }) %>% dplyr::ungroup()
  
  return(dft1)
}


#' getSegTable
#'
#' \code{getSegTable} return a bed style representation of a QDNAseq object's segmentation
#'
#' @param object QDNAseq object with segmentation slot
#' @param trimLength Numeric number of bins to trim from start and end of segments
#'
#' @return a GRanges List object with n_samples elements, representing the individual samples segmentation profiles
#'
#' @export
getSegTable <- function(object, trimLength = 1, slot="segmented"){
  require(GenomicRanges, quietly=TRUE, warn.conflicts = FALSE)
  require(IRanges, quietly=TRUE, warn.conflicts = FALSE)
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  
  export_genomic_ranges_call <- function(sample){
    
    gr = retrieveSegmentation(object, sample, slot=slot)
    counts = retrieveAssayData(object, sample, "copynumber")
    
    n_segments = length(gr@ranges)

    if(!is.null(Biobase::assayDataElement(object, "calls"))){
      gr$reads = smoothOperator(retrieveAssayData(object, sample, "calls"), gr, function(x){stats::median(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
      gr$var = smoothOperator(retrieveAssayData(object, sample, "calls"), gr, function(x){stats::var(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
      gr$mean = smoothOperator(retrieveAssayData(object, sample, "calls"), gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
      gr$sample_mean = base::mean(assayDataElement(object[, sample], "calls"), na.rm=TRUE)
      gr$sample_var = stats::var(assayDataElement(object[, sample], "calls"), na.rm=TRUE)
    }
    gr$copynumber = smoothOperator(retrieveAssayData(object, sample, "copynumber"), gr, function(x){stats::median(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
    gr$support = smoothOperator(retrieveAssayData(object, sample, "copynumber"), gr, function(x){length(na.omit(x))}, trimLength = 0, minLength=1)
    gr$length = gr@ranges@width
    gr$sample = sample
    
    ## convert to genomic coordinates (bp coordinates)
    gra = createGR(object[,sample])
    ir_new = do.call("c", lapply(seqlevels(gr), function(x){
      grat = gra[seqnames(gra) == x]
      grt = gr[seqnames(gr) == x]
      
      if(length(grt) == 0){
        return(IRanges())
      }
      
      stopifnot(all(length(grat@ranges) >= grt@ranges@start))
      binsize = grat@ranges@width[[1]]
      irStart = c();irEnd = c(); irWidth = c();
      for(idx in 1:length(grt)){
        start_grt = grt@ranges@start[idx]
        end_grt = grt@ranges@start[idx] + grt@ranges@width[idx] - 1
        start = grat@ranges@start[start_grt]
        end = grat@ranges@start[start_grt] + sum(grat@ranges@width[start_grt:end_grt]) - 1
        width = end - start + 1
        irStart = c(irStart, start); irEnd = c(irEnd, end);
        irWidth = c(irWidth, width)
      }
      
      irnew = IRanges::IRanges(start = irStart,
                               end = irEnd, 
                               width = irWidth)
      
      return(irnew)
    }))
    
    gr@ranges = ir_new
    # end coordinate transformation
    
    return(gr)
  }
  
  result = lapply(colnames(object), export_genomic_ranges_call)
  return(GRangesList(result))
}

#' writeSegTable
#'
#' \code{writeSegTable} combines multiple QDNAseq objects to a single object, useful for parallelisation across cells.
#'
#' @param object QDNAseq object
#' @param filename filename to write bed file to
#' @param trimLength Numeric number of bins to trim from start and end of segments
#'
#' @export
writeSegTable <- function(object, filename, trimLength=1, minLength=NULL){
  # minLength is ignored, keep for backward compatibility
  
  gr = getSegTable(object, trimLength=trimLength)
  
  df = do.call(rbind, lapply(gr, function(gr1){return(
    dplyr::bind_cols(
      data.frame(seqnames=seqnames(gr1),
                 start=start(gr1),
                 end=end(gr1),
                 width=width(gr1)),
      data.frame(gr1@elementMetadata)))}))
  
  write.table(df, file=filename, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
}

#' exportSignals 
#'
#' \code{exportSignals} export QDNAseq object to a tsv file, for processing with signals
#'
#' @param object QDNAseq object
#' @param filename filename to write file to
#'
#' @export
exportSignals <- function(object, file, ...){
  require(Biobase)
  
  valid = binsToUseInternal(object)
  object = object[valid]
  feature <- featureNames(object)
  chromosome <- fData(object)$chromosome
  start <- fData(object)$start
  end <- fData(object)$end
  dat <- assayDataElement(object, "copynumber")
  calls <- assayDataElement(object, "calls")
  name = pData(object)[["name"]]
  
  out <- data.frame(chr=chromosome, start=as.integer(start),
                    end=as.integer(end), state=as.integer(c(dat)), copy=c(calls), cell_id=rep(name, each=dim(dat)[[1]]))
  write.table(out, file=file,
              quote=FALSE, sep="\t", na="", row.names=FALSE, ...)
}

