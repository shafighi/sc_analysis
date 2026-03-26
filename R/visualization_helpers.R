
# plot copy number heatmap
plotCopynumberHeatmap <- function(object, file=NULL, 
                                          cluster_rows = FALSE, row_split=NULL,
                                          cutoff=10, show_unobserved_states=TRUE,
                                          har=NULL, useCopynumber=TRUE,
                                          show_cell_names=TRUE, abbreviate_cell_names=TRUE, show_chromosome_names=TRUE,
                                          use_cell_names="", fontsize_row=9, fontsize_col=9, fontsize_chr=12,
                                          fontsize_leg_title=18, fontsize_leg_label=14,
                                          row_title_gp=13, column_title_gp=13,
                                          column_title=NULL, bottom_annotation=NULL,
                                          show_heatmap_legend=TRUE, scale_copynumber=1.0,
                                          raster_device="tiff",
                                          colorMap="MSK", row_title, ...){

  if("character" %in% class(cluster_rows) && cluster_rows == "ploidy"){
    ploidyOrder = sort(apply(object@assayData$copynumber, 2, mean, na.rm=TRUE), index.return=TRUE)[["ix"]]
    object = object[, ploidyOrder]
    cluster_rows=FALSE
  }
  
  if("list" %in% class(cluster_rows)){
    n_clusters = unlist(lapply(cluster_rows, function(x) length(x) != 0))
    n_each = unlist(lapply(cluster_rows, function(x) length(x)))
    sortOrder = unlist(cluster_rows)
    row_split = rep(seq(1, sum(n_clusters)), times=n_each[n_clusters])
    object = object[, sortOrder]
    cluster_rows = FALSE
  }
    
  # DEBUGGING purposes - default parameters
  # file = "~/plot-outliers_30.pdf"
  # file = NULL
  # cutoff=10
  # row_split=NULL
  # cluster_rows=NULL
  # show_cell_names=TRUE
  # abbreviate_cell_names=TRUE
  # show_unobserved_states=TRUE
  # colorMap="MSK"
  # showUnique=TRUE
  # har=NULL
  
  require(ComplexHeatmap, quietly=TRUE)
  require(circlize, quietly=TRUE)
  require(viridis, quietly=TRUE)
  require(cluster, quietly=TRUE)
  require(RColorBrewer, quietly=TRUE)
  require(stringr, quietly=TRUE)
  require(QDNAseq, quietly=TRUE)
  require(Biobase, quietly=TRUE)
  # require(magick, quietly=TRUE)
  if(require(fastcluster, quietly=TRUE)){
    ht_opt("fast_hclust" = TRUE)
  }else{
    ht_opt("fast_hclust" = FALSE)
  }


  stopifnot("QDNAseqCopyNumbers" %in% class(object))
  # stopifnot(!is.null(assayDataElement(object, "segmented")))
  stopifnot(grepl("\\.pdf$", file) || is.null(file))
  if(!is.null(row_split)){
    stopifnot(is.integer(row_split) || is.factor(row_split))
    stopifnot("dendrogram" %in% class(cluster_rows) || "logical" %in% class(cluster_rows))
  }

  if(useCopynumber){
    Y = t(round(Biobase::assayDataElement(object, "copynumber"), digits=0)) * scale_copynumber
    Y[!(abs(Y%%1) < 1e-3)] = NA
  }else{
    warning("Using segmentation slot")
    Y = t(round(Biobase::assayDataElement(object, "segmented"), digits=0)) * scale_copynumber
    Y[!(abs(Y%%1) < 1e-3)] = NA
  }
  

  # Prepare 
  if(abbreviate_cell_names){
    rw = substr(stringr::str_extract(rownames(Y), "_00\\d{4}"), 2, 7)
    rownames(Y) = ifelse(is.na(rw), NA, as.character(as.numeric(rownames(rw))))
  }else{
    rw = use_cell_names
  }
  
  chr = factor(sapply(strsplit(colnames(Y), ":"), "[[", 1), levels=c(as.character(seq(1, 22)), "X", "Y"))

  # filter out NaN-only chromosomes
  subsets = lapply(split(seq(1, ncol(Y)), chr), function(m, ind) m[,ind], m = Y)
  subchr = split(chr, chr)

  list_of_nan_chromosomes = lapply(subsets, function(x){return(all(is.na(x)))})
  subsets = subsets[which(!unlist(list_of_nan_chromosomes))]
  subchr = subchr[which(!unlist(list_of_nan_chromosomes))]
  Y <- do.call("cbind", subsets)
  chr <- unlist(subchr)

  if (is.null(column_title)) {
    column_title <- as.character(levels(chr)[levels(chr) %in% as.character(chr)])
  }

  colnames(Y) = NULL
  ## create clustering of cells ====
  if(is.null(cluster_rows)){
      # cluster_rows = cluster::diana(Y)
      row_split = NULL
  }
  c_cols = FALSE # cluster::agnes(t(Y))
  
  
  ## plot heatmap ====
  if(any(Y > cutoff, na.rm=TRUE)){
    warning("Cuting off high copy numbers!")
  }
  out_of_range = any(Y > cutoff, na.rm=TRUE)
  Y_plot = ifelse(Y <= cutoff, Y, cutoff+1)
  
  ## Color scheme ####
  # colset = c("#999999", "#56B4E9", "#43BF71FF", "#0072B2", "#D55E00", "#414487FF", "#E69F00", "#440154FF", "#F0E442", "#CC79A7", "#CC79A7", "#CC79A7")
  if(colorMap == "MSK"){
    # custom color scheme used by cellranger and at MSK
    colset = c("#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A","#FC8D59","#E34A33","#B30000","#980043","#DD1C77","#DF65B0","#C994C7","#D4B9DA")
    NA_color = "#000000"
    UQ_color = "#FFFFFF"
  }else{
    colset = c("#000000", 
               "#56B4E9", "#009E73", #"#43BF71FF",
               #Middle greens
               # "#43BF71FF", "#009E73FF", 
               # "#21908CFF", "#1E9B8AFF", "#21A685FF", "#2BB07FFF", "#3BBB75FF", "#51C56AFF", "#43BF71FF", "#009E73FF", 
               "#0072B2", "#D55E00", 
               "#0D0887FF", "#E69F00", #"#FBBC21FF", #"#F48849FF",
               "#5402A3FF", "#F0E442", #"#F4E258FF", #"#FEBC2AFF",
               "#8B0AA5FF", "#FCFFA4FF", #"#F0F921FF",
               # "#B93289FF", "#DB5C68FF",
               # "#DA4E3CFF", "#47039FFF",
               # "#F06F20FF", "#7301A8FF",
               # "#FB9706FF", "#9C179EFF", 
               # "#FAC52BFF", "#BD3786FF",
               "#CC79A7")
    NA_color = "#D3D3D3"
    UQ_color = "#FFFFFF"
  }
  ## end color scheme
  
  ## Legend - copy number states ####
  cutoff_observed = min(cutoff, max(Y_plot, na.rm=TRUE))
  if(show_unobserved_states){
    innames = c(sapply(0:(cutoff_observed), as.character))
  }else{
    innames = names(table(Y_plot))
  }

  if(out_of_range){
    if(cutoff <= 9){
      labels = c(paste0(" ", sapply(0:cutoff_observed, as.character)), paste0(">", cutoff))
    }else{
      labels = c(paste0("  ", sapply(0:9, as.character)), paste0(" ",sapply(10:cutoff_observed, as.character)), paste0(">", cutoff))
    }
    names = c(innames, as.character(cutoff_observed+1))
  }else{
    labels = c(paste0("", sapply(0:cutoff_observed, as.character)))
    names = innames
  }
  ## end of legend part 
  
  colors = structure(c(UQ_color, colset), names = c("-1", names))
  
  if(is.null(bottom_annotation)){
    lst = as.character(unique(chr))
    #hac = columnAnnotation(foo = ComplexHeatmap::anno_block(setNames(as.list(1:length(lst)), lst), gp = gpar(fontsize = fontsize_chr, col="white"), labels_gp = gpar(fontsize = fontsize_chr-2)))
    hac = columnAnnotation(
      foo = ComplexHeatmap::anno_block(
        #setNames(as.list(1:length(lst)), lst),  # Make sure lst is a character vector or factor
        labels=lst,
        gp = gpar(fill = NA, col = "white"),  # 'col' controls the color of borders; ensure it's valid here
        labels_gp = gpar(col = "black",fontsize = fontsize_chr - 2)  # Adjust the font size for labels
      )
    )
    hac = NULL
    }else{
    hac = bottom_annotation
  }
  if(is.logical(bottom_annotation) && !bottom_annotation){
    hac = NULL
  }

  ht = Heatmap(Y_plot, name = "CNheatmap", cluster_columns = c_cols, 
               row_split=row_split, cluster_rows = cluster_rows,
               col = colors, na_col = NA_color,
               column_split = chr, bottom_annotation = hac, left_annotation = NULL, right_annotation = har, 
               row_title_side = "left", row_dend_side = "left", row_names_side = "left", row_title_gp = gpar(fontsize = row_title_gp),
               column_title=column_title, column_title_side = "bottom",  column_title_gp = gpar(fontsize = column_title_gp), show_column_names = show_chromosome_names,
               row_names_gp = gpar(fontsize = fontsize_row), column_names_gp = gpar(fontsize = fontsize_col), show_row_names = show_cell_names, row_labels = rw,
               use_raster = TRUE, raster_device = raster_device, raster_quality = 10,
               show_heatmap_legend=show_heatmap_legend,
               heatmap_legend_param = list(title = "Copy number", title_position = "lefttop-rot",
                                           at=names, labels=labels,  grid_height = unit(0.5, "cm"),
                                           title_gp = gpar(fontsize=fontsize_leg_title),
                                           labels_gp = gpar(fontsize=fontsize_leg_label)), row_title=row_title,...)
  
  
  if(!is.null(file)){
    pdf(file, width = 11, height = 4)
    draw(ht)
    dev.off()
  }else{
    # draw(ht)
    return(ht)
  }
  
}


plotUniqueEvents <- function(object, background, dt, indices=NULL, difference_only=FALSE, show_cell_names=FALSE, horizontal=FALSE, ...){
  require(ggpubr)
  require(grid)
  
  valid = binsToUseInternal(object)
  if(is.null(indices)){
    p1 = plotCopynumberHeatmap(object[, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names, ...)
    uniqueEventEncoding = matrix(-1, nrow=dim(object)[[1]], ncol=dim(object)[[2]])
    uniqueEventIndices = object@assayData$copynumber != background & !is.na(object@assayData$copynumber)
    uniqueEventEncoding[uniqueEventIndices] = object@assayData$copynumber[uniqueEventIndices]
    print(prop.table(table(uniqueEventIndices)))
    
    profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object),
                        copynumber = uniqueEventEncoding,
                        phenodata = Biobase::phenoData(object))
  }else{
    p1 = plotCopynumberHeatmap(object[indices, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names, ...)
    uniqueEventEncoding = matrix(-1, nrow=length(indices), ncol=dim(object)[[2]])
    uniqueEventIndices = object@assayData$copynumber[indices,] != background[indices,] & !is.na(object@assayData$copynumber[indices,])
    uniqueEventEncoding[uniqueEventIndices] = object[indices,]@assayData$copynumber[uniqueEventIndices]
    print(prop.table(table(uniqueEventIndices)))
    
    profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object[indices]),
                        copynumber = uniqueEventEncoding,
                        phenodata = Biobase::phenoData(object[indices]))
  }
  
  
  p2 = plotCopynumberHeatmap(profileUnique[, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names, ...)
  
  gb1 = grid::grid.grabExpr(draw(p1))
  gb2 = grid::grid.grabExpr(draw(p2))
  
  # grid.newpage()
  # pushViewport(viewport(y = 1, height = 0.5, just = "top"))
  # grid.draw(gb1)
  # popViewport()
  # pushViewport(viewport(x = 0, y = 0, height = 0.5, just = c("bottom")))
  # grid.draw(gb2)
  # popViewport()
  
  # gb1
  # is.grob(gb1)
  # pushViewport(viewport(width = 0.5, height = 0.5))
  # grid.draw(gb1)
  
  gtb1 <- gtable::gtable_matrix("hm_gtbl1", matrix(list(gb1)), unit(1, "null"), unit(1, "null"))
  gtb2 <- gtable::gtable_matrix("hm_gtbl2", matrix(list(gb2)), unit(1, "null"), unit(1, "null"))
  if(horizontal){
    pp = cowplot::plot_grid(gtb1, gtb2, ncol=2)
  }else{
    pp = cowplot::plot_grid(gtb1, gtb2, nrow=2)
  }
  # # hm_gtable <- gtable_matrix("hm_gtbl", matrix(list(gb_heatmap)), unit(1, "null"), unit(1, "null"))
  # pp = ggpubr::ggarrange(gtb1 + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), 
  #                gtb2, nrow=2, ncol=1, common.legend=TRUE)
  
  if(difference_only){
    return(p2)
  }
  
  return(pp)
}

plotBackgroundProfiles <- function(object, background, dt, indices=NULL, background_only=FALSE, show_cell_names=FALSE, ...){
  valid = binsToUseInternal(object)
  if(is.null(indices)){
    p1 = plotCopynumberHeatmap(object[, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names, ...)
    profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object),
                        copynumber = background,
                        phenodata = Biobase::phenoData(object))
  }else{
    p1 = plotCopynumberHeatmap(object[indices, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names, ...)
    profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object[indices,]),
                        copynumber = background[indices,],
                        phenodata = Biobase::phenoData(object[indices,]))
  }
  p2 = plotCopynumberHeatmap(profileUnique[, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names, ...)
  gb1 = grid::grid.grabExpr(draw(p1))
  gb2 = grid::grid.grabExpr(draw(p2))
  
  gtb1 <- gtable::gtable_matrix("hm_gtbl1", matrix(list(gb1)), unit(1, "null"), unit(1, "null"))
  gtb2 <- gtable::gtable_matrix("hm_gtbl2", matrix(list(gb2)), unit(1, "null"), unit(1, "null"))
  pb = cowplot::plot_grid(gtb1, gtb2, nrow=2)
  
  if(background_only){
    return(p2)
  }
  
  return(pb)
}


plotUniqueOverBackground <- function(object, background, dt, indices=NULL, show_cell_names=FALSE, ...){
  require(ggpubr)
  require(grid)
  
  valid = binsToUseInternal(object)
  if(is.null(indices)){
    uniqueEventEncoding = matrix(-1, nrow=dim(object)[[1]], ncol=dim(object)[[2]])
    uniqueEventIndices = object@assayData$copynumber != background & !is.na(object@assayData$copynumber)
    uniqueEventEncoding[uniqueEventIndices] = object@assayData$copynumber[uniqueEventIndices]
    print(prop.table(table(uniqueEventIndices)))
    
    profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object),
                        copynumber = uniqueEventEncoding,
                        phenodata = Biobase::phenoData(object))
    profileBackground = new("QDNAseqCopyNumbers",
                            bins = Biobase::featureData(object),
                            copynumber = background,
                            phenodata = Biobase::phenoData(object))
  }else{
    uniqueEventEncoding = matrix(-1, nrow=length(indices), ncol=dim(object)[[2]])
    uniqueEventIndices = object@assayData$copynumber[indices,] != background[indices,] & !is.na(object@assayData$copynumber[indices,])
    uniqueEventEncoding[uniqueEventIndices] = object[indices,]@assayData$copynumber[uniqueEventIndices]
    print(table(uniqueEventIndices))
    
    profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object[indices]),
                        copynumber = uniqueEventEncoding,
                        phenodata = Biobase::phenoData(object[indices]))
    
    profileBackground = new("QDNAseqCopyNumbers",
                            bins = Biobase::featureData(object[indices,]),
                            copynumber = background[indices,],
                            phenodata = Biobase::phenoData(object[indices,]))
  }
  p2 = plotCopynumberHeatmap(profileUnique[, labels(dt)], cluster_rows=dt, show_cell_names=show_cell_names, ...)
  p1 = plotCopynumberHeatmap(profileBackground[, labels(dt)], cluster_rows=dt, show_cell_names=show_cell_names, ...)
  
  gb1 = grid::grid.grabExpr(draw(p1))
  gb2 = grid::grid.grabExpr(draw(p2))
  
  gtb1 <- gtable::gtable_matrix("hm_gtbl1", matrix(list(gb1)), unit(1, "null"), unit(1, "null"))
  gtb2 <- gtable::gtable_matrix("hm_gtbl2", matrix(list(gb2)), unit(1, "null"), unit(1, "null"))
  pp = cowplot::plot_grid(gtb1, gtb2, nrow=2)
  # # hm_gtable <- gtable_matrix("hm_gtbl", matrix(list(gb_heatmap)), unit(1, "null"), unit(1, "null"))
  # pp = ggpubr::ggarrange(gtb1 + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), 
  #                gtb2, nrow=2, ncol=1, common.legend=TRUE)
  
  return(pp)
}


## visualize clones
plotClones <- function(cprofiles, cTree, object, K_set=NULL, include=NULL, exclude=NULL){
  cprofiles[["0"]] = NULL
  csize = table(cTree[cTree != 0])
  valid = binsToUseInternal(object)

  if(is.null(K_set)){
    K_set = 1:(length(cprofiles))
  }
  stopifnot(length(cprofiles) == length(csize))
  
  clone_copynumber = matrix(data = NA, ncol=dim(object)[[1]], nrow=sum(csize[K_set]))

  offset = 1
  for(k in K_set){
    clone_size = csize[[k]]
    clone_copynumber[offset:(offset+clone_size-1),valid] = matrix(rep(t(cprofiles[[k]]), times=clone_size), nrow=clone_size, byrow=TRUE)
    offset = offset + clone_size
  }

  profile = new("QDNAseqCopyNumbers",
                bins = Biobase::featureData(object),
                copynumber = t(clone_copynumber),
                phenodata = data.frame(cluster=paste0("Cluster ", base::rep(K_set, times=csize[K_set]))))
  Biobase::assayDataElement(profile, "segmented") = t(clone_copynumber)
  
  cn = plotCopynumberHeatmap(selectChromosomes(profile, include=include, exclude=exclude), cluster_rows = TRUE, show_cell_names = FALSE, 
                    row_split = as.factor(paste0("Cluster ", str_pad(base::rep(K_set, times=csize[K_set]), ceiling(log10(length(K_set))), pad = "0"))),
                    cluster_row_slices=FALSE, show_row_dend = FALSE, 
                    row_title = paste0(LETTERS[K_set], " - ", csize[K_set]), 
                    row_title_rot = 0)

  return(cn)
}


## Observe clusters
### Visualize clusters - size proportional to elements
#plotCopynumberCluster <- function(object, membership, cluster_profiles, K_set=NULL, include=NULL, exclude=NULL){
#
#  if(is.null(K_set)){
#    K_set = 1:(dim(cluster_profiles)[[2]])
#  }
#  stopifnot(length(membership) == dim(cluster_profiles)[[2]])
#  stopifnot(dim(object)[[1]] == dim(cluster_profiles)[[1]])
#  stopifnot(dim(object)[[2]] == sum(unlist(lapply(membership, length))))
#  
#  valid = binsToUseInternal(object)
#  
#  cluster_copynumber = matrix(data = NA, nrow=sum(unlist(lapply(K_set, function(x) length(membership[[x]])))),
#                              ncol=dim(object)[[1]])
#  
#  zeros = which(lapply(membership, length) == 0 & seq(1,(dim(cluster_profiles)[[2]])) %in% K_set)
#  K_prime = length(K_set) - length(zeros)
#  K_valid = setdiff(K_set, zeros)
#
#  offset = 1
#  for(k in K_valid){
#    members = membership[[k]]
#    cluster_copynumber[offset:(offset+length(members)-1), valid] = matrix(rep(t(cluster_profiles[valid, k]), times=length(members)), nrow=length(members), byrow=TRUE) 
#    offset = offset + length(members)
#  }
#  cluster_size = unlist(lapply(membership, length))
#
#  profile = new("QDNAseqCopyNumbers",
#                bins = Biobase::featureData(object),
#                copynumber = t(cluster_copynumber),
#                phenodata = data.frame(cluster=paste0("Cluster ", base::rep(seq(1, K_prime), times=cluster_size[K_valid]))))
#  Biobase::assayDataElement(profile, "segmented") = t(cluster_copynumber)
#  
#  cn = plotCopynumberHeatmap(selectChromosomes(profile, include=include, exclude=exclude), cluster_rows = TRUE, show_cell_names = FALSE, 
#                    row_split = as.factor(paste0("Cluster ", str_pad(base::rep(seq(1, K_prime), times=cluster_size[K_valid]), ceiling(log10(K_prime)), pad = "0"))),
#                    cluster_row_slices=FALSE, show_row_dend = FALSE, 
#                    row_title = paste0(LETTERS[1:K_prime], " - ", cluster_size[K_valid]), 
#                    row_title_rot = 0)
#
#  return(cn)
#}
#
### Visualize  unique events
#plotCopynumberUnique <- function(object, membership, unique_matrix, 
#                                 K_set=NULL, abbreviate_cell_names = TRUE, include=NULL, exclude=NULL){
#
#  stopifnot(dim(object) == dim(unique_matrix))
#  if(is.null(K_set)){
#    K_set = 1:length(membership)
#  }
#  
#  valid = binsToUseInternal(object)
#
#  # ## part to see unique copy numbers
#  # per_cell_copynumber = Biobase::assayDataElement(object, "copynumber")
#  # 
#  # n_rows = sum(unlist(lapply(K_set, function(x) length(membership[[x]])))) + 
#  #              sum(unlist(lapply(K_set, function(x) length(membership[[x]]) > 0)))
#  sorted_copynumber = matrix(data = NA, nrow=dim(object)[[1]], ncol=dim(object)[[2]])
#  
#  K_set = K_set[unlist(lapply(K_set, function(x) length(membership[[x]]) > 0))]
#  offset = 1
#  row_names = c()
#  row_split = c()
#  row_gaps = c()
#  counter = 0
#  for(k in K_set){
#    # print(k)
#    members = membership[[k]]
#    if(length(members) == 0){
#      next
#    }
#    
#    # set cluster profile - not needed anymore
#    # sorted_copynumber[valid, (offset+counter-1)] = cluster profile
#    # row_names = c(row_names, paste0("Cluster ", k))
#    # row_split = c(row_split, paste0("Cluster ", str_pad(k, ceiling(log10(length(K_set))), pad="0")))
#    # row_gaps = c(row_gaps, 1)
#    
#    # print("---")
#    # print(offset+k-1)
#    # print((offset+k))
#    # print((offset+length(members)+k-1))
#    
#    # # add individual profiles
#    sorted_copynumber[valid, (offset+counter):(offset+length(members)+counter-1)] = unique_matrix[valid, members]
#    if(abbreviate_cell_names){
#      row_names = c(row_names, stringr::str_extract(colnames(object)[members], "00\\d{4}"))
#    }else{
#      row_names = c(row_names, colnames(object)[members])
#    }
#    row_split = c(row_split, rep(paste0("Members ", str_pad(k, ceiling(log10(length(K_set))), pad="0")), times=length(members)))
#    # 0 sets row_gap to solid line, other values add whitespace
#    row_gaps = c(row_gaps, 3)
#    offset = offset + length(members)
#    # counter = counter + 1
#  }
#  cluster_size = unlist(lapply(membership, length))
#  
#  profile = new("QDNAseqCopyNumbers",
#                bins = Biobase::featureData(object),
#                copynumber = sorted_copynumber,
#                phenodata = data.frame(cluster = row_names))
#  Biobase::assayDataElement(profile, "segmented") = sorted_copynumber
#
#  cn = plotCopynumberHeatmap(selectChromosomes(profile, include=include, exclude=exclude), 
#                            cluster_rows = TRUE, show_cell_names = TRUE, 
#                            abbreviate_cell_names=FALSE, use_cell_names=row_names,
#                            row_split = factor(row_split, levels=unique(row_split)),
#                            cluster_row_slices=FALSE, show_row_dend = FALSE,
#                            row_gap = unit(row_gaps, "mm"),  border = FALSE, #border adds solid lines
#                            column_gap = unit(.9, "mm"),
#                            row_title = "", row_title_rot = 0)
#
#  return(cn)
#}
#
#



pre_scUnique_filtering <- function(df,OUTPUT,mapd_cutoff=1.5,mapd_density_control= TRUE,gini_norm_cutoff= 1.5,alpha_cutoff=2,alpha_hard_cutoff =0.02,gini_density_control=TRUE,rpc_cutoff=15){
  filtered_rpc <- df %>% dplyr::filter(!replicating, rpc >= rpc_cutoff)
  
  filtered_mapd <- qc_mapd(filtered_rpc,cutoff_percentile=0.05, cutoff_value = mapd_cutoff, symmetric=FALSE,
                           densityControl=mapd_density_control, densityCutoff=0.10)
  filtered_mapd$dmapd.outlier <- filtered_mapd$dmapd.outlier[,"outlier"]
  filtered_gini <- qc_gini_norm(filtered_mapd, cutoff_value = gini_norm_cutoff, 
                                densityControl = gini_density_control)
  filtered_gini$dgini.outlier <- filtered_gini$dgini.outlier[,"outlier"]
  iq = qc_alpha(filtered_gini, cutoff_value=alpha_cutoff, hard_cutoff = alpha_hard_cutoff)
  iq$outlier = iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier
  iq$outlier <- iq$outlier[,"outlier"]
  print(iq$dmapd.outlier )
  fig_outliers = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=alpha, colour=outlier), 
                                                      size=0.8, alpha=1.0) + facet_wrap(~SLX, scales = "free_x")+
    theme_pubclean()
  fig_outliers
  ggsave(paste0(OUTPUT,"/cell_qc.jpg"), width = 3, height = 4)
  return(iq)
}

general_analysis_plots <- function(df,iq,OUTPUT){
  sample_name = df$UID[1]
  
  pdf(paste0(OUTPUT,"/normal_dist_alpha.pdf"), width = 3.5, height = 3, pointsize=7)
  cutoff_value <- 1.5
  x <- df$hmm.alpha
  hist(x, breaks = 60, col = "lightblue", main = "Filtering Based on Alpha Distribution",xlab=expression(alpha))
  abline(v = mean(x), col = "blue", lwd = 1)
  text(mean(x), par("usr")[4], labels = sprintf("mean(x).   %.2f", mean(x)),
       pos = 2, offset = 0.5, col = "blue", srt = 90, cex = 0.7)
  condition_region <- x [log(x) > mean(log(x)) + (cutoff_value * sd(log(x)))]
  abline(v = exp(mean(log(x)) + (cutoff_value * sd(log(x)))), col = "red", lwd = 1)
  text(exp(mean(log(x)) + (cutoff_value * sd(log(x)))), par("usr")[4], labels = sprintf("exp(mean(log(x)) + (cutoff_value * sd(log(x)))).   %.2f", exp(mean(log(x)) + (cutoff_value * sd(log(x))))),
       pos = 2, offset = 0.5, col = "red", srt = 90, cex = 0.7)
  abline(v = exp(mean(log(x)) + (sd(log(x)))), col = "green", lwd = 1)
  text(exp(mean(log(x)) + (sd(log(x)))), par("usr")[4], labels = sprintf("exp(mean(log(x)) + (sd(log(x)))).   %.2f", exp(mean(log(x)) + (sd(log(x))))),
       pos = 2, offset = 0.5, col = "green", srt = 90, cex = 0.7)
  hist(condition_region, breaks = 30, main = "Distribution with Filtering Condition", col = rgb(1, 0, 0, 0.5), add = TRUE,alpha=0.5)
  legend("topright", legend = c("Data", "Mean", "Filtering Region"),
         fill = c("lightblue", "blue", rgb(1, 0, 0, 0.5)), border = NA)
  dev.off()
  
  
  pdf(paste0(OUTPUT,"/total.reads_all.pdf"), width = 2, height = 3, pointsize=7)  # Adjust width and height as needed
  
  # Change the size of the font and figure
  par(
    cex.lab = 1,      # Set the size of axis labels
    cex.main = 1      # Set the size of the main title
  )
  
  # Create a histogram
  hist(df$total.reads, 
       main = "Distribution of Total Reads\n(all cells)",  # Title of the plot
       xlab = "Total Reads",                  # Label for x-axis
       ylab = "Probability",                    # Label for y-axis
       col = "skyblue",                       # Color of bars
       border = "black",                      # Border color of bars
       xlim = c(0, max(df$total.reads)),     # X-axis limits
       breaks = 40,                           # Number of bins
       prob = TRUE)                           # Density scale
  
  # Add a smoothed density curve
  lines(density(df$total.reads), col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Histogram", "Density"), fill = c("skyblue", "red"), bty = "n")
  
  dev.off()
  
  
  
  
  pdf(paste0(OUTPUT,"/total.reads_non_replicating.pdf"), width = 2, height = 3, pointsize=7)  # Adjust width and height as needed
  
  # Change the size of the font and figure
  par(
    cex.lab = 1,      # Set the size of axis labels
    cex.main = 1      # Set the size of the main title
  )
  
  # Create a histogram
  hist(iq$total.reads, 
       main = "Distribution of Total Reads\n(non-replicating)",  # Title of the plot
       xlab = "Total Reads",                  # Label for x-axis
       ylab = "Probability",                    # Label for y-axis
       col = "skyblue",                       # Color of bars
       border = "black",                      # Border color of bars
       xlim = c(0, max(iq$total.reads)),     # X-axis limits
       breaks = 40,                           # Number of bins
       prob = TRUE)                           # Density scale
  
  # Add a smoothed density curve
  lines(density(df$total.reads), col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Histogram", "Density"), fill = c("skyblue", "red"), bty = "n")
  
  dev.off()
  
  
  
  
  pdf(paste0(OUTPUT,"/total.reads_non_outliers.pdf"), width = 2, height = 3, pointsize=7)  # Adjust width and height as needed
  
  # Change the size of the font and figure
  par(
    cex.lab = 1,      # Set the size of axis labels
    cex.main = 1      # Set the size of the main title
  )
  
  # Create a histogram
  hist(iq$total.reads[iq$outlier==FALSE], 
       main = "Distribution of Total Reads\n(non-replicating)",  # Title of the plot
       xlab = "Total Reads",                  # Label for x-axis
       ylab = "Probability",                    # Label for y-axis
       col = "skyblue",                       # Color of bars
       border = "black",                      # Border color of bars
       xlim = c(0, max(iq$total.reads)),     # X-axis limits
       breaks = 40,                           # Number of bins
       prob = TRUE)                           # Density scale
  
  # Add a smoothed density curve
  lines(density(df$total.reads), col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Histogram", "Density"), fill = c("skyblue", "red"), bty = "n")
  
  dev.off()
  
  
  
  
  
  fig_cellcycle = ggplot(data = df) + geom_quasirandom(aes(x=SLX, y=cycling_activity,
                                                           color=replicating), size=0.8, alpha=1.0) +
    facet_wrap(~UID, scales = "free_x") + theme_pubclean()
  fig_cellcycle
  ggsave(paste0(OUTPUT,"/cellcycle.jpg"), width = 3, height = 4)
  
  
  print(iq)
  fig_outliers = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=alpha, colour=outlier), 
                                                      size=0.8, alpha=1.0) + facet_wrap(~SLX, scales = "free_x")+theme_pubclean()
  fig_outliers
  ggsave(paste0(OUTPUT,"/cell_qc.jpg"), width = 3, height = 4)
  
  
  
  
  fig_cellcycle = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=dmapd,color=dmapd.outlier), size=0.8, alpha=1.0) +
    facet_wrap(~UID, scales = "free_x") + theme_pubclean()
  fig_cellcycle
  ggsave(paste0(OUTPUT,"/dmapd_cutoff.jpg"), width = 3, height = 4)
  
  
  
  fig_cellcycle = ggplot(data = iq) + geom_quasirandom(aes(x=SLX, y=dgini,
                                                           color=dgini.outlier), size=0.8, alpha=1.0) +
    facet_wrap(~UID, scales = "free_x") + theme_pubclean()
  fig_cellcycle
  ggsave(paste0(OUTPUT,"/dgini_cutoff.jpg"), width = 3, height = 4)
  
  
  
  

  ggplot(iq, aes(used.reads, rpc, color = ploidy)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Used Reads",
         y = "RPC",
         color="Ploidy")
  ggsave(paste0(OUTPUT,"/rpc_usedreads_ploidy_replicatings.jpg"), width = 8, height = 6)
  
  
  
  ggplot(iq, aes(used.reads, rpc, color = as.factor(outlier))) +
    geom_point(size=1) +
    
    # Discrete color scale
    scale_color_brewer(type = "qual", palette = "Set1") +
    
    theme_light() +
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Used Reads",
         y = "RPC",
         color = "Outlier")
  ggsave(paste0(OUTPUT,"/rpc_usedreads_outliers_replicatings.jpg"), width = 8, height = 6)
  
  
  
  ggplot(iq, aes(alpha, rpc, color = ploidy)) +
    geom_point(size=1) + 
    scale_color_viridis_c() +
    theme_light() + 
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Alpha",
         y = "RPC",
         color="Ploidy")
  ggsave(paste0(OUTPUT,"/alpha_rpc_ploidy_replicatings.jpg"), width = 8, height = 6)
  
  
  ggplot(df, aes(rpc,alpha , color = as.factor(replicating))) +
    geom_point(size=1) +
    
    # Discrete color scale
    scale_color_brewer(type = "qual", palette = "Set1") +
    
    theme_light() +
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "RPC",
         y = "Alpha",
         color = "Replicating")
  ggsave(paste0(OUTPUT,"/alpha_rpc_replicating_all.jpg"), width = 8, height = 6)
  
  
  ggplot(df, aes(alpha,used.reads , color = as.factor(replicating))) +
    geom_point(size=1) +
    
    # Discrete color scale
    scale_color_brewer(type = "qual", palette = "Set1") +
    
    theme_light() +
    theme(aspect.ratio = 0.5,  
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_line(color="grey90")) +
    labs(title = sample_name,
         x = "Alpha",
         y = "Used reads",
         color = "Replicating")
  ggsave(paste0(OUTPUT,"/alpha_used_reads.jpg"), width = 8, height = 6)
  
  

   
  
  
}


plot_heatmap <- function(object,OUTPUT,SAMPLENAME,post_name,col_cluster){
  cn = object@assayData$copynumber
  df_no_na <- na.omit(cn)
  
  png(paste0(OUTPUT,"/",SAMPLENAME,"_",post_name,".png"), width = 7000, height = 2000, res = 300)
  
  pheatmap(t(df_no_na), 
           cluster_rows=TRUE,
           cluster_cols=col_cluster,
           display_numbers = FALSE, 
           number_color="black",
           show_rownames = FALSE,
           labels_col = "",
           fontsize_row = 7, 
           fontsize_col = 7 )
  
  dev.off()
}

plot_heatmap_t <- function(object,OUTPUT,SAMPLENAME,post_name,row_cluster){
  cn = object@assayData$copynumber
  df_no_na <- na.omit(cn)
  df_no_na = t(df_no_na)
  
  png(paste0(OUTPUT,"/",SAMPLENAME,"_",post_name,".png"), width = 7000, height = 3500, res = 300)
  
  pheatmap(df_no_na, 
           cluster_rows=row_cluster,
           cluster_cols=TRUE,
           display_numbers = FALSE, 
           number_color="black",labels_row = "")
  dev.off()
}
