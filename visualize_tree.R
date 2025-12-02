library("reshape")
library(dplyr)
library(ComplexHeatmap)
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/helper_functions.R")
source("OneDrive - CRUK Cambridge Institute/CINSignatureDiscovery/CIN_Compendium_Discovery/scripts/main_functions.R")
source("OneDrive - CRUK Cambridge Institute/scCNV_analysis/scUnique-main/R/visualize.R")


sample="SLX-23965-50-qc_100"
sample="SLX-24077-50-qc_100"
sample="SLX-24130-50-qc_100"
sample="SLX-24130-50-qc_100"


sample="PEO50_qc_100"
sample="PEO1-evol-15_100"
sample="SLX-25393-15-qc_100"

sample="SLX-25394_D7-15-qc_100"
tree <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".medicc_tree.RDS"))
cn <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".finalCN.RDS"))
tree_profiles <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".tree_profiles.RDS"))
HEATMAP_PATH <- paste0("../../Volumes/Fl/sc_analysis/all_samples_23july2024/heatmap_dendrogram_",sample,".pdf")


plot(tree)
library(dendextend)  # Install with install.packages("dendextend") if needed
library(ape)
tree_colored <- color_branches(tree, k = 4)  # Color branches into 4 clusters
plot(tree_colored)


labels_old <- labels(tree)



# Function to rename labels
rename_labels <- function(label) {
  # Replace SINCEL-xxx with corresponding PEO
  label <- gsub("^SINCEL-219", "PEO4", label)
  label <- gsub("^SINCEL-212", "PEO6", label)
  label <- gsub("^SINCEL-211", "PEO1", label)
  label <- gsub("^SINCEL-249", "PEO1 M", label)
  label <- gsub("^SINCEL_252", "PEO1 S", label)
  
  
  label <- gsub("^SINCEL-281", "G3", label)
  label <- gsub("^SINCEL-282", "D5", label)
  label <- gsub("^SINCEL-283", "D7", label)
  label <- gsub("^SINCEL-277", "D11", label)
  label <- gsub("^SINCEL-280", "C9", label)
  label <- gsub("^SINCEL-279", "B4", label)
  label <- gsub("^SINCEL-274", "PEO1 Parent", label)
  
  # Extract last part after the last underscore "_"
  suffix <- sub(".*_", "", label)
  
  # Keep only the PEO prefix and the last part
  prefix <- sub("_.*", "", label)  # Extract PEOx prefix
  new_label <- paste0(prefix, "_", suffix)
  
  return(new_label)
}

# Apply function to labels
labels_new <- sapply(labels_old, rename_labels)

labels(tree) <- labels_new

# Create a vector for label colors: set default to black
label_colors <- rep("black", length(labels(tree)))

# Assign specific colors for PEO1, PEO4, PEO6
label_colors[grepl("^G3", labels(tree))] <- "#FF5733"  # Fiery Red
label_colors[grepl("^D5", labels(tree))] <- "#33FF57"  # Vibrant Green
label_colors[grepl("^D7", labels(tree))] <- "#FF8C33"  # Bright Blue
label_colors[grepl("^D11", labels(tree))] <- "#FF33A1"  # Hot Pink
label_colors[grepl("^C9", labels(tree))] <- "#33FFF5"  # Aqua Blue
label_colors[grepl("^B4", labels(tree))] <- "#F5FF33"  # Sunny Yellow
label_colors[grepl("^PEO1 Parent", labels(tree))] <- "#3357FF"  # Orange


labels_colors(tree) <- label_colors
# Color the labels in the dendrogram using color_labels
tree <- tree %>% set("labels_col", label_colors)
# Open a PDF device
pdf(paste0("../../Volumes/Fl/sc_analysis/all_samples_23july2024/dendrogram_colored_labels_",sample,".pdf"), width = 25, height = 10)
# Plot the colored dendrogram
plot(tree, main = "MEDDIC2 tree", cex = 0.3)

# Add a legend at the bottom right
legend("bottomright", legend = c("G3", "D5", "D7", "D11", "C9", "B4", "PEO1 Parent"),
       fill = c("#FF5733", "#33FF57", "#FF8C33", "#FF33A1", "#33FFF5", "#F5FF33", "#3357FF"),
       title = "Labels")

# Close the PDF device to save the file
dev.off()



object = cn
file =HEATMAP_PATH
cluster_rows = prune(tree, "diploid")
cluster_rows = prune(cluster_rows, "diploid_diploid")
row_split=NULL
cutoff=10
show_unobserved_states=TRUE
har=NULL
useCopynumber=TRUE
show_cell_names=TRUE
abbreviate_cell_names=FALSE
show_chromosome_names=FALSE
use_cell_names=labels(cluster_rows)

fontsize_row=9
fontsize_col=9 
fontsize_chr=12
fontsize_leg_title=18
fontsize_leg_label=14 
row_title_gp=13
column_title_gp=13
column_title=""
bottom_annotation=NULL
show_heatmap_legend=TRUE
scale_copynumber=1.0
raster_device="tiff"
colorMap="MSK"
row_title=paste0("tree")

  
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
  
  row_labels <- labels(cluster_rows)
  row_order <- match(row_labels, rownames(Y_plot))
  Y_plot_ordered <- Y_plot[row_order, ]
  ht = Heatmap(Y_plot_ordered, name = "CNheatmap", cluster_columns = c_cols, 
               row_split=row_split, cluster_rows = cluster_rows,
               col = colors, na_col = NA_color,
               column_split = chr, bottom_annotation = hac, left_annotation = NULL, right_annotation = har, 
               row_title_side = "left", row_dend_side = "left",row_order = order.dendrogram(cluster_rows), row_names_side = "left", row_title_gp = gpar(fontsize = row_title_gp),
               column_title=column_title, column_title_side = "bottom",  column_title_gp = gpar(fontsize = column_title_gp), show_column_names = show_chromosome_names,
               row_names_gp = gpar(fontsize = fontsize_row), column_names_gp = gpar(fontsize = fontsize_col), show_row_names = show_cell_names, row_labels = rw[1:15],
               use_raster = TRUE, raster_device = raster_device, raster_quality = 10,
               show_heatmap_legend=show_heatmap_legend,
               heatmap_legend_param = list(title = "Copy number", title_position = "lefttop-rot",
                                           at=names, labels=labels,  grid_height = unit(0.5, "cm"),
                                           title_gp = gpar(fontsize=fontsize_leg_title),
                                           labels_gp = gpar(fontsize=fontsize_leg_label)), row_title=row_title,row_dend_width = unit(5, "cm"))
  
  
  if(!is.null(file)){
    pdf(file, width = 20, height = 6)
    draw(ht)
    dev.off()
  }else{
    # draw(ht)
    return(ht)
  }
  









