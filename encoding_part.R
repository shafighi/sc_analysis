
sample <- "SLX-25394_B4-15-qc_100"

tree <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".tree.RDS"))
object <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".finalCN.RDS"))
background <- readRDS(paste0("../../Volumes/Fl/scUnique_results/",sample,"/",sample,".profiles_background_freq.RDS"))
dt <- as.dendrogram(tree)
require(ggpubr)
require(grid)

indices=NULL
difference_only=FALSE
show_cell_names=FALSE
horizontal=FALSE

valid = binsToUseInternal(object)

p1 = plotCopynumberHeatmap(object[, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names)
uniqueEventEncoding = matrix(-1, nrow=dim(object)[[1]], ncol=dim(object)[[2]])
uniqueEventIndices = object@assayData$copynumber != background & !is.na(object@assayData$copynumber)
uniqueEventEncoding[uniqueEventIndices] = object@assayData$copynumber[uniqueEventIndices]
print(prop.table(table(uniqueEventIndices)))
    
profileUnique = new("QDNAseqCopyNumbers",
                        bins = Biobase::featureData(object),
                        copynumber = uniqueEventEncoding,
                        phenodata = Biobase::phenoData(object))

  
p2 = plotCopynumberHeatmap(profileUnique[, labels(dt)], cluster_rows=dt, show_cell_names = show_cell_names)
  
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
