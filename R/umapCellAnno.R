# Takes a Seurat object and returns a UMAP embedding with cell counts as labels

umapAnno <- function(seurat.obj,
                     point.size = 1,
                     label.size = 8,
                     title = "",
                     counts.as.title = FALSE,
                     legend.pos = "right",
                     use.colors = ""){  
  
  `%>%` <- magrittr::`%>%`
  extract.clusters <- data.table::setDT(Seurat::FetchData(seurat.obj, vars = c("seurat_clusters")), 
                                        keep.rownames = TRUE)
  
  cluster.counts <- extract.clusters %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::tally() %>%
    dplyr::pull(n) %>%
    data.table::setattr("names",levels(seurat.obj))
    
  if (counts.as.title == TRUE){ 
    title = paste(scales::comma(sum(cluster.counts$n)),"Cells")
  }
  
  
  if (use.colors == ""){
    use.colors <- hcl(h = seq(15, 375, length = length(unique(extract.clusters$seurat_clusters)) + 1), 
                      c = 100,
                      l = 65)[1:length(unique(extract.clusters$seurat_clusters))]
  }
    
  return(Seurat::DimPlot(object = RenameIdents(object = seurat.obj, cluster.counts), 
                reduction = 'umap', 
                pt.size = point.size, 
                label = TRUE, 
                label.size = label.size,
                cols = use.colors) + ggtitle(title) + theme(legend.position = legend.pos))    
} 
