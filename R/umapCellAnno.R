# Takes a Seurat object and returns a UMAP embedding with cell counts as labels

umapAnno <- function(seurat.obj,
                     point.size = 1,
                     label.size = 8,
                     title = "",
                     counts.as.title = FALSE,
                     use.colors = ""){  
  
  `%>%` <- magrittr::`%>%`
  extract.clusters <- data.table::setDT(FetchData(seurat.obj, vars = c("seurat_clusters")), 
                                        keep.rownames = TRUE)
  
  cluster.counts <- extract.clusters %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::tally()
  
  names(x = cluster.counts$n) <- levels(x = seurat.obj)
  seurat.obj <- RenameIdents(object = seurat.obj, cluster.counts$n)
  
  if (counts.as.title == TRUE){ 
    title = paste(scales::comma(sum(cluster.counts$n)),"Cells")
  }
  
  if (use.colors == ""){
    use.colors <- hcl(h = seq(15, 375, length = length(unique(extract.clusters$seurat_clusters)) + 1), 
                      c = 100,
                      l = 65)[1:length(unique(extract.clusters$seurat_clusters))]
  }
  
  print(DimPlot(object = seurat.obj, 
                reduction = 'umap', 
                pt.size = point.size, 
                label = TRUE, 
                label.size = label.size,
                cols = use.colors) + ggtitle(title))
  seurat.obj@active.ident <- seurat.obj$seurat_clusters
  
} 





