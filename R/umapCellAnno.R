# Takes a Seurat object and returns a UMAP embedding with cell counts as labels

umapCellAnno <- function(seurat.obj,
                     point.size = 1,
                     label.size = 8,
                     title = ""){  
  
  `%>%` <- magrittr::`%>%`
  extract.clusters <- data.table::setDT(FetchData(seurat.obj, vars = c("seurat_clusters")), 
                                        keep.rownames = TRUE)
  
  cluster.counts <- extract.clusters %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::tally()
  
  new.cluster.ids <- condition.extract.count$n
  names(x = new.cluster.ids) <- levels(x = seurat.obj)
  seurat.obj <- RenameIdents(object = seurat.obj, new.cluster.ids)
  print(DimPlot(object = seurat.obj, 
          reduction = 'umap', 
          pt.size = point.size, 
          label = TRUE, 
          label.size = label.size) + ggtitle(title))
  seurat.obj@active.ident <- seurat.obj$seurat_clusters
} 
