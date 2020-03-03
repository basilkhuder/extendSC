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
  
  new.cluster.ids <- cluster.counts$n
  names(x = new.cluster.ids) <- levels(x = seurat.obj)
  seurat.obj <- RenameIdents(object = seurat.obj, new.cluster.ids)
  
  if (counts.as.title == TRUE){ 
    title <- paste(scales::comma(sum(cluster.counts$n)),"cells")
    }
  
  print(DimPlot(object = seurat.obj, 
                reduction = 'umap', 
                pt.size = point.size, 
                label = TRUE, 
                label.size = label.size,
                cols = use.colors) + ggtitle(title))
  seurat.obj@active.ident <- seurat.obj$seurat_clusters
} 
