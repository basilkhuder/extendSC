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
    title = paste(scales::comma(sum(cluster.counts$n)),"Cells")
  }
  
  if (use.colors == ""){
    hues <- seq(15, 375, length = length(unique(extract.clusters$seurat_clusters)) + 1)
    use.colors <- hcl(h = hues, l = 65, c = 100)[1:length(unique(extract.clusters$seurat_clusters))]
  }
    
  print(DimPlot(object = seurat.obj, 
                reduction = 'umap', 
                pt.size = point.size, 
                label = TRUE, 
                label.size = label.size,
                cols = use.colors) + ggtitle(title))
  seurat.obj@active.ident <- seurat.obj$seurat_clusters
} 
