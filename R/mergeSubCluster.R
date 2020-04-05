mergeSubCluster <- function(seurat.obj, subcluster.obj, cluster_replace){ 
  seurat.obj.ident <- WhichCells(object = seurat.obj, ident = cluster_replace)
  seurat.obj <- SetIdent(object = seurat.obj, cells = seurat.obj.ident,
                                              value = Idents(object = subcluster.obj))
  seurat.obj@meta.data$Seurat_Assignment <- Idents(object = seurat.obj)
  return(seurat.obj)
}
  
