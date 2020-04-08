mergeSubCluster <- function(seurat.obj, subcluster.obj, cluster_replace){ 
  seurat.obj.ident <- Seurat::WhichCells(object = seurat.obj, ident = cluster_replace)
  seurat.obj <- Seurat::SetIdent(object = seurat.obj, cells = seurat.obj.ident,
                         value = Seurat::Idents(object = subcluster.obj))
  seurat.obj@meta.data$Seurat_Assignment <- Idents(object = seurat.obj)
  return(seurat.obj)
}
  
