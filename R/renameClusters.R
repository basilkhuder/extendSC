renameClusters <- function(seurat.obj, cluster_names){ 
  seurat.obj.cluster.ids <- cluster_names
  names(seurat.obj.cluster.ids) <- levels(seurat.obj)
  seurat.obj <- RenameIdents(seurat.obj, seurat.obj.cluster.ids)
  seurat.obj@meta.data$Seurat_Assignment <- Idents(seurat.obj)
  return(seurat.obj)
  }
