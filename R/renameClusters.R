renameClusters <- function(seurat.obj, cluster_names){ 
  names(cluster_names) <- levels(x = seurat.obj)
  seurat.obj <- RenameIdents(seurat.obj, seurat.obj.cluster.ids)
  seurat.obj@meta.data$Seurat_Assignment <- Idents(seurat.obj)
  return(seurat.obj)
}

