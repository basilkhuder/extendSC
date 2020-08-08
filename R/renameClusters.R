#' Easily rename/annotate Seurat clusters
#' @param seurat.obj A Seurat object. 
#' @param cluster.names Names for the renamed/annotated clusters
#' @param annotation.name Name to store new clusters in Seurat object
#' @return A Seurat object 
#' @examples
#' renameClusters(seurat.obj, cluster.names = c("Cluster1","Cluster2","Cluster3"), annotation.name = "Seurat_Annotations")
#' @export
#' @export
renameClusters <- function(seurat.obj, 
                           cluster.names,
                           annotation.name = NULL) { 
  
  if(is.null(annotation.name)){ 
    annotation.name <- "Seurat_Assignment"
    }
  
  names(cluster.names) <- levels(x = seurat.obj)
  seurat.obj <- RenameIdents(seurat.obj, cluster.names)
  seurat.obj@meta.data[[annotation.name]] <- Idents(seurat.obj)
  return(seurat.obj)
  
}




