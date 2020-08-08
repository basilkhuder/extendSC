#' Easily rename/annotate Seurat clusters
#' @param seurat.obj A Seurat object. 
#' @param cluster.names Names for the renamed/annotated clusters
#' @param annotation.name Name to store new clusters in Seurat object
#' @return A Seurat object 
#' @examples
#' crToSeurat(directory = "Files", min_cells = 3, min_features = 3)
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




