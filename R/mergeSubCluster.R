#' Merge a reclusted sub-sample back into main object
#' @param seurat.obj A Seurat object. 
#' @param subcluster.obj A Seurat sub-cluster object
#' @param cluster.repace Clusters to replace
#' @param annotation.name Name to store cluster annotations
#' @return Merged seurat object
#' @examples
#' mergedSubCluster(seurat.obj, subcluster.obj = seurat.sub.obj, cluster.replace = "T Cells", annotation.name = "Seurat_Assignment") 
#' @export
mergeSubCluster <- function(seurat.obj, 
                            subcluster.obj, 
                            cluster.replace,
                            annotation.name){ 
  
  if(is.null(annotation.name)) {
    annotation.name <- "Seurat_Assignment"
  }
  seurat.obj.ident <- WhichCells(object = seurat.obj, ident = cluster.replace)
  seurat.obj <- SetIdent(object = seurat.obj, cells = seurat.obj.ident,
                         value = Idents(object = subcluster.obj))
  seurat.obj@meta.data[[annotation.name]] <- Idents(object = seurat.obj)
  return(seurat.obj)
}
  
