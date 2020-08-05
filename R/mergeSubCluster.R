#' Merge a reclusted sub-sample back into main object
#' @param seurat.obj A Seurat object. 
#' @param subcluster.obj A Seurat sub-cluster object
#' @param cluster.repace Clusters to replace
#' @param annotation.name Name to store cluster annotations
#' @return UMAP and a Feature Plot
#' @examples
#' chooseClusterRes(seurat.obj, cluster.res = list(.1,.2,.5), feature.plot = c("gene1","gene2","gene3"))
#' @export
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
  
