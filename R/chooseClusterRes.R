chooseClusterRes <- function(seurat.obj.list){ 

  umap.list <- lapply(seq_along(seurat.obj.list), function(x) umapCellAnno(seurat.obj.list[[x]]) 



} 
