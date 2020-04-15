chooseClusterRes <- function(seurat.obj, res.range){
  umap <- try(Embeddings(seurat.obj, reduction = "umap"))
  if(class(umap) == "try-error"){ 
    stop("No UMAP coordinates. Run the RunUMAP() function")
    } 
    res.range <- seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
    cluster.list <- lapply(seq_along(res.range), function(x)
        Idents(FindClusters(seurat.obj, 
                            resolution = res.range[[x]], 
                            verbose = FALSE)))
    return(clusterPlots(umap, cluster.list))
} 

