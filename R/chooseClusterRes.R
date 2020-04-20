chooseClusterRes <- function(seurat.obj, cluster.res){
    `%>%` <- magrittr::`%>%`
    umap <- try(Embeddings(seurat.obj, reduction = "umap"))
    if(class(umap) == "try-error"){ 
        stop("No UMAP coordinates. Run the RunUMAP() function")
    } 
    res.range <- seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
    cluster.list <- lapply(seq_along(res.range), function(x)
        Seurat::Idents(Seurat::FindClusters(seurat.obj, 
                            resolution = res.range[[x]], 
                            verbose = FALSE)))
    cluster.list <- lapply(seq_along(cluster.list), function(x)
      as.data.frame(cbind(umap, cluster.list[[x]])))
    
    cluster.counts <- vector(mode = "list", length = length(cluster.list))
    for(i in 1:length(cluster.counts)){
      cluster.counts[[i]] <- cluster.list[[i]] %>%
      dplyr::group_by_at(3) %>%
      dplyr::tally() %>%
      dplyr::arrange(desc(n))
      cluster.list[[i]]$V3 <- factor(cluster.list[[i]]$V3 , levels = cluster.counts[[i]][[1]])
      cluster.list[[i]] <- cluster.list[[i]][order(cluster.list[[i]]$V3), ]
    }
    
    clusterPlots(cluster.list = cluster.list,
                 cluster.res = res.range)
} 
