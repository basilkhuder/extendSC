chooseClusterRes <- function(seurat.obj, 
                             cluster.res,
                             feature.plot = NULL,
                             plot.cols = NULL){
  
  `%>%` <- magrittr::`%>%`
  umap <- try(as.data.frame(Seurat::Embeddings(seurat.obj, 
                                               reduction = "umap")))
  if(class(umap) == "try-error"){ 
    stop("No UMAP coordinates. Run the RunUMAP() function")
  } 
  umap <- data.table::setDT(umap, keep.rownames = TRUE)
  res.range <- seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
  
  cluster.list <- lapply(seq_along(res.range), function(x)
    tibble::enframe(Seurat::Idents(Seurat::FindClusters(seurat.obj, 
                                                        resolution = res.range[[x]], 
                                                        verbose = FALSE)),
                    name = "rn", value = "Clusters"))
  cluster.list <- lapply(seq_along(res.range), function(x)
    dplyr::full_join(umap, cluster.list[[x]]))
  
  cluster.counts <- vector(mode = "list", length = length(cluster.list))
  for(i in 1:length(cluster.counts)){
    cluster.counts[[i]] <- cluster.list[[i]] %>%
      dplyr::group_by_at(4) %>%
      dplyr::tally() %>%
      dplyr::arrange(desc(n))
    cluster.list[[i]]$Clusters <- factor(cluster.list[[i]]$Clusters , 
                                         levels = cluster.counts[[i]][[1]])
    cluster.list[[i]] <- cluster.list[[i]][order(cluster.list[[i]]$Clusters), ]
  }
  
  clusterPlots(seurat.obj = seurat.obj,
               cluster.list = cluster.list,
               cluster.res = res.range,
               feature.plot = feature.plot,
               plot.cols = plot.cols)
} 
