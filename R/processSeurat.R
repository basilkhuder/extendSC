processSeurat <- function(seurat.obj,
                          dims = 1:50,
                          cluster.res = 0.2,
                          vars_to_regress = c("orig.ident","percent.mt"),
                          seed.use = 24,
                          n.cores = NULL) {
  `%>%` <- magrittr::`%>%`
  if(!is.null(n.cores)){
    ifelse(n.cores > future::availableCores()[[1]], 
           stop("n.cores is greater than available system cores"),
           future::plan("multicore", workers = n.cores))
  }
  seurat.obj <- seurat.obj %>%
    SCTransform(vars.to.regress = vars_to_regress,
                verbose = FALSE) %>%
    RunPCA() %>%
    FindNeighbors(dims = dims) 
  if(class(cluster.res) == "list"){ 
    res.range <- seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
    seurat.list <- lapply(seq_along(seurat.list), function(x)  
      FindClusters(seurat.obj, resolution = res.range[[x]]))
    seurat.list <- lapply(seq_along(seurat.list), function(x)  
      RunUMAP(seurat.list[[x]], reduction = "pca", dims = dims, seed.use = seed.use)) 
  } else { 
    seurat.obj <- FindClusters(seurat.obj, resolution = cluster.res)
    seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", dims = dims, seed.use = seed.use)
  }
  return(seurat.obj)
}
