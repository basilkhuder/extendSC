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
    Seurat::SCTransform(vars.to.regress = vars_to_regress,
                verbose = FALSE) %>%
    Seurat::RunPCA() %>%
    Seurat::FindNeighbors(dims = dims) 
  if(class(cluster.res) == "list"){ 
    res.range <- seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
    seurat.list <- lapply(seq_along(res.range), function(x)  
      Seurat::FindClusters(seurat.obj, resolution = res.range[[x]]))
    seurat.list <- lapply(seq_along(seurat.list), function(x)  
      Seurat::RunUMAP(seurat.list[[x]], reduction = "pca", dims = dims, seed.use = seed.use)) 
  } else { 
    seurat.obj <- Seurat::FindClusters(seurat.obj, resolution = cluster.res)
    seurat.obj <- Seurat::RunUMAP(seurat.obj, reduction = "pca", dims = dims, seed.use = seed.use)
  }
  return(seurat.obj)
}
