processSeurat <- function(seurat.obj,
                          dims = 1:50,
                          cluster.res = 0.2,
                          vars_to_regress = c("orig.ident","percent.mt"),
                          seed.use = 24,
                          n.cores = NULL) {
  
  if(!is.null(n.cores)){
    ifelse(n.cores > future::availableCores()[[1]], 
           stop("n.cores is greater than available system cores"),
           future::plan("multicore", workers = n.cores))
  }
  
  seurat.obj <- seurat.obj %>%
                SCTransform(vars.to.regress = vars_to_regress,
                            verbose = FALSE) %>%
                RunPCA() %>%
                FindNeighbors(dims = dims) %>%
                FindClusters(resolution = cluster.res) %>%
                RunUMAP(reduction = "pca", dims = dims, seed.use = seed.use)

  return(seurat.obj)
}
                  
