processSeurat <- function(seurat.obj,
                          dims = 1:50,
                          cluster.res = 0.2,
                          seed.use = 24,
                          n.cores = NULL) {
  `%>%` <- magrittr::`%>%`
  if(!is.null(n.cores)){
    ifelse(n.cores > future::availableCores()[[1]], 
           stop("n.cores is greater than available system cores"),
           future::plan("multicore", workers = n.cores))
  }
  seurat.obj <- seurat.obj %>%
                Seurat::SCTransform(vars.to.regress = c("percent.mt","orig.ident"), verbose = FALSE) %>%
                Seurat::RunPCA() %>%
                Seurat::FindNeighbors(dims = dims) %>%
                Seurat::RunUMAP(reduction = "pca", dims = dims, seed.use = seed.use)

  return(seurat.obj)
}
                                     
