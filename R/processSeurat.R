#' Pipes together several downstream Seurat steps including variance-stabilizing transformation, PCA, clustering and nonlinear dimensionality reduction 
#' @param seurat.obj A seurat object. 
#' @param dims Dimensions of reduction to use as input for FindNeighbors and RunUMAP
#' @param vars_to_regress Variables to regress out in a second non-regularized linear regression for SCTransform
#' @param cluster.res Cluster resolution for FindClusters
#' @param seed.use Seed to set for RunUMAP
#' @param n.cores Number of cores to use
#' @return A processed Seurat object 
#' @examples
#' processSeurat(seurat.obj, dims = 1:50, cluster.res = .2, vars_to_regress = c("orig.ident","percent.mt"))
#' @export

processSeurat <- function(seurat.obj,
                          dims = 1:50,
                          cluster.res = 0.2,
                          vars_to_regress = c("orig.ident","percent.mt"),
                          seed.use = 24,
                          n.cores = NULL) {
  if(!is.null(n.cores)){
    if_else(n.cores > availableCores()[[1]], 
            stop("n.cores is greater than available system cores"),
            plan("multicore", workers = n.cores))
  }
  
  seurat.obj <- seurat.obj %>% SCTransform(vars.to.regress = vars_to_regress, 
                                           verbose = FALSE) %>%
    RunPCA() %>% FindNeighbors(dims = dims) 
    seurat.obj <- FindClusters(seurat.obj, resolution = cluster.res)
    seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", dims = dims, seed.use = seed.use)
  
  return(seurat.obj)
}
