#' Allows experimentation of different cluster resolutions on a Seurat object. 
#' @param seurat.obj A seurat object. 
#' @param cluster.res Cluster resolutions interested in 
#' @param ident.plot Produce an identity plot
#' @param feature.plot Produce a feature plot based on given genes
#' @param plot.cols Choose the number of columns for the figures.
#' @return UMAP and identity/feature plots
#' @examples
#' chooseClusterRes(seurat.obj, cluster.res = list(.1,.2,.5), feature.plot = c("gene1","gene2","gene3"))
#' @export

chooseClusterRes <- function(seurat.obj, 
                             cluster.res,
                             ident.plot = TRUE,
                             feature.plot = NULL,
                             plot.cols = NULL){
  
  umap <- try(as_tibble(Embeddings(seurat.obj,reduction = "umap"), 
                        rownames = "Cells"), silent = TRUE)
  
  if(class(umap) == "try-error"){ 
    stop("No UMAP coordinates. Run the RunUMAP() function")
  } 
  
  res.range <- seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
  cluster.list <- map(res.range, ~ Idents(FindClusters(seurat.obj, 
                                                        resolution = res.range[[.x]], 
                                                        verbose = FALSE)))
  
  cluster.list <- map(cluster.list, ~tibble(Cells = names(.x), 
                                            Cluster = .x)) %>%
    map( ~full_join(umap, .x))
  
  clusterPlots(seurat.obj = seurat.obj,
               cluster.list = cluster.list,
               cluster.res = res.range,
               feature.plot = feature.plot,
               plot.cols = plot.cols,
               ident.plot = ident.plot)
} 
