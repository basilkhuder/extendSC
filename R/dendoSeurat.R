#' Produce hierarchical clustering for a sub-cluster of a downsampled Seurat object and return a dendrogram. 
#' @param seurat.obj A seurat object. 
#' @param cluster Cluster interested in
#' @param annotation.name Variable name given to Seurat cluster assignments 
#' @param down.sample Amount of cells to sample 
#' @param seed Value for the seed to set
#' @param variable.genes (Optional) Subet counts data to this many variable genes for distance matrix calculation
#' @param return.clusters (Optional) Put the height you want to cut the dendrogram at and return an object that contains cells and hierarchical clusters
#' dendoSeurat(seurat.obj, cluster = "5", annotation.name = "seurat_clusters", down.sample = 50, seed = 1, return.clusters = 4) 
#' @export

dendoSeurat <- function(seurat.obj,
                        cluster,
                        annotation.name,
                        down.sample,
                        variable.genes = NULL,
                        seed = 1,
                        return.clusters = NULL) { 
  
  set.seed(seed)
  cell.extract <- as_tibble(FetchData(seurat.obj, vars = annotation.name),
                            rownames = "Cells") %>%
    filter(!!as.name(annotation.name) == cluster) %>%
    slice_sample(n = down.sample, replace = FALSE)
  
  if(!is.null(variable.genes)) { 
    if(class(try(VariableFeatures(seurat.obj), silent = TRUE)) == "try-error"){
      stop("Seurat object does not have any variable features.")
    }
    counts <- as_tibble(GetAssayData(seurat.obj)[VariableFeatures(seurat.obj)[seq(variable.genes)],cell.extract$Cells], 
                        rownames = "Genes") %>% 
      pivot_longer(cols = -(Genes), names_to = "Cells") %>%
      pivot_wider(names_from = Genes, values_from = value) %>% 
      right_join(y = cell.extract, by = "Cells")
  } else { 
    counts <- as_tibble(GetAssayData(seurat.obj)[,cell.extract$Cells], 
                        rownames = "Genes") %>% 
      pivot_longer(cols = -(Genes), names_to = "Cells") %>%
      pivot_wider(names_from = Genes, values_from = value) %>% 
      right_join(y = cell.extract, by = "Cells")
  }
  
  dist <- dist(counts %>% select(-c(Cells, !!as.name(annotation.name))))
  tree <- hcut(dist, hang = .1)
  tree$labels<- extract2(counts, annotation.name)
  
  if(!is.null(return.clusters)){ 
    print(fviz_dend(tree))
    tree$labels<- extract2(counts, "Cells")
    tree.cut <- cutree(tree, h = return.clusters)
    tree.cut <- tibble(Cells = names(tree.cut), Clusters = tree.cut)
    return(counts %>% right_join(y = tree.cut, by = "Cells"))
  } else { 
    return(fviz_dend(tree))
    }
} 
