#' Produce hierachial clustering for a sub-cluster of a downsampled Seurat object. 
#' @param seurat.obj A seurat object. 
#' @param clusters Cluster to cluster
#' @param annotation.name Variable name given to Seurat cluster assignments 
#' @param down.sample Amount of cells to sample 
#' @param seed Value for the seed to set

#' @param variable.genes (Optional) Subet counts data to this many variable genes for distance matrix calculation
#' hierachialSeurat(seurat.obj, clusters = "Cluster1", annotation.name = "Seurat_Assignment", down.sample = 50) 
#' @export

hierachialSeurat <- function(seurat.obj,
                             clusters,
                             annotation.name,
                             down.sample,
                             variable.genes = NULL,
                             seed = 1) { 
  
  set.seed(seed)
  cell.extract <- as_tibble(FetchData(seurat.obj, vars = annotation.name),
                            rownames = "Cells") %>%
    filter(!!as.name(annotation.name) == clusters) %>%
    slice_sample(n = down.sample, replace = FALSE)
  
  if(!is.null(variable.genes)) { 
    
    if(class(try(VariableFeatures(seurat.obj), silent = TRUE)) == "try-error"){
      stop("Seurat object does not have any variable features.")
    }
    
    counts <- as_tibble(GetAssayData(seurat.obj)[seq(variable.genes),cell.extract$Cells], 
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
  return(fviz_dend(tree))
  
} 
