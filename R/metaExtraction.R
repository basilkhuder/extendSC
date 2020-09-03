#' Easily extract Seurat meta-data into a tibble
#' @param seurat.obj A seurat object. 
#' @param types The types of meta-data you want to extract. Options are c("Embeddings","Clusters","Module Scores")
#' @param vars The names within the Seurat object of the meta-data you want to extract. 
#' @param merge Whether you want all data to merged. If set to TRUE, data will be merged by Cells. 
#' @return Meta-data in tibble format
#' @examples
#' extractMeta(seurat.obj, types = c("Clusters","Module Scores","Embeddings"), vars = c("Seurat_Assignment, "Module_Scores1", "umap"))
#' @export


extractMeta <- function(seurat.obj, 
                        types, 
                        vars, 
                        merge = TRUE) { 
  
  if(length(types) != length(vars)){ 
    stop("The amount of types need to be equal to the amount of variables.")
  }
  
  types <- str_to_lower(types)
  meta.list <- vector(mode = "list", length = length(types))

  if(any(types %in% "clusters")) { 
    clusters <- vars[which(types %in% "clusters")]
    clusters <- as_tibble(FetchData(seurat.obj, vars = clusters), rownames = "Cells")
    meta.list[[which(types %in% "clusters")]] <- clusters
  } 
  
  if(any(types %in% "module scores")) { 
    module.scores <- vars[which(types %in% "module scores")]
    module.scores.df <- map(module.scores, ~as_tibble(FetchData(seurat.obj, vars = .x), rownames = "Cells"))
    for(i in seq_along(module.scores.df)){ 
      meta.list[[which(types %in% "module scores")[[i]]]] <- module.scores.df[[i]]
  } 
  }
    
  if(any(types %in% "embeddings")){ 
    embeddings <- vars[which(types %in% "embeddings")]
    embeddings <- as_tibble(Embeddings(seurat.obj, reduction = embeddings), rownames = "Cells")
    meta.list[[which(types %in% "embeddings")]] <- embeddings
  }
  
  if(any(types %in% "counts")){ 
    counts <- vars[which(types %in% "counts")]
    counts <- GetAssayData(seurat.obj, assay = counts)
    meta.list[[which(types %in% "counts")]] <- counts
  }
  
  if(isTRUE(merge)){ 
    return(purrr::reduce(meta.list, full_join))
  } 
  return(meta.list)
  
} 

#' Easily extract counts from a Seurat object
#' @param seurat.obj A seurat object. 
#' @param assay The assay to pull from
#' @param genes Optional argument to extract certain genes only
#' @param tibble Whether you want counts returned as a tibble
#' @return Counts
#' @examples
#' extractCounts(seurat.obj, assay = "RNA", genes = c("TBC1D3D","LINC00470"), tibble = TRUE) 
#' @export

extractCounts <- function(seurat.obj,
                          assay = "SCT",
                          genes = NULL,
                          tibble = FALSE) { 
  
  
  counts <- GetAssayData(seurat.obj, assay = assay)
  if(!is.null(genes)){
    counts <- counts[map_dbl(genes, ~ which(rownames(counts) %in% .x)),]
    }
  
  if(isTRUE(tibble)){ 
    if(length(genes) == 1){ 
      counts <- as_tibble(counts, rownames = "Cells") %>%
        mutate(Genes = genes) %>%
        dplyr::select(Genes, everything())
      return(counts)
      }
    counts <- as_tibble(counts, rownames = "Genes") %>%
      pivot_longer(cols = -Genes)
    }
  return(counts)
 }
