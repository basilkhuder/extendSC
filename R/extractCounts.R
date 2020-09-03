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
