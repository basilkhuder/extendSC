#' @export
#' Produces a variable gene correlation heatmap and returns a correlation matrix. 
#' @param seurat.obj A seurat object. 
#' @param variable.genes The amount of variables genes to include. 
#' @param downsample Amount of cells to optionally downsample. 
#' @param return.table Optional argument to return the correlation matrix
#' @param fontsize The font size for the gene names on the heatmap. 
#' @return A correlation heatmap and (optionally) a correlation matrix. 
#' @examples
#' variableGeneMatrix(seurat.obj = seurat.obj, variable.genes = 200, downsample = 1000, return.table = FALSE, fontsize = 10) 

variableGeneMatrix <- function(seurat.obj,
                               variable.genes = 200,
                               downsample = NULL,
                               return.table = FALSE,
                               fontsize = 10) { 
  
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, nfeatures = variable.genes)
  features <- Seurat::VariableFeatures(seurat.obj)
  
  if(!is.null(downsample)){ 
    seurat.obj = subset(seurat.obj, cells = sample(Seurat::Cells(seurat.obj), downsample))
  }
  expr.matrix <- Matrix::t(GetAssayData(seurat.obj)[features,])
  gene.cor <- cor(as.matrix(expr.matrix))
  print(pheatmap::pheatmap(gene.cor, 
                           fontsize = fontsize))
  
  if(isTRUE(return.table)){ 
    return(gene.cor)
    } 
} 


