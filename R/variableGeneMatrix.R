#' @export
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


