variableGeneMatrix <- function(seurat.obj,
                               variable.genes = 200,
                               downsample = NULL,
                               return.table = FALSE) { 
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, nfeatures = variable.genes)
  features <- Seurat::VariableFeatures(seurat.obj)
  
  if(!is.null(downsample)){ 
    seurat.obj = subset(seurat.obj, cells = sample(Seurat::Cells(seurat.obj), downsample))
  }
  expr.matrix <- GetAssayData(seurat.obj)[features,]
  gene.cor <- cor(t(expr.matrix))
  print(pheatmap::pheatmap(gene.cor))
  
  if(isTRUE(return.table)){ 
    return(gene.cor)
    } 
} 


