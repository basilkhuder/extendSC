variableGeneMatrix <- function(seurat.obj,
                               variable.genes = 200,
                               downsample = NULL) { 
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, nfeatures = variable.genes)
  features <- Seurat::VariableFeatures(seurat.obj)
  
  if(!is.null(downsample)){ 
    seurat.obj = subset(seurat.obj, cells = sample(Seurat::Cells(seurat.obj), downsample))
  }
  expr.matrix <- GetAssayData(seurat.obj)[features,]
  gene.cor <- cor(t(gene.expr))
  print(pheatmap::pheatmap(gene.cor))
} 


