#Produce SingleR annotations on a custom reference made from a Seurat object 

seuratToSingleR <- function(ref = seurat.obj,
                            object = sce.obj
                            method = method,
                            heatmap = TRUE){ 
  
  reference_avg <- as.matrix(Seurat::AverageExpression(seurat.obj, assay = "RNA")$RNA)
  labels <- colnames(seurat.obj)
  SR.output <- SingleR::SingleR(test = sce.obj, 
                   ref = reference_avg, labels = labels, 
                   method = method)
  return(pheatmap::pheatmap(SR.output))
  } 
