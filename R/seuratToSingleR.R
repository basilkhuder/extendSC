#Produce SingleR annotations on a custom reference made from a Seurat object 

seuratToSingleR <- function(object = seurat.obj,
                            reference = reference.obj,
                            method = method,
                            heatmap = TRUE){ 
  
  reference_avg <- as.matrix(Seurat::AverageExpression(reference_obj, assay = "RNA")$RNA)
  labels <- colnames(reference_avg)
  SR.output <- SingleR::SingleR(test = as.SingleCellExperiment(seurat_obj), 
                   ref = reference_avg, labels = labels, 
                   method = method)
  return(pheatmap::pheatmap(SR.output))
  } 
