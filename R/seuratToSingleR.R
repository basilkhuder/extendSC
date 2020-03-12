#Produce SingleR annotations on a custom reference made from a Seurat object 

seuratToSingleR <- function(object = seurat.obj,
                            reference = reference.obj,
                            method = method){ 
  
  seurat_obj_avg <- Seurat::AverageExpression(reference_obj, assay = "RNA")
  labels <- colnames(reference_obj_avg)
  SingleR::SingleR(test = as.SingleCellExperiment(seurat_obj), 
                   ref = seurat_obj_ag, labels = labels, 
                   method = method)
 


                            
                            
                            
                            
} 
