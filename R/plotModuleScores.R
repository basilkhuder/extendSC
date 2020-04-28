#' @export

plotModuleScores <- function(seurat.obj,
                             features,
                             name = NULL){ 
  
  if(!is.null(name)){ 
    name <- name
  } else { 
    name <- "Module_scores"
  }
  m.obj <- Seurat::AddModuleScore(object = seurat.obj,
                                  features = features,
                                  name = "Module_Scores")
  print(Seurat::RidgePlot(m.obj, features = features))
  print(Seurat::FeaturePlot(m.obj, features = features)) 
 } 
                             
