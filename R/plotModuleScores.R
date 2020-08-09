#' @export

plotModuleScores <- function(seurat.obj,
                             features,
                             name = NULL){ 
  
  if(!is.null(name)){ 
    name <- name
  } else { 
    name <- "Module_scores"
  }
  m.obj <- AddModuleScore(object = seurat.obj,
                                  features = features,
                                  name = "Module_Scores")
  print(RidgePlot(m.obj, features = features))
  print(FeaturePlot(m.obj, features = features)) 
 } 
                             
