#' @export
produceMarkers <- function(seurat.obj,
                           cells.per.ident = Inf,
                           top.gene.plot = TRUE,
                           file.name = NULL ){ 
  
  markers <- FindAllMarkers(object = seurat.obj, 
                            only.pos = TRUE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25, 
                            max.cells.per.ident = cells.per.ident)
  
  print(markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC))
  
  if(!is.null(file.name)){ 
    file.name <- file.name
  } else { 
    file.name <- paste0(gsub("-","",Sys.Date()),"markers.txt")
  } 
  write.table(markers %>% group_by(cluster), 
              file = file.name,
              row.names=FALSE, sep="\t")
  if(top.gene.plot == TRUE){
    cluster.averages <- AverageExpression(object = seurat.obj, return.seurat = TRUE)
    top3 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
    print(pheatmap(GetAssayData(cluster.averages)[unique(top3$gene),]))
    }
  }
