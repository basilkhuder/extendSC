#' @export
produceMarkers <- function(seurat.obj,
                           cells.per.ident = Inf,
                           top.gene.plot = TRUE,
                           output.name = NULL ){ 
  
  markers <- FindAllMarkers(object = seurat.obj, 
                            only.pos = TRUE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25, 
                            max.cells.per.ident = cells.per.ident)
  
  top3 <- markers %>% 
    group_by(cluster) %>% 
    slice_max(n = 3, order_by = avg_logFC) %>%
    select(gene, everything())
  print(top3)
  
  file.name <- if_else(is.null(output.name),
                       str_c(str_replace_all(Sys.Date(),"-","_"),"_markers.txt"),
                       output.name)
  write_tsv(markers, path = file.name)
  
  if(top.gene.plot == TRUE){
    cluster.averages <- AverageExpression(object = seurat.obj, return.seurat = TRUE)
    print(pheatmap(GetAssayData(cluster.averages)[unique(top3$gene),]))
  }
}
