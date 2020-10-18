annotateBySingleR <- function(seurat.obj,
                              database = "hpca",
                              cluster.slot = "seurat_clusters",
                              sce.assay = NULL,
                              annotate.type = c("new.annotation.name", "ident"),
                              annotation.name = NULL) {
  sce.obj <- as.SingleCellExperiment(seurat.obj, assay = sce.assay)
  sce.database <- HumanPrimaryCellAtlasData()
  sce.obj <- SingleR(
    test = sce.obj,
    ref = sce.database,
    labels = sce.database$label.main,
    method = "cluster",
    clusters = sce.obj[[cluster.slot]]
  )
  singler.annot <- sce.obj$pruned.labels
  
  
  if (any(annotate.type %in% "new.annotation.name")) {
    seurat.obj <- renameClusters(
      seurat.obj,
      cluster.names = singler.annot,
      rename.type = "new.annotation.name",
      new.annotation.name = annotation.name,
      pull.from = cluster.slot
    )
    
  }
  
  if (any(annotate.type %in% "ident")) {
    seurat.obj <- renameClusters(
      seurat.obj,
      cluster.names = singler.annot,
      rename.type = "ident",
      pull.from = cluster.slot
    )
  }
  
  return(seurat.obj)
  
}
