#' Easily rename/annotate Seurat clusters or meta-data slots
#' @param seurat.obj A Seurat object.
#' @param cluster.names Names for the renamed/annotated clusters
#' @param rename.type The type of renaming you'd like. Possible values include "ident" which
#' replaces the current identity class of the seurat object, "annotation.name", which replaces a
#' specified annotation slot of the seurat.obj or "new.annotation", which creates a new slot for
#' the renamed clusters. More than one renaming type is allowed
#' @param annotation.name The slot to add renamed clusters to. Only has to be specified if
#' "annotation.name" is included within rename.type
#' @param pull.from By default, renameClusters pulls the current identity of the seurat.obj as what
#' will be renamed and then stored back (either stored back into the main identity of the object, into
#' the specified "annotation.name" or as a "new.annotation.") However, if you instead want to rename
#' clusters or annotations that are stored elsewhere, use pull.from to specify the name of the slot
#' that this information is stored under
#' @examples
#' renameClusters(seurat.obj, cluster.names = c("Cluster1","Cluster2","Cluster3"), rename.type = c("ident","annotation.name"),
#' annotation.name = "Seurat_Annotations", pull.from = NULL)
#' @export

renameClusters <- function(seurat.obj,
                           cluster.names,
                           rename.type = c("ident", "annotation.name"),
                           annotation.name = NULL,
                           new.annotation.name = NULL,
                           pull.from = NULL) {
  if (!is.null(pull.from)) {
    current.idents <-
      seurat.obj[[pull.from]] %>%
      as_tibble(rownames = "Cells") %>%
      set_colnames(c("Cells", "Idents"))
  } else {
    current.idents <-
      tibble(Cells = names(Idents(seurat.obj)),
             Idents = Idents(seurat.obj))
  }
  
  annotated.ident <- tibble(Idents = levels(current.idents$Idents),
                            Annotations = cluster.names) %>%
    full_join(current.idents) %>%
    mutate(Annotations = factor(Annotations, levels = unique(Annotations))) %>%
    select(-Idents) %>%
    relocate(Cells, .before = Annotations) %>%
    deframe()
  
  if (any(rename.type %in% "ident")) {
    Idents(seurat.obj) <- annotated.ident
  }
  
  if (any(rename.type %in% "annotation.name")) {
    seurat.obj[[annotation.name]] <- annotated.ident
    
  }
  
  if (any(rename.type %in% "new.annotation.name")) {
    seurat.obj[[new.annotation.name]] <- annotated.ident
    
  }
  return(seurat.obj)
  
}


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
