#' Pipes together several downstream Seurat steps including variance-stabilizing transformation, PCA, clustering and nonlinear dimensionality reduction
#' @param seurat.obj A seurat object.
#' @param dims Dimensions of reduction to use as input for FindNeighbors and RunUMAP
#' @param vars_to_regress Variables to regress out in a second non-regularized linear regression for SCTransform
#' @param cluster.res Cluster resolution for FindClusters
#' @param seed.use Seed to set for RunUMAP
#' @param n.cores Number of cores to use
#' @return A processed Seurat object
#' @examples
#' processSeurat(seurat.obj, dims = 1:50, cluster.res = .2, vars_to_regress = c("orig.ident","percent.mt"))
#' @export

processSeurat <- function(seurat.obj,
                          dims = 1:50,
                          cluster.res = 0.2,
                          vars_to_regress = c("orig.ident", "percent.mt"),
                          seed.use = 24,
                          n.cores = NULL) {
  if (!is.null(n.cores)) {
    if_else(n.cores > availableCores()[[1]],
            stop("n.cores is greater than available system cores"),
            plan("multicore", workers = n.cores))
  }
  
  seurat.obj <-
    seurat.obj %>% SCTransform(vars.to.regress = vars_to_regress,
                               verbose = FALSE) %>%
    RunPCA() %>% FindNeighbors(dims = dims)
  seurat.obj <- FindClusters(seurat.obj, resolution = cluster.res)
  seurat.obj <-
    RunUMAP(
      seurat.obj,
      reduction = "pca",
      dims = dims,
      seed.use = seed.use
    )
  
  return(seurat.obj)
  
}

#' Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes
#' @param seurat.obj A seurat object.
#' @param mito.low Filter cells that have total mitochondrial genes low than this quantile
#' @param mito.high Filter cells that have total mitochondrial genes higher than this quantile
#' @param feature.cut Filter cells that have total features higher than this quantile
#' @param multiple.samples Whether the seurat.obj is made up of multiple sames (differentiated by orig.ident)
#' @examples
#' featureFiltration(seurat.obj, mito.low = .5, mito.high = .975, feature.cut = .975)
#' @export

featureFiltration <- function(seurat.obj,
                              mito.low = .05,
                              mito.high = .975,
                              feature.cut = .975,
                              multiple.samples = TRUE) {
  seurat.obj <- PercentageFeatureSet(seurat.obj,
                                     pattern = "^MT-",
                                     col.name = "percent.mt")
  
  if (isFALSE(multiple.samples)) {
    ident <- as_tibble(FetchData(
      seurat.obj,
      vars = c("percent.mt",
               "nFeature_RNA",
               "nCount_RNA")
    ), rownames = "Cells")
    
    quantile.amounts <- ident %>%
      summarize(
        mito.low = quantile(percent.mt, probs = mito.low),
        mito.high = quantile(percent.mt, probs = mito.high),
        feature.cut = quantile(nFeature_RNA, probs = feature.cut)
      )
    
    ident.filtered <-
      ident %>% filter(
        percent.mt > quantile.amounts$mito.low,
        percent.mt < quantile.amounts$mito.high,
        nFeature_RNA < quantile.amounts$feature.cut
      )
    
    if (!identical(ident.filtered$Cells, character(0))) {
      seurat.obj <- subset(seurat.obj, cells = ident.filtered$Cells)
    }
    
  } else {
    ident <- as_tibble(FetchData(
      seurat.obj,
      vars = c("orig.ident",
               "percent.mt",
               "nFeature_RNA",
               "nCount_RNA")
    ), rownames = "Cells")
    
    quantile.amounts <- ident %>%
      group_by(orig.ident) %>%
      summarize(
        mito.low = quantile(percent.mt, probs = mito.low),
        mito.high = quantile(percent.mt, probs = mito.high),
        feature.cut = quantile(nFeature_RNA, probs = feature.cut)
      )
    
    ident.filtered <- ident %>% group_split(orig.ident) %>%
      imap_dfr(
        ~ filter(
          .x,
          percent.mt > quantile.amounts$mito.low[[.y]],
          percent.mt < quantile.amounts$mito.high[[.y]],
          nFeature_RNA < quantile.amounts$feature.cut[[.y]]
        )
      )
    
    if (!identical(ident.filtered$Cells, character(0))) {
      seurat.obj <- subset(seurat.obj, cells = ident.filtered$Cells)
    }
  }
  return(seurat.obj)
}

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

#' Merge a reclusted sub-sample back into main object
#' @param seurat.obj A Seurat object.
#' @param subcluster.obj A Seurat sub-cluster object
#' @param cluster.repace Clusters to replace
#' @param annotation.name Name to store cluster annotations
#' @return Merged seurat object
#' @examples
#' mergeSubCluster(seurat.obj, subcluster.obj = seurat.sub.obj, cluster.replace = "T Cells", annotation.name = "Seurat_Assignment")
#' @export

mergeSubCluster <- function(seurat.obj,
                            subcluster.obj,
                            cluster.replace,
                            annotation.name) {
  if (is.null(annotation.name)) {
    annotation.name <- "Seurat_Assignment"
  }
  seurat.obj.ident <-
    WhichCells(object = seurat.obj, ident = cluster.replace)
  seurat.obj <-
    SetIdent(
      object = seurat.obj,
      cells = seurat.obj.ident,
      value = Idents(object = subcluster.obj)
    )
  seurat.obj@meta.data[[annotation.name]] <-
    Idents(object = seurat.obj)
  return(seurat.obj)
}

