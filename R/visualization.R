clusterPlots <- function(seurat.obj,
                         res.range,
                         cluster.list,
                         point.size = .1,
                         label.size = 10,
                         title = "",
                         legend.title.size = 0,
                         legend.text.size = 15,
                         axis.title.x.size = 15,
                         axis.title.y.size = 15,
                         axis.text.y.left.size = 15,
                         axis.text.x.bottom.size = 15,
                         legend = TRUE,
                         cell.legend.size = 5,
                         counts.in.legend = TRUE) {
  l.coord <- cluster.list %>%
    group_by(Clusters) %>%
    summarize(UMAP1 = median(UMAP_1), UMAP2 = median(UMAP_2))
  
  use.cols <-
    hcl(h = seq(15, 375, length = length(unique(
      cluster.list$Clusters
    )) + 1),
    c = 100,
    l = 65)[seq_along(unique(cluster.list$Clusters))]
  
  p1 <- ggplot(data = cluster.list,
               mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Clusters), size = point.size) +
    scale_color_manual(values = use.cols) +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(color = 'black'),
      legend.title = element_text(size = legend.title.size),
      legend.text = element_text(size = legend.text.size),
      axis.title.x = element_text(size = axis.title.x.size),
      axis.title.y = element_text(size = axis.title.y.size),
      axis.text.y.left = element_text(size = axis.text.y.left.size),
      axis.text.x.bottom = element_text(size = axis.text.x.bottom.size)
    ) +
    guides(colour = guide_legend(override.aes = list(size = cell.legend.size)))
  
  geom_text_repel(
    data = l.coord,
    mapping = aes(x = UMAP1, y = UMAP2, label = Clusters),
    size = label.size,
    direction = "y"
  )
  ggtitle(glue(res.range, " Cluster Resolution"))
  print(p1)
  
}

#' Produces a UMAP plot from a Seurat object that can put cluster cell amounts in the legend, title, or on the plot.
#' @param annotation.name The name of the slot that contains the Seurat clusters. Defaults to "Seurat_Assignment" first, then "seurat_clusters"
#' @param point.size Size of the UMAP plot points
#' @param title Title of the plot
#' @param legend.title.size Size of the legend title
#' @param legend.text.size Size of the legend text
#' @param axis.title.x.size X-axis title text size
#' @param axis.title.y.size Y-axis title text size
#' @param axis.text.x.left.size X-axis sub-text size
#' @param axis.text.x.bottom.size  Y-axis sub-text size
#' @param counts.as.title Put the total amount of cells as the title
#' @param legend Show the legend
#' @param cell.legend.size Show the legend
#' @param counts.in.legend Put the amount of cells alongside cluster name in the legend
#' @param use.cols Use custom colors for the cell clusters
#' @return A UMAP plot
#' @examples
#' umapCellAnno(seurat.obj, annotation.name = "Seurat_Clusters")
#' @export

umapCellAnno <- function(seurat.obj,
                         annotation.name = NULL,
                         point.size = 1,
                         label.size = 10,
                         title = "",
                         legend.title.size = 0,
                         legend.text.size = 15,
                         axis.title.x.size = 15,
                         axis.title.y.size = 15,
                         axis.text.y.left.size = 15,
                         axis.text.x.bottom.size = 15,
                         counts.as.title = FALSE,
                         legend = TRUE,
                         cell.legend.size = 10,
                         counts.in.legend = TRUE,
                         use.cols = NULL) {
  if (!is.null(annotation.name)) {
    vars <- annotation.name
  } else if (!class(try(seurat.obj$Seurat_Assignment, silent = TRUE)
  ) == "try-error") {
    vars <- "Seurat_Assignment"
  } else {
    vars <- "seurat_clusters"
  }
  
  umap <-
    as_tibble(Embeddings(seurat.obj, reduction = "umap"), rownames = "Cell") %>%
    mutate(Clusters = FetchData(seurat.obj, vars = vars)[[1]])
  
  cluster.counts <- umap %>%
    group_by_at(4) %>%
    tally() %>%
    arrange(desc(n))
  
  umap$Clusters <-
    factor(umap$Clusters , levels = cluster.counts[[1]])
  umap <- umap[order(umap$Clusters),]
  
  if (is.null(use.cols)) {
    use.cols <-
      hcl(h = seq(15, 375, length = length(unique(umap[[4]])) + 1),
          c = 100,
          l = 65)[seq_along(unique(umap[[4]]))]
  }
  
  if (counts.as.title == TRUE) {
    title = paste(comma(sum(cluster.counts$n)), "Cells")
  }
  
  labels <- cluster.counts[[1]]
  if (isTRUE(counts.in.legend)) {
    labels <-
      as.character(glue("{cluster.counts[[1]]}\n ({cluster.counts[[2]]} Cells)"))
    labels <- paste0(labels, "\n")
  }
  
  l.coord <-
    umap %>% group_by(Clusters) %>% summarize(UMAP1 = median(UMAP_1),
                                              UMAP2 = median(UMAP_2))
  
  p1 <- ggplot(data = umap, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Clusters), size = point.size) +
    scale_color_manual(values = use.cols, labels = labels) +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(color = 'black'),
      legend.title = element_text(size = legend.title.size),
      legend.text = element_text(size = legend.text.size),
      axis.title.x = element_text(size = axis.title.x.size),
      axis.title.y = element_text(size = axis.title.y.size),
      axis.text.y.left = element_text(size = axis.text.y.left.size),
      axis.text.x.bottom = element_text(size = axis.text.x.bottom.size)
    ) +
    guides(colour = guide_legend(override.aes = list(size = cell.legend.size))) +
    geom_text_repel(
      data = l.coord,
      mapping = aes(x = UMAP1, y = UMAP2, label = Clusters),
      size = label.size,
      direction = "y"
    )
  return(p1)
}

#' Allows experimentation of different cluster resolutions on a Seurat object.
#' @param seurat.obj A seurat object.
#' @param cluster.res Cluster resolutions interested in
#' @param feature.plot Produce a feature plot based on given genes
#' @return UMAP and a Feature Plot
#' @examples
#' chooseClusterRes(seurat.obj, cluster.res = list(.1,.2,.5), feature.plot = c("gene1","gene2","gene3"))
#' @export

chooseClusterRes <- function(seurat.obj,
                             cluster.res,
                             ident.plot = TRUE,
                             feature.plot = NULL,
                             plot.cols = NULL) {
  umap <- try(as_tibble(Embeddings(seurat.obj, reduction = "umap"),
                        rownames = "Cells"), silent = TRUE)
  
  if (class(umap) == "try-error") {
    stop("No UMAP coordinates. Run the RunUMAP() function")
  }
  
  res.range <-
    seq(cluster.res[[1]], cluster.res[[2]], cluster.res[[3]])
  cluster.list <- map(res.range, ~ Idents(
    FindClusters(seurat.obj,
                 resolution = res.range[[.x]],
                 verbose = FALSE)
  ))
  
  cluster.list <- map(cluster.list, ~ tibble(Cells = names(.x),
                                             Clusters = .x)) %>%
    map(~ full_join(umap, .x))
  
  walk2(
    cluster.list,
    res.range,
    ~ clusterPlots(
      seurat.obj = seurat.obj,
      cluster.list = .x,
      res.range = .y
    )
  )
  
  if (!is.null(feature.plot)) {
    print(
      FeaturePlot(
        seurat.obj,
        cols = c("grey", "red"),
        features = feature.plot,
        reduction = "umap",
        ncol = plot.cols
      )
    )
  }
}

#' Produces a variable gene correlation heatmap and returns a correlation matrix.
#' @param seurat.obj A seurat object.
#' @param variable.genes The amount of variables genes to include.
#' @param downsample Amount of cells to optionally downsample.
#' @param return.table Optional argument to return the correlation matrix
#' @param fontsize The font size for the gene names on the heatmap.
#' @return A correlation heatmap and (optionally) a correlation matrix.
#' @examples
#' variableGeneMatrix(seurat.obj = seurat.obj, variable.genes = 200, return.table = FALSE, fontsize = 10)
#' @export

variableGeneMatrix <- function(seurat.obj,
                               variable.genes = 200,
                               return.table = FALSE,
                               fontsize = 10) {
  seurat.obj <-
    FindVariableFeatures(seurat.obj, nfeatures = variable.genes)
  features <- VariableFeatures(seurat.obj)
  expr.matrix <- Matrix::t(GetAssayData(seurat.obj)[features, ])
  gene.cor <- cor(as.matrix(expr.matrix))
  print(pheatmap(gene.cor, fontsize = fontsize))
  
  if (isTRUE(return.table)) {
    return(gene.cor)
  }
}

#' Produce hierarchical clustering for a sub-cluster of a downsampled Seurat object and return a dendrogram.
#'@param seurat.obj A seurat object.
#' @param cluster Cluster interested in
#' @param annotation.name Variable name given to Seurat cluster assignments
#' @param down.sample Amount of cells to sample
#' @param seed Value for the seed to set
#' @param variable.genes (Optional) Subet counts data to this many variable genes for distance matrix calculation
#' @param return.clusters (Optional) Put the height you want to cut the dendrogram at and return an object that contains cells and hierarchical clusters
#' dendoSeurat(seurat.obj, cluster = "5", annotation.name = "seurat_clusters", down.sample = 50, seed = 1, return.clusters = 4)
#' @export

dendoSeurat <- function(seurat.obj,
                        cluster,
                        annotation.name,
                        down.sample,
                        variable.genes = NULL,
                        seed = 1,
                        return.clusters = NULL) {
  set.seed(seed)
  cell.extract <-
    as_tibble(FetchData(seurat.obj, vars = annotation.name),
              rownames = "Cells") %>%
    filter(!!as.name(annotation.name) == cluster) %>%
    slice_sample(n = down.sample, replace = FALSE)
  
  if (!is.null(variable.genes)) {
    if (class(try(VariableFeatures(seurat.obj), silent = TRUE)
    ) == "try-error") {
      stop("Seurat object does not have any variable features.")
    }
    counts <-
      as_tibble(GetAssayData(seurat.obj)[VariableFeatures(seurat.obj)[seq(variable.genes)], 
                                         cell.extract$Cells],
                rownames = "Genes") %>%
      pivot_longer(cols = -(Genes), names_to = "Cells") %>%
      pivot_wider(names_from = Genes, values_from = value) %>%
      right_join(y = cell.extract, by = "Cells")
  } else {
    counts <- as_tibble(GetAssayData(seurat.obj)[, cell.extract$Cells],
                        rownames = "Genes") %>%
      pivot_longer(cols = -(Genes), names_to = "Cells") %>%
      pivot_wider(names_from = Genes, values_from = value) %>%
      right_join(y = cell.extract, by = "Cells")
  }
  
  dist <-
    hclust(dist(counts %>% select(-c(
      Cells,!!as.name(annotation.name)
    ))))
  dist$labels <- extract2(counts, annotation.name)
  
  if (!is.null(return.clusters)) {
    print(fviz_dend(tree))
    dist$labels <- extract2(counts, "Cells")
    tree.cut <- cutree(dist, h = return.clusters)
    tree.cut <- tibble(Cells = names(tree.cut), Clusters = tree.cut)
    return(counts %>% right_join(y = tree.cut, by = "Cells"))
  } else {
    return(fviz_dend(tree))
  }
}
