umapCellAnno <- function(seurat.obj,
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
                         counts.in.legend = TRUE,
                         use.cols = NULL){
  
  `%>%` <- magrittr::`%>%`
  cnames <- try(seurat.obj$Seurat_Assignment, silent = TRUE)
  if(!class(cnames) == "try-error"){ 
    vars <- "Seurat_Assignment"
  } else {
      vars <- "seurat_clusters"
    }

  umap <- as.data.frame(Embeddings(seurat.obj, reduction = "umap")) %>%
    dplyr::mutate(Clusters = Seurat::FetchData(seurat.obj, vars = vars)[[1]]) 
  extract.clusters <- data.table::setDT(Seurat::FetchData(seurat.obj, vars = vars),keep.rownames = TRUE)
  cluster.counts <- extract.clusters %>%
    dplyr::group_by_at(2) %>%
    dplyr::tally() %>%
    dplyr::arrange(desc(n))
  umap$Clusters <- factor(umap$Clusters , 
                          levels = cluster.counts$Seurat_Assignment)
  umap <- umap[order(umap$Clusters), ]
  
  
  if (is.null(use.cols)){
    use.cols <- hcl(h = seq(15, 375, length = length(unique(extract.clusters[[2]])) + 1), 
                    c = 100,
                    l = 65)[1:length(unique(extract.clusters[[2]]))]
  }
  
  if (counts.as.title == TRUE){ 
    title = paste(scales::comma(sum(cluster.counts$n)),"Cells")
  }
  
  p1 <- ggplot2::ggplot(data = umap, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    ggplot2::geom_point(aes(color = Clusters), size = point.size) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
    ggplot2::theme(axis.line = element_line(color = 'black'),
                   legend.title = element_text(size = legend.title.size),
                   legend.text = element_text(size = legend.text.size),
                   axis.title.x = element_text(size = axis.title.x.size),
                   axis.title.y = element_text(size = axis.title.y.size),
                   axis.text.y.left = element_text(size = axis.text.y.left.size),
                   axis.text.x.bottom = element_text(size = axis.text.x.bottom.size)) +
    ggplot2::guides(colour = guide_legend(override.aes = list(size=10))) 
  
  p1 <- Seurat::LabelClusters(p1, id = "Clusters", size = label.size, repel = TRUE) + 
    ggplot2::ggtitle(title)
  
  if (isTRUE(counts.in.legend)){ 
    labels <- as.character(glue::glue("{cluster.counts[[1]]}\n ({cluster.counts[[2]]} Cells)"))
    labels <- paste0(labels,"\n")
    p1 <- p1 + ggplot2::scale_color_manual(values = use.cols, labels = labels)
  }
  return(p1)
} 
