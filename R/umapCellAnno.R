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
                         use.cols = NULL){
  

  if(!is.null(annotation.name)) { 
    vars <- annotation.name
  } else if (!class(try(seurat.obj$Seurat_Assignment, silent = TRUE)) == "try-error") { 
    vars <- "Seurat_Assignment"
  } else {
    vars <- "seurat_clusters"
  }
  
  umap <- as_tibble(Embeddings(seurat.obj, reduction = "umap"), rownames = "Cell") %>%
    mutate(Clusters = FetchData(seurat.obj, vars = vars)[[1]])
  
  cluster.counts <- umap %>% 
    group_by_at(4) %>% 
    tally() %>% 
    arrange(desc(n))
  
  umap$Clusters <- factor(umap$Clusters , levels = cluster.counts[[1]])
  umap <- umap[order(umap$Clusters), ]
  
  if (is.null(use.cols)){
    use.cols <- hcl(h = seq(15, 375, length = length(unique(umap[[4]])) + 1), 
                    c = 100,
                    l = 65)[seq_along(unique(umap[[4]]))]
  }
  
  if (counts.as.title == TRUE){ 
    title = paste(comma(sum(cluster.counts$n)),"Cells")
  }
  
  labels <- cluster.counts[[1]]
  if (isTRUE(counts.in.legend)){ 
    labels <- as.character(glue("{cluster.counts[[1]]}\n ({cluster.counts[[2]]} Cells)"))
    labels <- paste0(labels,"\n")
  } 
  
  l.cord <- umap %>% group_by(Clusters) %>% summarize(UMAP1 = median(UMAP_1), 
                                                      UMAP2 = median(UMAP_2))
  
  p1 <- ggplot(data = umap, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Clusters), size = point.size) +
    scale_color_manual(values = use.cols, labels = labels) + 
    theme_bw() + 
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black'),
                   legend.title = element_text(size = legend.title.size),
                   legend.text = element_text(size = legend.text.size),
                   axis.title.x = element_text(size = axis.title.x.size),
                   axis.title.y = element_text(size = axis.title.y.size),
                   axis.text.y.left = element_text(size = axis.text.y.left.size),
                   axis.text.x.bottom = element_text(size = axis.text.x.bottom.size)) +
    guides(colour = guide_legend(override.aes = list(size=cell.legend.size))) +
    geom_text_repel(data = l.cord, mapping = aes(x = UMAP1,y=UMAP2,label = Clusters), size = label.size,
                    direction = "y")  
  return(p1)
} 

