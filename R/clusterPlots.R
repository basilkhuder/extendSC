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
                         counts.in.legend = TRUE){
  
  l.coord <- cluster.list %>% 
    group_by(Clusters) %>% 
    summarize(UMAP1 = median(UMAP_1), UMAP2 = median(UMAP_2))
  
  use.cols <- hcl(h = seq(15, 375, length = length(unique(cluster.list$Clusters)) + 1), 
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
        panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black'),
                     legend.title = element_text(size = legend.title.size),
                     legend.text = element_text(size = legend.text.size),
                     axis.title.x = element_text(size = axis.title.x.size),
                     axis.title.y = element_text(size = axis.title.y.size),
                     axis.text.y.left = element_text(size = axis.text.y.left.size),
                     axis.text.x.bottom = element_text(size = axis.text.x.bottom.size)) +
      guides(colour = guide_legend(override.aes = list(size=cell.legend.size))) 
    
    geom_text_repel(data = l.coord, mapping = aes(x = UMAP1, y = UMAP2,label = Clusters), size = label.size,
                    direction = "y")  
      ggtitle(glue(res.range, " Cluster Resolution"))
      print(p1)
    
} 

   
