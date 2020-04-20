clusterPlots <- function(cluster.list,
                         cluster.res,
                         point.size = 1,
                         label.size = 10,
                         title = "",
                         legend.title.size = 0,
                         legend.text.size = 15,
                         axis.title.x.size = 15,
                         axis.title.y.size = 15, 
                         axis.text.y.left.size = 15,
                         axis.text.x.bottom.size = 15,
                         legend = TRUE,
                         cell.legend.size = 10, 
                         counts.in.legend = TRUE){
    
    `%>%` <- magrittr::`%>%`
    for(i in 1:length(cluster.list)){ 
        title <- glue::glue(cluster.res[[i]], " Cluster Resolution")
        use.cols <- hcl(h = seq(15, 375, length = length(levels(cluster.list[[i]]$V3)) + 1), 
                        c = 100,
                        l = 65)[1:length(levels(cluster.list[[i]]$V3))]
        p1 <- ggplot2::ggplot(data = cluster.list[[i]], 
                              mapping = aes(x = UMAP_1, y = UMAP_2)) +
            ggplot2::geom_point(aes(color = V3), size = point.size) +
            scale_color_manual(values = use.cols) + 
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
            ggplot2::guides(colour = guide_legend(override.aes = list(size=cell.legend.size))) 
        
        print(Seurat::LabelClusters(p1, id = "V3", size = label.size, repel = TRUE) + 
                   ggplot2::ggtitle(title))
    } 
} 
