clusterPlots <- function(seurat.obj,
                         cluster.list,
                         cluster.res,
                         feature.plot,
                         ident.plot,
                         plot.cols,
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
    
    `%>%` <- magrittr::`%>%`
    for(i in 1:length(cluster.list)){ 
        title <- glue::glue(cluster.res[[i]], " Cluster Resolution")
        use.cols <- hcl(h = seq(15, 375, length = length(levels(cluster.list[[i]]$Clusters)) + 1), 
                        c = 100,
                        l = 65)[1:length(levels(cluster.list[[i]]$Clusters))]
        p1 <- ggplot2::ggplot(data = cluster.list[[i]], 
                              mapping = aes(x = UMAP_1, y = UMAP_2)) +
            ggplot2::geom_point(aes(color = Clusters), size = point.size) +
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
        
        print(Seurat::LabelClusters(p1, id = "Clusters", size = label.size, repel = TRUE) + 
                   ggplot2::ggtitle(title))
        
        new.ident <- cluster.list[[i]]$Clusters
        names(new.ident) <- cluster.list[[i]]$rn
        Idents(seurat.obj) <- new.ident
        
        produceMarkers(seurat.obj,
                       cells.per.ident = 50, 
                       top.gene.plot = TRUE)
        
        if(isTRUE(feature.plot)){ 
          print(FeaturePlot(seurat.obj, cols = c("grey", "red"),
                      features = feature.plot,
                      reduction = "umap", ncol = plot.cols))
        }
        
        if(isTRUE(ident.plot)){ 
          ident.df <- data.frame(table(seurat.obj@active.ident, 
                                           seurat.obj@meta.data[, "orig.ident"]))
          
          ggplot(data=ident.df, aes(x=Var2, y=Freq, fill = Var1)) +
            geom_bar(stat="identity", color="black", position = 'fill') + 
            labs(x="Condition", y="Proportion") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
          
          }
        
    } 
} 
