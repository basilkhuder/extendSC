#' @export
produceQCPlots <- function(ident,
                           mito.high,
                           mito.low,
                           feat.cut) {
    ident.names <- unique(ident$orig.ident)
    ident.list <- lapply(seq_along(ident.names),
                         function(x)
                             dplyr::filter(ident, orig.ident == ident.names[x]))
    
    for (i in 1:length(unique(ident$orig.ident)) {
        plot1 <-
            ggplot2::ggplot(data = ident.list[[i]],
                            mapping = ggplot2::aes(x = nCount_RNA, y = percent.mt)) +
            ggplot2::geom_point(color = "#f8766d") +
            ggplot2::geom_hline(yintercept = mito.low[[i]],
                                linetype = "dashed",
                                color = "red") +
            ggplot2::geom_hline(yintercept = mito.high[[i]],
                                linetype = "dashed",
                                color = "red") +
            ggplot2::ggtitle(ident.names[[i]]))
    
        plot2 <- ggplot2::ggplot(data = ident.list[[i]], 
                                 mapping = ggplot2::aes(x = nCount_RNA, y = nFeature_RNA)) +
            ggplot2::geom_point(color = "#f8766d") +
            ggplot2::geom_hline(yintercept = feat.cut[[i]], linetype = "dashed", color = "red") +
            ggplot2::ggtitle(ident.names[[i]]))
    }
    
    print(plot1)
    print(plot2)
}

