#' @export
produceQCPlots <- function(ident,
                           mito.high,
                           mito.low,
                           feat.cut){ 
    
    if(class(ident) == "Seurat") {
        ident <- data.table::setDT(Seurat::FetchData(ident,
                                                     vars = c(
                                                         "orig.ident",
                                                         "percent.mt",
                                                         "nCount_RNA",
                                                         "nFeature_RNA"
                                                     )),
                                   keep.rownames = TRUE)
    }
    ident.names <- unique(ident$orig.ident)
    ident.list <- lapply(seq_along(ident.names), 
                         function(x) dplyr::filter(ident, orig.ident == ident.names[x]))
    plot1 <- lapply(seq_along(ident.list), function(x)
        ggplot2::ggplot(data = ident.list[[x]], mapping = ggplot2::aes(x = nCount_RNA, y = percent.mt)) +
            ggplot2::geom_point(color = "#f8766d") + 
            ggplot2::geom_hline(yintercept=mito.low[[x]], linetype="dashed", color = "red") +
            ggplot2::geom_hline(yintercept=mito.high[[x]], linetype="dashed", color = "red") +
            ggplot2::ggtitle(ident.names[[x]]))
    plot2 <- lapply(seq_along(ident.list), function(x)
        ggplot2::ggplot(data = ident.list[[x]], mapping = ggplot2::aes(x = nCount_RNA, y = nFeature_RNA)) +
            ggplot2::geom_point(color = "#f8766d") + 
            ggplot2::geom_hline(yintercept=feat.cut[[x]], linetype="dashed", color = "red") +
            ggplot2::ggtitle(ident.names[[x]]))
    print(plot1)
    print(plot2)
}
