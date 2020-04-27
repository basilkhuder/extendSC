featureFiltration <- function(seurat.obj,
                              mito.low = .05,
                              mito.high = .975,
                              feature.cut = .975,
                              produce_plots = TRUE) {
    seurat.obj <- Seurat::PercentageFeatureSet(seurat.obj,
                                               pattern = "^MT-",
                                               col.name = "percent.mt")
    ident <- data.table::setDT(Seurat::FetchData(
        seurat.obj,
        vars = c("orig.ident",
                 "percent.mt",
                 "nFeature_RNA",
                 "nCount_RNA")
    ), keep.rownames = TRUE)
    
    ident.list <- as.list(unique(ident$orig.ident))
    mito.low <- lapply(ident.list, function(x) 
        quantile(ident[ident$orig.ident == x, ]$percent.mt, probs = mito.low))
    
    mito.high <- lapply(ident.list, function(x)
        quantile(ident[ident$orig.ident == x, ]$percent.mt, probs = mito.high))
    
    feat.cut <- lapply(ident.list, function(x)
        quantile(ident[ident$orig.ident == x, ]$nFeature_RNA,
                              probs = feature.cut))
    
    if (isTRUE(produce_plots)) {
        produceQCPlots(ident, mito.low, mito.high, feat.cut)
    }
    
    filter.list <- vector(mode = "list", length = length(ident.list))
    for (i in 1:length(ident.list)) {
        ident.current <- ident[ident$orig.ident == ident.list[[i]]]
        ident.current <-
            ident.current[ident.current$percent.mt > mito.low[[i]]]
        ident.current <-
            ident.current[ident.current$percent.mt < mito.high[[i]]]
        ident.current <-
            ident.current[ident.current$nFeature_RNA < feat.cut[[i]]]
        filter.list[[i]] <- ident.current
    }
    return(Seurat::subset(seurat.obj, cells = do.call(rbind, filter.list)$rn))
}
