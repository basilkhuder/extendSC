# Filters a Seurat object based upon mitochondrial and feature percentages. Instead of the traditional approach
# of arbitrarily filtering on a fixed %, filtering is based upon quantile percentages. 
# If seurat.obj is made up of multiple samples, each sample must have unique orig.ident

featureFiltration <- function(seurat.obj,
                              mito.low = .02,
                              mito.high = .975,
                              feature.cut = .975){
  seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^MT-", col.name = "percent.mt")
  ident <- data.table::setDT(Seurat::FetchData(seurat.obj, vars = c("orig.ident","percent.mt","nFeature_RNA")),
                             keep.rownames = TRUE) 
  ident.list <- as.list(unique(ident$orig.ident))
  mito.low <- purrr::map(ident.list, ~quantile(ident[ident$orig.ident == .x,]$percent.mt,
                                               probs = mito.low))
  mito.high <- purrr::map(ident.list, ~quantile(ident[ident$orig.ident == .x,]$percent.mt,
                                                probs = mito.high))
  feat.cut <- purrr::map(ident.list, ~quantile(ident[ident$orig.ident == .x,]$nFeature_RNA,
                                               probs = feature.cut))
  filter.list <- vector(mode = "list", length = length(ident.list))
  for (i in 1:length(ident.list)){ 
    ident.current <- ident[ident$orig.ident == ident.list[[i]]]
    ident.current <- ident.current[ident.current$percent.mt > mito.low[[i]]]
    ident.current <- ident.current[ident.current$percent.mt < mito.high[[i]]]
    ident.current <- ident.current[ident.current$nFeature_RNA < feat.cut[[i]]]
    filter.list[[i]] <- ident.current
  }
  return(subset(seurat.obj,cells = do.call(rbind, filter.list)$rn))
}
