#' Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes
#' @param seurat.obj A seurat object. 
#' @param mito.low Filter cells that have total mitochondrial genes low than this quantile
#' @param mito.high Filter cells that have total mitochondrial genes higher than this quantile 
#' @param feature.cut Filter cells that have total features higher than this quantile 
#' @examples
#' featureFiltration(seurat.obj, mito.low = .5, mito.high = .975, feature.cut = .975)
#' @export

featureFiltration <- function(seurat.obj,
                              mito.low = .05,
                              mito.high = .975,
                              feature.cut = .975) {
  seurat.obj <- PercentageFeatureSet(seurat.obj,
                                     pattern = "^MT-",
                                     col.name = "percent.mt")
  
  ident <- as_tibble(FetchData(seurat.obj,vars = c("orig.ident",
                                         "percent.mt",
                                         "nFeature_RNA",
                                         "nCount_RNA")), rownames = "Cells")
  
  quantile.amounts <- ident %>% 
    group_by(orig.ident) %>% 
    summarize(mito.low = quantile(percent.mt, probs = mito.low),
              mito.high = quantile(percent.mt, probs = mito.high),
              feature.cut = quantile(nFeature_RNA, probs = feature.cut))
  
  ident.filtered <- ident %>% group_split(orig.ident) %>%
    imap_dfr(~ filter(.x, percent.mt > quantile.amounts$mito.low[[.y]], 
                      percent.mt < quantile.amounts$mito.high[[.y]],
                      nFeature_RNA < quantile.amounts$feature.cut[[.y]]))
  
  return(subset(seurat.obj, cells = ident.filtered$Cells))
}
