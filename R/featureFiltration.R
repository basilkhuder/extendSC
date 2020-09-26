#' Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes
#' @param seurat.obj A seurat object. 
#' @param mito.low Filter cells that have total mitochondrial genes low than this quantile
#' @param mito.high Filter cells that have total mitochondrial genes higher than this quantile 
#' @param feature.cut Filter cells that have total features higher than this quantile 
#' @param multiple.samples Whether the seurat.obj is made up of multiple sames (differentiated by orig.ident)
#' @examples
#' featureFiltration(seurat.obj, mito.low = .5, mito.high = .975, feature.cut = .975)
#' @export

featureFiltration <- function(seurat.obj,
                              mito.low = .05,
                              mito.high = .975,
                              feature.cut = .975,
                              multiple.samples = TRUE) {
  seurat.obj <- PercentageFeatureSet(seurat.obj,
                                     pattern = "^MT-",
                                     col.name = "percent.mt")
  
  if (isFALSE(multiple.samples)) {
    ident <- as_tibble(FetchData(
      seurat.obj,
      vars = c("percent.mt",
               "nFeature_RNA",
               "nCount_RNA")
    ), rownames = "Cells")
    
    quantile.amounts <- ident %>%
      summarize(
        mito.low = quantile(percent.mt, probs = mito.low),
        mito.high = quantile(percent.mt, probs = mito.high),
        feature.cut = quantile(nFeature_RNA, probs = feature.cut)
      )
    
    ident.filtered <-
      ident %>% filter(
        percent.mt > quantile.amounts$mito.low,
        percent.mt < quantile.amounts$mito.high,
        nFeature_RNA < quantile.amounts$feature.cut
      )
    
    if(!identical(ident.filtered$Cells, character(0))){
      seurat.obj <- subset(seurat.obj, cells = ident.filtered$Cells)
    }
  
  } else {
    ident <- as_tibble(FetchData(
      seurat.obj,
      vars = c("orig.ident",
               "percent.mt",
               "nFeature_RNA",
               "nCount_RNA")
    ), rownames = "Cells")
    
    quantile.amounts <- ident %>%
      group_by(orig.ident) %>%
      summarize(
        mito.low = quantile(percent.mt, probs = mito.low),
        mito.high = quantile(percent.mt, probs = mito.high),
        feature.cut = quantile(nFeature_RNA, probs = feature.cut)
      )
    
    ident.filtered <- ident %>% group_split(orig.ident) %>%
      imap_dfr(
        ~ filter(
          .x,
          percent.mt > quantile.amounts$mito.low[[.y]],
          percent.mt < quantile.amounts$mito.high[[.y]],
          nFeature_RNA < quantile.amounts$feature.cut[[.y]]
        )
      )
    
    if(!identical(ident.filtered$Cells, character(0))){
      seurat.obj <- subset(seurat.obj, cells = ident.filtered$Cells)
    }
  }
  return(seurat.obj)
}
