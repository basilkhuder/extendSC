# Annotate clusters on UMAP embedding with cell numbers

function <- umapAnno(seurat.obj) { 

extract.clusters <- setDT(FetchData(seurat_obj, vars = c("seurat_clusters")), 
                          keep.rownames = TRUE)
                          
cluster.counts <- extract.clusters %>%
                  group_by(seurat_clusters) %>%
                  tally()

new.cluster.ids <- condition.extract.count$n
names(x = new.cluster.ids) <- levels(x = seurat_obj)
seurat_obj <- RenameIdents(object = seurat_obj, new.cluster.ids)
DimPlot(object = seurat_obj, reduction = 'umap', pt.size =1, label = TRUE, label.size = 8)
seurat_obj@active.ident <- seurat_obj$seurat_clusters
} 
