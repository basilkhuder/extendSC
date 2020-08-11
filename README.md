extendSC
================
*By: Basil Khuder*

extendSC adds functionalities to scRNA-Seq workflows using the Seurat package. 

### Installation:
```r
devtools::install_github("basilkhuder/extendSC")
```

### ***crToSeurat.R***:
**Only Cell Ranger >= 3.0 output is supported**

Takes a directory with the standard CellRanger counts output
(raw/filtered/analysis) and returns a list of Seurat objects or a merged Seurat object with ```merge = TRUE```. Each
sample should have an individual folder, as provided in the following directory structure example:
- directory
    - sample.one
        - analysis
        - filtered_gene_bc_matrices
        - raw_gene_bc_matrices
    - sample.two
         - analysis
         - filtered_gene_bc_matrices
         - raw_gene_bc_matrices
         
``` r
crToSeurat(directory = "directory",
           sample.names = c("sample.one","sample.two"),
           merge = TRUE,
           min_cells = 3,
           min_features = 200)
```

### ***featureFiltration.R***:

Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes. If an object is made up of multiple samples, each sample is filtered separately based upon specific probability quantiles. Default quantile parameters are: 

``` r
featureFiltration(object = seurat.obj, 
                  mito.low = .02,
                  mito.high = .975
                  feature.cut = .975)                      
```
### ***processSeurat.R***:

Pipes together several downstream Seurat steps including variance-stabilizing transformation, PCA, clustering and nonlinear dimensionality reduction. Also allows for parallelization using the ```future``` package.

``` r
processSeurat(object = seurat.obj,
              dims = 1:50,
              cluster.res = .2,
              seed.use = 24,
              n.cores = NULL)               
```

### ***umapCellAnno.R***:

Takes a processed Seurat object and returns a ggplot UMAP embedding with cell counts as labels or counts in the legend. Specify the name of the slot that contains your Seurat clusters with ```annotation.name```. Default is *Seurat_Assignment* followed by *seurat_clusters*

``` r
umapCellAnno(seurat.obj,
             annotation.name,
             point.size = 1,
             label.size = 8,
             title = "",
             counts.as.title = FALSE,
             counts.as.labels = FALSE,
             legend = TRUE,
             counts.in.legend = TRUE,
             use.cols = "")       
```

### ***chooseClusterRes.R***:
```chooseClusterRes()``` allows you to explore various cluster resolution with the outputs being a UMAP, feature plots and identity plots. Unlike ```processSeurat()```, ```chooseClusterRes()``` does not add the range of cluster resolutions to the seurat object.  Cluster resolutions must be in a list with the first element being "from", the second being "to" and the third being "by." 
```r
chooseClusterRes(seurat.obj, 
                 cluster.res = list(.1,.5,.1),
                 ident.plot = TRUE,
                 feature.plot = c("Gene1","Gene2","Gene3"),
                 plot.cols = 1)                           
```

### ***dendoSeurat.R***:
```dendoSeurat()``` allows you to hierarchically cluster a subset of a Seurat object and return a dendrogram. Due to the sheer computational demands of hierarchical clustering, down-sampling of your Seurat object is required. If you'd like both a dendrogram and the hierarchical clustering labels returned, set ```return.clusters``` to the height of the tree you want cut to produce the groupings. 
```r
chooseClusterRes(seurat.obj, 
                 cluster = "T Cells",
                 annotation.name = "Seurat_Assignment",
                 down.sample = 50,
                 variable.genes = NULL,
                 seed = 1,
                 return.clusters = 4)                           
```
