# extendSC
#### An array of R functions that extend the functionalities of Seurat, SingleR, and other scRNA-Seq packages. 
#### Basil Khuder
---------------------------
extendSC adds functionalities to scRNA-Seq workflows for the Seurat and SingleR packages. 

### ***crToSeurat.R***:

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
           sample.names = c("sample.one,"sample.two"),
           merge = TRUE,
           min_cells = 3,
           min_features = 200)
```

### ***featureFiltration.R***:

Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes. If sample is made up of individuals, each individual is filtered separately based upon specific probability quantiles. Default quantile parameters are: 

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
              cluster.res = 0.2,
              seed.use = 24,
              n.cores = NULL)               
```

### ***umapCellAnno.R***:

Takes a processed Seurat object and returns a ggplot UMAP embedding with cell counts as labels or counts in the legend. 

``` r
umapCellAnno(seurat.obj,
             point.size = 1,
             label.size = 8,
             title = "",
             counts.as.title = FALSE,
             counts.as.labels = FALSE,
             legend = TRUE,
             counts.in.legend = TRUE,
             use.cols = "")       
```
### ***seuratToSingleR.R***:
```seuratToSingleR()``` takes an annotated Seurat object to be used as a SingleR reference, any ```SingleCellExperiment``` object you're interested in annotating, and returns either a ```PlotScoreHeatmap()``` or a table with annotations. 

``` r
seuratToSingleR(reference = seurat.obj,
                sce.obj = sce.obj,
                heatmap = TRUE,
                table = TRUE)
```


