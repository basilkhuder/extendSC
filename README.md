#### An array of R functions that extend the functionalities of Seurat, SingleR, and other scRNA-Seq packages. 
#### Basil Khuder
---------------------------
extendSC adds functionalities to scRNA-Seq workflows for Seurat and SingleR. For some of the functions (```crToSeurat()```), parameters are added to a *parameters.json* file
that is read by invoking the following:

``` r
parameters <- rjson::fromJSON(file = "parameters.json")
parameters[["min.cells"]]
```

    ## [1] 3

The default parameters are:

``` r
parameters
```

    ## $min.cells
    ## [1] 3
    ## 
    ## $min.features
    ## [1] 200
    ## 
    ## $mito.low
    ## [1] 0.02
    ## 
    ## $mito.high
    ## [1] 0.975
    ## 
    ## $feature.cut
    ## [1] 0.975
    ## 
    ## $cluster_resolution
    ## [1] 0.2

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
           parameters = parameters)
```

### ***featureFiltration.R***:

**_Update (3/11/20):_** *```featureFiltration()```no longer uses the JSON parameters file. Instead, ```mito.low```, ```mito.high``` and ```feature.cut``` are parameters set within the function.*

Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes. If sample is made up of individuals, each individual is filtered separately based upon specific probability quantiles. Default quantile parameters are: 

``` r
featureFiltration(object = seurat.obj, 
                  mito.low = .02,
                  mito.high = .975
                  feature.cut = .975)                      
```

### ***umapCellAnno.R***:

```umapCellAnno()``` takes a processed Seurat object and returns a ggplot UMAP embedding with cell counts as labels. 

``` r
umapCellAnno(object = seurat.obj,
             point.size = 1,
             label.size = 8,
             counts.as.title = TRUE,
             legend.pos = "right")            
```
