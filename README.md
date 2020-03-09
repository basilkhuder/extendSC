scRNA-Seq R Functions
================
#### Basil Khuder

An array of R functions that extend the functionalities of Seurat
and scRNA-Seq. All parameters are added to a *parameters.json* file
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
               merge = TRUE)
```

### ***featureFiltration.R***:

Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes. If sample is made up of individuals, each individual is filtered separately based upon specific probability quantiles. Default parameters for quantiles are: 

```
    ## $mito.low
    ## [1] 0.02
    ## 
    ## $mito.high
    ## [1] 0.975
    ## 
    ## $feature.cut
    ## [1] 0.975
    ## 
```

``` r
featureFiltration(object = seurat.obj, 
                  parameters = parameters)                      
```
