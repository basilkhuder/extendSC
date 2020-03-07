scRNA-Seq R Functions
================

An array of R functions that extend the functionalities of popular
scRNA-Seq packages. All parameters are added to a *parameters.json* file
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
    ## [1] 0.05
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
(raw/filtered/analysis) and returns a list of Seurat objects. Each
sample should have an individual folder.

    crToSeurat(directory = "directory",
               sample.names = c("sample.one,"sample.two")
               merge = TRUE)
