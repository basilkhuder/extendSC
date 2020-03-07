# scRNA-Seq R-Functions
An array of R functions that extend the functionalities of popular scRNA-Seq packages. All parameters are added to a *parameters.json* file that is read by invoking the following:
```
parameters <- rjson::fromJSON(file = "parameters.json")
parameters[["insert_name"]]
```
The default parameters are: 
```
{ 
"min.cells": 3,
"min.features": 200,
"mito.low": 0.050,
"mito.high": 0.975,
"feature.cut": 0.975,
"cluster_resolution": 0.2,
} 

```

### **_crToSeurat.R_**: 
Takes a directory with the standard CellRanger counts output (raw/filtered/analysis) and returns a list of Seurat objects. Each sample should have an individual folder. 
```
crToSeurat(directory = "directory",
           sample.names = c("sample.one,"sample.two")
           merge = TRUE)
```

