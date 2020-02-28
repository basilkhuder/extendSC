# Takes each folder of CellRanger counts output (raw/filtered/analysis) and generates a Seurat object

crToSeurat <- function(directory){ 
  
  folders <- list.dirs(directory,recursive = FALSE)
  matrix.list <- vector(mode = "list", length = length(folders))
  features.list <- vector(mode = "list", length = length(folders))
  barcode.list <- vector(mode = "list", length = length(folders))
  
  for (i in 1:length(folders)){ 
    filtered.folder <- list.files(path = folders[[i]], pattern = "filtered")
    full.dir <- paste(folders[[i]],filtered.folder,sep = "/")
    
    matrix.list[[i]] <- readMM(paste(full.dir,list.files(path = full.dir, pattern = "matrix"),sep = "/"))
    
    features.list[[i]] <- read.delim(paste(full.dir,list.files(path = full.dir, pattern = "features"),sep = "/"), 
                                     header = FALSE,
                                     stringsAsFactors = FALSE)
    
    barcode.list[[i]] <- read.delim(paste(full.dir,list.files(path = full.dir, pattern = "barcode"),sep = "/"),
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
    
    colnames(matrix.list[[i]]) <- barcode.list[[i]]$V1
    rownames(matrix.list[[i]]) <- features.list[[i]]$V1
    matrix.list[[i]] <- CreateSeuratObject(matrix.list[[i]])
  }
  
  return(matrix.list)
  
}

