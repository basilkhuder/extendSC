# Takes a directory with standard CellRanger counts output (raw/filtered/analysis) and generates a 
# Seurat object based upon the filtered CR counts

crToSeurat <- function(directory){ 
  
  folders <- list.dirs(directory, recursive = FALSE)
  matrix.list <- vector(mode = "list", length = length(folders))
  
  for (i in 1:length(folders)){ 
    filtered.folder <- list.files(path = folders[[i]], pattern = "filtered")
    full.dir <- paste(folders[[i]],filtered.folder,sep = "/")
    matrix.list[[i]] <- readMM(paste(full.dir,list.files(path = full.dir, pattern = "matrix"),sep = "/"))
    features <- read.delim(paste(full.dir,list.files(path = full.dir, pattern = "features"),sep = "/"), 
                                     header = FALSE,
                                     stringsAsFactors = FALSE)
    barcodes <- read.delim(paste(full.dir,list.files(path = full.dir, pattern = "barcode"),sep = "/"),
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
    matrix.list[[i]] <- matrix.list[[i]] %>%
                        set_colnames(barcode$V1) %>%
                        set_rownames(features$V1) %>%
                        CreateSeuratObject()
  }
  return(matrix.list)
}
