# Takes a directory with standard CellRanger counts output (raw/filtered/analysis) and returns a list of Seurat objects.
# With the merge parameter, object list can be merged with sample ids provided

crToSeurat <- function(directory, 
                       cells = 3,
                       features = 200,
                       sample.names = NULL,
                       merge = FALSE){ 
  
  `%>%` <- magrittr::`%>%`
  folders <- list.dirs(directory, recursive = FALSE)
  matrix.list <- vector(mode = "list", length = length(folders))
  
  for (i in 1:length(folders)){ 
    filtered.folder <- list.files(path = folders[[i]], pattern = "filtered")
    full.dir <- paste(folders[[i]],filtered.folder,sep = "/")
    matrix.list[[i]] <- Matrix::readMM(paste(full.dir,list.files(path = full.dir, pattern = "matrix"),sep = "/"))
    features <- read.delim(paste(full.dir,list.files(path = full.dir, pattern = "features"),sep = "/"), 
                           header = FALSE,
                           stringsAsFactors = FALSE)
    barcodes <- read.delim(paste(full.dir,list.files(path = full.dir, pattern = "barcode"),sep = "/"),
                           header = FALSE,
                           stringsAsFactors = FALSE)
    matrix.list[[i]] <- matrix.list[[i]] %>%
      set_colnames(barcodes$V1) %>%
      set_rownames(features$V2) %>%
      CreateSeuratObject(min.cells = cells, 
                         min.features = features,
                         project = sample.names[i])
  }
  if (isTRUE(merge)){
    matrix.list <- merge(matrix.list[[1]], y = matrix.list[-1], add.cell.ids = sample.names)
  } 
  return(matrix.list)
}
