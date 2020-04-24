#' @export

crToSeurat <- function(directory, 
                       min_cells = 3,
                       min_features = 200, 
                       sample_names = c("sample"), 
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
      magrittr::set_colnames(barcodes$V1) %>%
      magrittr::set_rownames(features$V2) 
  }
  
  matrix.list <- lapply(seq_along(matrix.list), function(x) 
    Seurat::CreateSeuratObject(matrix.list[[x]], min.cells = min_cells,
                               min.features = min_features,
                               project = sample_names[x]))
  
  if (length(matrix.list) == 1){ 
    matrix.list <- matrix.list[[1]]
    }
  
  
  if (isTRUE(merge)){
    matrix.list <- merge(matrix.list[[1]], 
                         matrix.list[-1],
                         add.cell.ids = sample_names)
  }
  return(matrix.list)
}
