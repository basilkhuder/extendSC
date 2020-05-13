#' @export

crToMTX <- function(folders) {
  path <- paste(folders, list.files(path = folders,
                                    pattern = "filtered"), sep = "/")
  mtx <- paste(path, list.files(path = path,
                                pattern = "matrix"), sep = "/")
  mtx <- Matrix::readMM(mtx)
  features <- list.files(path = path, pattern = "features")
  if(isTRUE(identical(features, character(0)))){ 
    features <- list.files(path = path, pattern = "gene")
  }
  if(isTRUE(identical(features, character(0)))){ 
    stop("Features/Gene file not found. If file exists, rename as 'features.tsv' or 
         'genes.tsv'")
  }
  features <- paste(path, features, sep = "/")
  features <- read.delim(features, header = FALSE, stringsAsFactors = FALSE)$V2
  barcodes <- paste(path, list.files(path = path,
                                     pattern = "barcodes"), sep = "/")
  barcodes <- read.delim(barcodes, header = FALSE, stringsAsFactors = FALSE)$V1
  rownames(mtx) <- features
  colnames(mtx) <- barcodes
  return(mtx)
} 
