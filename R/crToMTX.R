#' @export

crToMTX <- function(folders) {
  
  path <- paste(folders, list.files(path = folders,
                                    pattern = "filtered"), sep = "/")
  mtx <- paste(path, list.files(path = path,
                                pattern = "matrix"), sep = "/")
  mtx <- Matrix::readMM(mtx)
  
  features <- paste(path, list.files(path = path,
                                     pattern = "features"), sep = "/")
  features <- read.delim(features, header = FALSE, stringsAsFactors = FALSE)$V2
  barcodes <- paste(path, list.files(path = path,
                                     pattern = "barcodes"), sep = "/")
  barcodes <- read.delim(barcodes, header = FALSE, stringsAsFactors = FALSE)$V1
  rownames(mtx) <- features
  colnames(mtx) <- barcodes
  return(mtx)
} 
