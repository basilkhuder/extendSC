#' Takes a CellRanger filtered output and returns a matrix
#' @param directory The path to the directory that contains the filtered_gene_matrices folder
#' @return A filtered gene matrix
#' @examples
#' directory <- "/home/users/documents/CRoutput"
#' matrix <- crToMTX(directory) 
#' @export

crToMTX <- function(directory) {
  path <- paste(directory, list.files(path = directory,
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
