#' Takes a CellRanger filtered output and returns a gene matrix
#' @param directory The path to the directory that contains the filtered_gene_matrices folder
#' @return A filtered gene matrix
#' @examples
#' directory <- "/home/users/documents/CRoutput"
#' matrix <- crToMTX(directory) 
#' @export

crToMTX <- function(directory) {
  path <- list.files(path = directory, pattern = "filtered", full.names = TRUE)
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
  features <- read_tsv(features, col_names = TRUE)$X2
  barcodes <- paste(path, list.files(path = path,
                                     pattern = "barcodes"), sep = "/")
  barcodes <- read_tsv(barcodes, col_names = TRUE)$X1
  rownames(mtx) <- features
  colnames(mtx) <- barcodes
  return(mtx)
} 
