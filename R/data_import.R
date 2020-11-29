#' Takes a directory CellRanger counts output and returns a list of Seurat objects or a merged Seurat object
#' @param seurat.obj A Seurat object. 
#' @param min_cells Include features detected in at least this many cells.
#' @param min_features Include cells where at least this many features are detected 
#' @param sample.names Names for each of the samples
#' @param gz Whether or not the files are gz compressed
#' @return A Seurat object 
#' @examples
#' crToSeurat(directory = "Files", min_cells = 3, min_features = 3)
#' @export

crToSeurat <- function(directory, 
                       min.cells = 3,
                       min.features = 200, 
                       sample.names = c("sample"), 
                       gz = FALSE, 
                       merge = TRUE){ 
  
  components <- c("matrix.mtx", "features.tsv", "barcodes.tsv")
  
  if (gz) {
    components <- map_chr(components, ~str_c(.,".gz"))
  } 
  
  path <- str_c(list.files(directory, full.names = TRUE), "filtered_gene_bc_matrices",
                sep = "/")
  
  matrix <- map(str_c(path,components[[1]],sep = "/"), readMM)
  
  features <- map(str_c(path,components[[2]],sep = "/"), ~
                    read_tsv(.x, col_names = FALSE, col_types = cols())$X2)
  
  barcodes <- map(str_c(path,components[[3]],sep = "/"), ~
                    read_tsv(.x, col_names = FALSE, col_types = cols())$X1) 
  
  matrix <- imap(matrix, ~ set_colnames(.x, barcodes[[.y]])) 
  matrix <- imap(matrix, ~ set_rownames(.x, features[[.y]]))
  matrix <- imap(matrix, ~ CreateSeuratObject(.x, 
                                      min.cells = min.cells,
                                      project = sample.names[[.y]],
                                      min.features = min.features))
  
  if (length(matrix) == 1) { 
    matrix <- matrix[[1]]
  }
  
  if (merge) {
    matrix <- merge(matrix[[1]], 
                    matrix[-1],
                    add.cell.ids = sample.names)
  }
  return(matrix)
}

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
