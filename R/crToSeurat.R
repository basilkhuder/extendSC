#' Takes a directory CellRanger counts output and returns a list of Seurat objects or a merged Seurat object
#' @param seurat.obj A Seurat object. 
#' @param min_cells Include features detected in at least this many cells.
#' @param min_features Include cells where at least this many features are detected 
#' @return A Seurat object 
#' @examples
#' crToSeurat(directory = "Files", min_cells = 3, min_features = 3)
#' @export

crToSeurat <- function(directory, 
                       min_cells = 3,
                       min_features = 200, 
                       sample_names = c("sample"), 
                       merge = FALSE){ 
  
  path <- str_c(directory, list.files(directory), "filtered_gene_bc_matrices",
        sep = "/")
  
  matrix <- map(str_c(path,"matrix.mtx",sep = "/"), readMM)
  features <- map(str_c(path,"features.tsv",sep = "/"), ~
                  read_tsv(.x, col_names = FALSE)$X2)
  
  barcodes <- map(str_c(path,"barcodes.tsv",sep = "/"), ~
                       read_tsv(.x, col_names = FALSE)$X1) 
  
  matrix <- imap(matrix, ~ set_colnames(.x, barcodes[[.y]]))
  matrix <- imap(matrix, ~ set_rownames(.x, features[[.y]]))
  matrix <- imap(matrix, ~ CreateSeuratObject(.x, min_cells = min_cells,
                                              sample_names = sample_names[[.y]],
                                              min_features = min_features))
  
  if (length(matrix) == 1){ 
    matrix <- matrix[[1]]
  }
  
  if (isTRUE(merge)){
    matrix <- merge(matrix[[1]], matrix[-1], add.cell.ids = sample_names)
  }
  return(matrix)
}
