extractMeta <- function(seurat.obj, 
                        types, 
                        vars, 
                        merge = TRUE) { 
  
  if (length(types) != length(vars)) { 
    stop("The amount of types need to be equal to the amount of variables.")
  }
  
  types <- str_to_lower(types)
  meta.list <- vector(mode = "list", length = length(types))
  
  if (any(types %in% "identity")) { 
    identity <- vars[which(types %in% "identity")]
    identity.df <- map(identity, ~as_tibble(FetchData(seurat.obj, vars = .x), rownames = "Cells"))
    for(i in seq_along(identity.df)){ 
      meta.list[[which(types %in% "identity")[[i]]]] <- identity.df[[i]]
    } 
    
  } 
  
  if (any(types %in% "clusters")) { 
    clusters <- vars[which(types %in% "clusters")]
    clusters.df <- map(clusters, ~as_tibble(FetchData(seurat.obj, vars = .x), rownames = "Cells"))
    for(i in seq_along(clusters.df)) { 
      meta.list[[which(types %in% "clusters")[[i]]]] <- clusters.df[[i]]
    } 
    
  } 
  
  if (any(types %in% "module scores")) { 
    module.scores <- vars[which(types %in% "module scores")]
    module.scores.df <- map(module.scores, ~as_tibble(FetchData(seurat.obj, vars = .x), rownames = "Cells"))
    for(i in seq_along(module.scores.df)) { 
      meta.list[[which(types %in% "module scores")[[i]]]] <- module.scores.df[[i]]
    } 
  }
  
  if (any(types %in% "embeddings")) { 
    embeddings <- vars[which(types %in% "embeddings")]
    embeddings.df <- map(embeddings, ~as_tibble(Embeddings(seurat.obj, reduction = .x), rownames = "Cells"))
    for(i in seq_along(embeddings.df)) { 
      meta.list[[which(types %in% "embeddings")[[i]]]] <- embeddings.df[[i]]
    } 
    
  }
  
  if (merge) { 
    return(purrr::reduce(meta.list, full_join))
  } 
  return(meta.list)
  
} 

#' Easily extract counts from a Seurat object
#' @param seurat.obj A seurat object. 
#' @param assay The assay to pull from
#' @param genes Optional argument to extract certain genes only
#' @param tibble Whether you want counts returned as a tibble
#' @return Counts
#' @examples
#' extractCounts(seurat.obj, assay = "RNA", genes = c("TBC1D3D","LINC00470"), tibble = TRUE) 
#' @export

extractCounts <- function(seurat.obj,
                          assay = "SCT",
                          genes = NULL,
                          tibble = FALSE) { 
  
  counts <- GetAssayData(seurat.obj, assay = assay)
  if(!is.null(genes)) {
    counts <- counts[map_dbl(genes, ~ which(rownames(counts) %in% .x)), ]
  }
  
  if (tibble) { 
    if (length(genes) == 1) { 
      counts <- as_tibble(counts, rownames = "Cells") %>%
        mutate(Genes = genes, .before = Cells)
      return(counts)
    }
    counts <- as_tibble(counts, rownames = "Genes") %>%
      pivot_longer(cols = -Genes) %>%
      rename(Cells = name)
  }
  return(counts)
}

#' @export
produceMarkers <- function(seurat.obj,
                           cells.per.ident = Inf,
                           top.gene.plot = TRUE,
                           output.name = NULL,
                           test.use = "wilcox") { 
  
  markers <- FindAllMarkers(object = seurat.obj, 
                            only.pos = TRUE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25, 
                            max.cells.per.ident = cells.per.ident,
                            test.use = test.use)
  
  top3 <- markers %>% 
    group_by(cluster) %>% 
    slice_max(n = 3, order_by = avg_logFC) %>%
    select(gene, everything())
  print(top3)
  
  file.name <- if_else(is.null(output.name),
                       str_c(str_replace_all(Sys.Date(), "-", "_"), "_markers.txt"),
                       output.name)
  write_tsv(markers, path = file.name)
  
  if (top.gene.plot) {
    cluster.averages <- AverageExpression(object = seurat.obj, return.seurat = TRUE)
    print(pheatmap(GetAssayData(cluster.averages)[unique(top3$gene), ]))
  }
  return(markers)
}
