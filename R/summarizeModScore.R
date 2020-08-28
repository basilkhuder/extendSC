summarizeModScore <- function(seurat.obj,
                              features,
                              annotation.name) { 
  
  mod.scores <- extract2(AddModuleScore(object = bal.mp,features = list(features),
                                        name = "module_score"),
                         "module_score1") %>%
  as_tibble(rownames = "Cells")
  
  clusters <- extract2(bal.mp, annotation.name) %>%
    as_tibble(rownames = "Cells") %>%
    full_join(mod.scores)
  
  print(clusters %>% 
          group_by(!!as.name(annotation.name)) %>% 
          rstatix::get_summary_stats())
  } 



