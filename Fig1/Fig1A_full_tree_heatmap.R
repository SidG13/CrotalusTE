  library(tidyverse)
  library(ggtree)
  library(ape)
  library(RColorBrewer)
  library(ggnewscale)
  
setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig1')
  
  #### All SVSP Tip100 (Cv1-hAT-Tip100 + Cat1-hAT-Tip100) Tree ####
  # Load the tree
  tree <- read.tree('../z_data/All_Cat1_and_Cv1_hAT-Tip100.nwk')
  
  # Root at consensus
  tree <- root(tree, outgroup = "EMBOSS0001_Cv1-Cat1_hAT-Tip100_Consensus", resolve.root = TRUE)
  
  bootstrap_values <- as.numeric(sub(".*/", "", tree$node.label))
  tree$bootstrap <- bootstrap_values
  
  # test collapse at bootstrap 70
  tree$node.label <- ifelse(bootstrap_values < 70, NA, tree$node.label)
  
  # Collapse nodes with low bootstrap values into polytomies
  collapsed_tree <- di2multi(tree)
  
  ggtree_obj <- ggtree(collapsed_tree, branch.length='none', open.angle = 45, size = 0.1)
  
  ggtree_obj$data <- ggtree_obj$data %>%
    mutate(bootstrap = ifelse(isTip, NA, bootstrap_values))  # Only for internal nodes
  
  # Convert names in the tree
  name_conv <- read.table('../z_data/Ultimate_hAT-Tip100_name_conversion.txt', header = T)
  
  ggtree_obj$data <- ggtree_obj$data %>% 
    left_join(name_conv %>% select(id, short_name), by = c('label' = 'short_name'))
  
  
  ggtree_obj$data <- ggtree_obj$data %>%
    mutate(group = case_when(
      grepl("Cv1", id) ~ "croVir",
      grepl("Cat1", id) ~ "croAtr",
    )) %>% 
    mutate(group2 = case_when(
      grepl("SVSP", id) ~ "SVSP",
      TRUE ~ 'nonSVSP'
    )) %>% 
    mutate(group3 = case_when(
      grepl("_E", id) ~ "enhancer",
      grepl("_P", id) ~ "promoter",
      grepl("_O", id) ~ "other",
      TRUE ~ "NonVenom"
    ))
  
  tip_ann <- ggtree_obj$data %>%
    filter(isTip) %>%
    select(label, group, group2, group3) %>%
    column_to_rownames("label")
  
  
  ann_colors <- list(
    group  = c("croVir" = "#199D77",
               "croAtr" = "#E4AB23"),
    
    group2 = c("SVSP" = "#70CACA",
               "nonSVSP" = "grey60"),
    
    group3 = c("enhancer" = brewer.pal(9, "Set1")[4],
               "promoter" = brewer.pal(9, "Set1")[2],
               "other" = "grey80",
               "NonVenom" = "white")
  )
  
p <- ggtree(collapsed_tree,
         branch.length = "none",
         open.angle = 45,
         size = 0.1) +
    theme_tree()
gheatmap(p, tip_ann) +
  scale_fill_manual(
    values = c(
      "croVir"   = "#199D77",
      "croAtr"   = "#E4AB23",
      "SVSP"     = "#70CACA",
      "nonSVSP"  = "grey60",
      "enhancer" = brewer.pal(9, "Set1")[4],
      "promoter" = brewer.pal(9, "Set1")[2],
      "other"    = "grey80",
      "NonVenom" = "white"
    ),
    na.value = "white"
  )
