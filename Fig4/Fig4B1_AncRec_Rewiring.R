library(tidyverse)
library(ape)
library(phytools)
library(ggtree)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig4')

tree <- read.tree("../z_data/All_Cat1_and_Cv1_hAT-Tip100.nwk")
matrix_data <- read.table('../z_data/NEW_CRE_hATs_rewired_trinary_matrix_for_modelling_Oct2025.txt', header = T)

#### Writing a loop to run simmap for the re-wired TFBSs: RUN ON EXTERNAL SERVER AT UTA ####
trait_cols <- colnames(matrix_data)[5:(5+21-1)] # all 21 re-wired TFBSs
simmap_list <- vector("list", length = length(trait_cols))
names(simmap_list) <- trait_cols

# 
# for(trait in trait_cols) {
#   cat("Starting trait:", trait, "\n")
# 
#   tfbs_vec <- matrix_data %>%
#     select(Tip100, all_of(trait)) %>%
#     deframe()
#   
#   # Recode as binary!!
#   tfbs_vec <- ifelse(tfbs_vec > 0, 1, 0)
#   tfbs_vec <- factor(tfbs_vec, levels = c(0, 1))
#   
#   # Drop tips from the tree not in the trait vector
#   to_drop <- setdiff(tree$tip.label, names(tfbs_vec))
#   tree_sub <- drop.tip(tree, to_drop)
#   
#   # Run SIM map 
#   simmap_list[[trait]] <- make.simmap(tree = tree_sub, x = tfbs_vec,
#                                       model = 'SYM', Q = 'mcmc')
#   
#   cat("Completed trait:", trait, "\n")
# }
#### End Loop ####
# Save all simmaps
#saveRDS(simmap_list, file = "simmaps_55prewired_TFBSs.rds")


#### Plotting individual TFBSs ####
simmap_list <- readRDS('../z_data/Z_simmaps_21rewired_TFBSs_for_rewiring_x100_reps.rds') # Find on Zenodo

# # Pick a specific simmap, e.g., the 3rd TFBS, 1st simulation out of 100
# sim_num <- 51
# simmap_obj <- simmap_list[[which(names(simmap_list) == 'Creb3l2.2')]][[sim_num]]
# trait_name <- names(simmap_list)[which(names(simmap_list) == 'Creb3l2.2')]
# tree_label <- paste0(trait_name, ' sim ', as.character(sim_num))
# 
# svsp_tips <- grep("SVSP", simmap_obj$tip.label, value = TRUE)
# 
# tree_sub_simmap <- drop.tip(simmap_obj, setdiff(simmap_obj$tip.label, svsp_tips))
# 
# # Count transitions in the subtree
# transitions <- countSimmap(tree_sub_simmap)
# print(transitions)
# 
# maps <- tree_sub_simmap$maps
# edges <- 1:length(maps)
# 
# # Get 0->1 gains
# gain_edges <- edges[sapply(maps, function(edge_map) {
#   states <- names(edge_map)
#   states[1] == "0" & any(states == "1")
# })]
# 
# # Get 1->0 losses
# loss_edges <- edges[sapply(maps, function(edge_map) {
#   states <- names(edge_map)
#   states[1] == "1" & any(states == "0")
# })]
# 
# # Plot simmap tree
# plotSimmap(tree_sub_simmap,
#            colors = setNames(c("steelblue", "orange"), c("0", "1")),
#            fsize = 0.7, lwd = 2)
# title(main = paste0('State transitions for "', tree_label, '"'))
# 
# # Add dots for 0->1 gains
# #edgelabels(pch = 21, bg = "red", cex = 1.5, edge = loss_edges)
# edgelabels(pch = 21, bg = "green", cex = 1.5, edge = gain_edges)
# 
# tip_order <- tree_sub_simmap$tip.label
# ann <- matrix_data[match(tip_order, matrix_data$Tip100), ]
# 
# 
# dot_cols <- ifelse(ann$type == "promoter", "blue",
#                    ifelse(ann$type == "enhancer", "orange", NA))
# 
# tiplabels(pch = 25, bg = "blue", cex = 1.2,
#           tip = which(ann$type == "promoter"))
# 
# # enhancer dots
# tiplabels(pch = 25, bg = "orange", cex = 1.2,
#           tip = which(ann$type == "enhancer"))
# 
# trait_vec <- matrix_data[match(tree_sub_simmap$tip.label, matrix_data$Tip100), trait_name]
# trait_vec <- ifelse(trait_vec > 0, 1, 0)
# 
# trait_vec <- as.numeric(trait_vec)
# 
# x_tip <- max(nodeHeights(tree_sub_simmap)) + 0.05  # shift right
# 
# y_tip <- 1:length(tree_sub_simmap$tip.label)
# 
# box_width  <- 0.02
# box_height <- 0.5
# 
# for (i in seq_along(trait_vec)) {
#   col_box <- ifelse(trait_vec[i] == 0, "white",
#                     ifelse(trait_vec[i] == 1, "red", NA))
#   if (!is.na(col_box)) {
#     rect(xleft   = x_tip,
#          ybottom = y_tip[i] - box_height/2,
#          xright  = x_tip + box_width,
#          ytop    = y_tip[i] + box_height/2,
#          col     = col_box,
#          border  = "black")
#   }
# }

# Count the # of CREs affected by motif gains

#### Run a loop or something to get all num_upstream_gain nodes AND total CREs per sim per tfbs ####

get_descendant_tips <- function(tree, edge_index) {
  node <- tree$edge[edge_index, 2]
  desc <- phytools::getDescendants(tree, node)
  desc <- desc[desc <= length(tree$tip.label)]
  tree$tip.label[desc]
}

### count upstream gain nodes AND total CREs for ONE simmap
count_upstream_gains <- function(simmap_obj, trait_name, matrix_data) {
  
  
  svsp_tips <- grep("SVSP", simmap_obj$tip.label, value = TRUE)
  tree_sub <- drop.tip(simmap_obj, setdiff(simmap_obj$tip.label, svsp_tips))
  
  maps <- tree_sub$maps
  edges <- seq_along(maps)
  
  gain_edges <- edges[sapply(maps, function(edge_map) {
    states <- names(edge_map)
    states[1] == "0" & any(states == "1")
  })]
  
  tip_order <- tree_sub$tip.label
  ann <- matrix_data[match(tip_order, matrix_data$Tip100), ]
  trait_vec <- as.numeric(matrix_data[match(tip_order, matrix_data$Tip100), trait_name] > 0)
  
  upstream_gain_nodes <- c()
  cre_counts <- numeric()  # store CRE counts sequentially
  
  for (edge_index in gain_edges) {
    
    desc_tips <- get_descendant_tips(tree_sub, edge_index)
    desc_ann <- ann[match(desc_tips, ann$Tip100), ]
    desc_trait <- trait_vec[match(desc_tips, tree_sub$tip.label)]
    
    # Check if any descendant is promoter/enhancer with trait==1
    cre_mask <- desc_ann$type %in% c("promoter", "enhancer") & desc_trait == 1
    
    if (any(cre_mask, na.rm = TRUE)) {
      
      upstream_gain_nodes <- c(upstream_gain_nodes, edge_index)
      
      # Count all CREs downstream of this gain (add them up)
      cre_counts <- c(cre_counts, sum(cre_mask, na.rm = TRUE))
    }
  }
  
  list(
    num_upstream_gain_nodes = length(upstream_gain_nodes),
    total_CREs = sum(cre_counts, na.rm = TRUE)
  )
}

### MASTER LOOP, run om all TFBSs sequentially  ###
results <- list()

for (trait in names(simmap_list)) {
  cat("Processing trait:", trait, "\n")
  
  sims <- simmap_list[[trait]]
  n_sims <- length(sims)
  
  res_trait <- tibble(
    trait = trait,
    sim   = seq_len(n_sims),
    num_upstream_gain_nodes = NA_real_,
    total_CREs = NA_real_
  )
  
  for (i in seq_len(n_sims)) {
    simmap_obj <- sims[[i]]
    
    # Count upstream gains and total CREs for this simulation
    result <- count_upstream_gains(simmap_obj, trait_name = trait, matrix_data)
    
    res_trait$num_upstream_gain_nodes[i] <- result$num_upstream_gain_nodes
    res_trait$total_CREs[i] <- result$total_CREs
  }
  
  results[[trait]] <- res_trait
}

results_df <- bind_rows(results)


trait_summary <- results_df %>%
  group_by(trait) %>%
  summarise(
    mean_upstream_gains = mean(num_upstream_gain_nodes, na.rm = TRUE),
    var_upstream_gains  = var(num_upstream_gain_nodes, na.rm = TRUE),
    mean_CREs_affected = mean(total_CREs, na.rm = TRUE),
    var_CREs_affected  = var(total_CREs, na.rm = TRUE),
    .groups = "drop"
  )

trait_summary <- trait_summary %>%
  arrange(desc(mean_upstream_gains)) %>%
  mutate(trait = factor(trait, levels = trait))

p1 <- ggplot(trait_summary, aes(x = mean_upstream_gains, y = trait)) + 
  geom_point(aes(size = mean_CREs_affected), color = "steelblue", alpha = 0.7) +
  scale_y_discrete(limits = rev(levels(trait_summary$trait))) +
  theme_minimal() +
  ylab("TFBS") +
  xlab("Mean Upstream Gains") +
p1

# save(p1, file = "p1_plot.RData")
