## Try to quantitatively model a TF having multiple TFBSs, and count appropriately ##
##

library(tidyverse)
library(ape)
library(phytools)
library(ggtree)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig4')

tree <- read.tree("../z_data/All_Cat1_and_Cv1_hAT-Tip100.nwk")
matrix_data <- read.table('../z_data/NEW_CRE_hATs_prewired_trinary_matrix_for_modelling_QUANT_TFBS_Oct2025.txt', header = T)

#### Writing a loop to run simmap for the re-wired TFBSs: RUN ON EXTERNAL SERVER AT UTA ####
trait_cols <- colnames(matrix_data)[5:59] # all pre-wired TFBSs
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

#### End Loop ####
# Save all simmaps as a single RDS file
#saveRDS(simmap_list, file = "simmaps_55prewired_TFBSs.rds")


#### Plotting individual TFBSs ####
simmap_list <- readRDS('../z_data/Z_simmaps_13prewired_TFs_QUANTITATIVE_for_unwiring_x100_reps.rds') # Find on Zenodo

# sim_num <- sample(1:100, 1) # random number for the simulation
# trait_name <- "Ddit3Cebpa"
# 
# simmap_obj <- simmap_list[[trait_name]][[sim_num]]
# 
# svsp_tips <- grep("SVSP", simmap_obj$tip.label, value = TRUE)
# tree_sub_simmap <- drop.tip(simmap_obj, setdiff(simmap_obj$tip.label, svsp_tips))
# 
# tree_label <- paste0(trait_name, " sim ", sim_num)
# 
# get_edge_transitions <- function(edge_map) {
#   s <- as.numeric(names(edge_map))
#   if (length(s) < 2) return(NULL)
#   data.frame(
#     from = s[-length(s)],
#     to   = s[-1],
#     stringsAsFactors = FALSE
#   )
# }
# 
# loss_edges <- which(
#   sapply(tree_sub_simmap$maps, function(m) {
#     tr <- get_edge_transitions(m)
#     !is.null(tr) && any(tr$to < tr$from)
#   })
# )
# 
# gain_edges <- which(
#   sapply(tree_sub_simmap$maps, function(m) {
#     tr <- get_edge_transitions(m)
#     !is.null(tr) && any(tr$to > tr$from)
#   })
# )
# 
# states <- sort(unique(as.numeric(names(unlist(tree_sub_simmap$maps)))))
# 
# state_cols <- setNames(
#   viridis::viridis(length(states)),
#   as.character(states)
# )
# 
# plotSimmap(tree_sub_simmap,
#            colors = state_cols,
#            fsize = 0.8,
#            lwd = 2)
# title(main = tree_label)
# 
# edgelabels(pch = 21, bg = "red",   cex = 1.4, edge = loss_edges)
# ann <- matrix_data[match(tree_sub_simmap$tip.label, matrix_data$Tip100), ]
# 
# tiplabels(
#   pch = 25, bg = "blue",   cex = 1.2,
#   tip = which(ann$type == "promoter")
# )
# 
# tiplabels(
#   pch = 25, bg = "orange", cex = 1.2,
#   tip = which(ann$type == "enhancer")
# )
# 
# tip_states <- matrix_data[match(tree_sub_simmap$tip.label,
#                                 matrix_data$Tip100), trait_name]
# tip_states <- as.numeric(tip_states)
# 
# x_tip <- max(nodeHeights(tree_sub_simmap)) + 0.05
# y_tip <- 1:length(tree_sub_simmap$tip.label)
# 
# box_width  <- 0.03
# box_height <- 0.5
# 
# for (i in seq_along(tip_states)) {
#   state_i <- tip_states[i]
#   rect(
#     xleft   = x_tip,
#     ybottom = y_tip[i] - box_height / 2,
#     xright  = x_tip + box_width,
#     ytop    = y_tip[i] + box_height / 2,
#     col     = state_cols[as.character(state_i)],
#     border  = "black"
#   )
# }
# 
# legend("topleft",
#        legend = names(state_cols),
#        fill   = state_cols,
#        border = "black",
#        bty = "n",
#        title = "Trait state",
#        cex = 0.7)


#### Run a loop or something to get all num_upstream_loss nodes per sim per tfbs ####

get_descendant_tips <- function(tree, edge_index) {
  node <- tree$edge[edge_index, 2]
  desc <- phytools::getDescendants(tree, node)
  tips <- desc[desc <= length(tree$tip.label)]
  tree$tip.label[tips]
}

get_loss_amount <- function(edge_map) {
  st <- as.numeric(names(edge_map))
  diffs <- diff(st)
  sum(abs(diffs[diffs < 0]))
}

# Count upstream loss nodes and total loss magnitude for ONE simmap, do individually in loop
count_upstream_losses <- function(simmap_obj, matrix_data) {
  
  svsp_tips <- grep("SVSP", simmap_obj$tip.label, value = TRUE)
  tree_sub <- drop.tip(simmap_obj, setdiff(simmap_obj$tip.label, svsp_tips))
  
  maps <- tree_sub$maps
  edges <- seq_along(maps)
  
  loss_edges <- edges[sapply(maps, function(edge_map) {
    st <- as.numeric(names(edge_map))
    any(diff(st) < 0)
  })]
  
  ann <- matrix_data[match(tree_sub$tip.label, matrix_data$Tip100), ]
  
  upstream_loss_nodes <- c()
  loss_amounts <- c()
  
  for (edge_index in loss_edges) {
    desc_tips <- get_descendant_tips(tree_sub, edge_index)
    desc_ann <- ann[match(desc_tips, ann$Tip100), ]
    
    if (any(desc_ann$type %in% c("promoter","enhancer"), na.rm=TRUE)) {
      upstream_loss_nodes <- c(upstream_loss_nodes, edge_index)
      loss_amounts <- c(loss_amounts, get_loss_amount(tree_sub$maps[[edge_index]]))
    }
  }
  
  list(
    num_upstream_loss_nodes = length(upstream_loss_nodes),
    total_loss_magnitude    = sum(loss_amounts)
  )
}

results <- list()

for (trait in names(simmap_list)) {
  cat("Processing trait:", trait, "\n")
  
  sims <- simmap_list[[trait]]
  n_sims <- length(sims)
  
  res_trait <- tibble(
    trait = trait,
    sim   = seq_len(n_sims),
    num_upstream_loss_nodes = NA_real_,
    total_loss_magnitude    = NA_real_
  )
  
  for (i in seq_len(n_sims)) {
    simmap_obj <- sims[[i]]
    res <- count_upstream_losses(simmap_obj, matrix_data)
    
    res_trait$num_upstream_loss_nodes[i] <- res$num_upstream_loss_nodes
    res_trait$total_loss_magnitude[i]    <- res$total_loss_magnitude
  }
  
  results[[trait]] <- res_trait
}

results_df <- bind_rows(results)
results_df


trait_summary <- results_df %>%
  group_by(trait) %>%
  summarise(
    mean_upstream_losses = mean(num_upstream_loss_nodes, na.rm = TRUE),
    var_upstream_losses  = var(num_upstream_loss_nodes, na.rm = TRUE),
    mean_loss_magnitude = mean(total_loss_magnitude, na.rm = TRUE),
    var_loss_magnitude  = var(total_loss_magnitude, na.rm = TRUE),
    .groups = "drop"
  )

motif_depletion_df <- data.frame( # Taken from Supplement
                     stringsAsFactors = FALSE,
                          check.names = FALSE,
                     `TFBS` = c("ATF4","JUN","Ddit3Cebpa.1",
                                          "TBX3.1","Arid3a.1","Dlx4","ZBTB26",
                                          "NFATC1.3","Arid3a","NFATC1.1",
                                          "Ddit3Cebpa","Dlx4.2","Ddit3Cebpa.2",
                                          "NFATC1.2","NFIX","ZBTB26.1",
                                          "MEIS1","Arid3a.2","SOX9.1","Dlx4.1",
                                          "Arid3a.4","Arid3a.3","SOX10.1",
                                          "NFATC1","TBX3","SOX10","SOX9",
                                          "Arid3a.5","Arid3a.6","TFAP4"),
              Motif.frequency.in.CREs = c(0.76,0.68,0.63,0.68,0.76,
                                          0.68,0.8,0.59,0.93,0.59,0.63,0.56,
                                          0.61,0.85,0.93,0.76,0.9,0.85,
                                          0.8,0.88,0.9,0.9,0.71,0.66,0.73,
                                          0.66,0.71,0.76,0.78,0.78),
        `Motif.frequency.in.non-CREs` = c(0.23,0.24,0.26,0.31,0.39,
                                          0.35,0.5,0.29,0.64,0.3,0.36,0.3,
                                          0.35,0.6,0.68,0.51,0.68,0.64,
                                          0.61,0.73,0.75,0.76,0.59,0.56,0.64,
                                          0.58,0.63,0.69,0.71,0.76),
  `Frequency.difference` = c(0.53,0.45,0.37,0.37,0.37,
                                          0.33,0.3,0.3,0.29,0.29,0.27,0.26,
                                          0.26,0.25,0.25,0.24,0.23,0.22,
                                          0.19,0.15,0.15,0.14,0.12,0.1,0.09,
                                          0.08,0.08,0.07,0.07,0.02)
)

motif_depletion_df <- motif_depletion_df %>% 
  mutate(TF = str_replace(TFBS, "\\..*$", "")) %>% 
  group_by(TF) %>% 
  summarize(mean_depletion = mean(Frequency.difference))

trait_summary <- trait_summary %>% 
  left_join(motif_depletion_df, by = c('trait'='TF')) %>% 
  arrange(mean_upstream_losses) %>%
  mutate(trait = factor(trait, levels = trait))
mean(trait_summary$mean_upstream_losses)


ggplot(trait_summary, aes(x = mean_upstream_losses,
                          y = trait,
                          size = mean_loss_magnitude)) +
  geom_point(aes(color = mean_depletion)) +
  scale_color_viridis_c(option = 'magma') +
  scale_size_continuous(range = c(3, 12)) + 
  labs(
    x = "Mean upstream losses",
    y = "TFBS",
    size = "Mean loss magnitude",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title  = element_text(hjust = 0.5),
    legend.position = "right"
  )
