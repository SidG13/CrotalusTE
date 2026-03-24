library(tidyverse)
library(purrr)
library(broom)
library(cowplot)
library(ape)
library(phytools)
library(phylolm)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig4')

#### Begin phylo modelling (rewired) ####
matrix_data_re <- read.table('../z_data/NEW_CRE_hATs_rewired_trinary_matrix_for_modelling_Oct2025.txt', header = T)
tree <- read.tree("../z_data/All_Cat1_and_Cv1_hAT-Tip100.nwk")
trait_cols <- colnames(matrix_data_re)[5:ncol(matrix_data_re)]

phylo_matrix_data <- matrix_data_re %>% 
  mutate(species = if_else(grepl('Cv1', Tip100), 'croVir3', 'croAtr2')) %>% 
  relocate(species, .after = 'type') %>% 
  column_to_rownames('Tip100') %>% 
  filter(type %in% c('promoter', 'enhancer')) %>% 
  mutate(across(.cols = -(1:4),  # Recode to bound and unbound (0 = unbound, 1/2 = bound)
                .fns = ~ case_when(
                  . %in% c(0, 1) ~ 0,
                  . == 2        ~ 1
                )))

fit <- phylolm( # warning is fine... ignore
  formula = log(tpm) ~ MEIS1.2 + species,
  data = phylo_matrix_data,
  phy = tree,
  model = "BM"
)

# summary(fit)$coefficients['MEIS1.1', "p.value"]
# summary(fit)

phylolm_results <- data.frame(
  trait = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (TFBS in trait_cols) {
  # Fit the model
  fit <- phylolm(
    formula = as.formula(paste("log(tpm) ~", TFBS, '+ species')),
    data = phylo_matrix_data,
    phy = tree,
    model = "BM"
  )
  
  # Get p-value
  p <- summary(fit)$coefficients[as.character(TFBS), "p.value"]
  
  phylolm_results <- rbind(phylolm_results, data.frame(trait = TFBS, p_value = p))
}

phylolm_results <- phylolm_results %>% 
  mutate(signif = if_else(p_value < 0.05, T, F)) %>% 
  arrange(p_value)
phylolm_results

phylolm_results$FDR_BH <- p.adjust(phylolm_results$p_value, method = "BH")
phylolm_results$signif_FDR <- phylolm_results$FDR_BH < 0.05
phylolm_results

# Now check all individual models for re-wired tfbss
boxwhisker_plots <- list()

for (tfbs in subset(phylolm_results, signif == TRUE) %>% pull(trait)) {
  
  p <- ggplot(phylo_matrix_data, aes(
    x = factor(.data[[tfbs]], levels = c(0,1), labels = c("unbound","bound")),
    y = log(tpm)
  )) +
    geom_boxplot(fill = "skyblue", alpha = 0.5, outlier.shape = F) +
    geom_point(aes(shape = type, color = species), alpha = 0.6, size = 3.5, position = position_jitter(width = 0.15)) +
    scale_color_manual(values = c("croVir3" = "#189B77", "croAtr2" = "#E2AC26")) +
    labs(
      x = NULL,
      y = "log(TPM)",
      title = tfbs
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5)
    )
  
  boxwhisker_plots[[tfbs]] <- p
}

plot_grid(plotlist = boxwhisker_plots)

rm(list=ls())
dev.off()

#### Begin NON-phylo modelling (rewired) ####
matrix_data_re <- read.table('../z_data/NEW_CRE_hATs_rewired_trinary_matrix_for_modelling_Oct2025.txt', header = T)
tree <- read.tree("../z_data/All_Cat1_and_Cv1_hAT-Tip100.nwk")
trait_cols <- colnames(matrix_data_re)[5:ncol(matrix_data_re)]

nonphylo_matrix_data <- matrix_data_re %>% 
  mutate(species = if_else(grepl('Cv1', Tip100), 'croVir3', 'croAtr2')) %>% 
  relocate(species, .after = 'type') %>% 
  column_to_rownames('Tip100') %>% 
  filter(type %in% c('promoter', 'enhancer')) %>% 
  mutate(across(.cols = -(1:4),  # Recode to bound and unbound (0 = unbound, 1/2 = bound)
                .fns = ~ case_when(
                  . %in% c(0,1) ~ 0,
                  . == 2        ~ 1
                )))

#fit <- lm(formula = log(tpm) ~ EHF + species, data = nonphylo_matrix_data) # NOT phylo-corrected, just lm
#summary(fit)
#summary(fit)$coefficients[2, "Pr(>|t|)"]

ggplot(nonphylo_matrix_data, aes(x = factor(EHF, levels = c(0,1), labels = c("unbound","bound")),
                              y = log(tpm))) +
  geom_boxplot(fill = "skyblue", alpha = 0.5) +
  geom_point(aes(color = species, shape = type), alpha = 0.6, size = 3.5) +
  scale_color_manual(values = c("croVir3" = "#189B77", "croAtr2" = "#E2AC26")) +
  labs(
    x = NULL,
    y = "log(TPM)"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5))

nonphylolm_results <- data.frame(
  trait = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (TFBS in trait_cols) {
  fit <- lm(
    formula = as.formula(paste("log(tpm) ~", TFBS)),
    data = nonphylo_matrix_data
  )
  
  # Get p-value
  p <- summary(fit)$coefficients[2, "Pr(>|t|)"]
  
  nonphylolm_results <- rbind(nonphylolm_results, data.frame(trait = TFBS, p_value = p))
}

nonphylolm_results <- nonphylolm_results %>% 
  mutate(signif = if_else(p_value < 0.05, T, F)) %>% 
  arrange(p_value)
nonphylolm_results

# Now check all individual models for re-wired tfbss
boxwhisker_plots <- list()

for (tfbs in subset(nonphylolm_results, signif == TRUE) %>% pull(trait)) {
  
  p <- ggplot(nonphylo_matrix_data, aes(
    x = factor(.data[[tfbs]], levels = c(0,1), labels = c("unbound","bound")),
    y = log(tpm)
  )) +
    geom_boxplot(fill = "skyblue", alpha = 0.5, outlier.shape = NA) +
    geom_point(aes(shape = type, color = species), alpha = 0.6, size = 3.5, position = position_jitter(width = 0.15)) +
    scale_color_manual(values = c("croVir3" = "#189B77", "croAtr2" = "#E2AC26")) +
    labs(
      x = NULL,
      y = "log(TPM)",
      title = tfbs
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5)
    )
  
  boxwhisker_plots[[tfbs]] <- p
}

plot_grid(plotlist = boxwhisker_plots)