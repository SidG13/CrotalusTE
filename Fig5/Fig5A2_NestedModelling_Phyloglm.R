library(tidyverse)
library(ape)
library(phylolm)

set.seed(1234)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig5')

## Test Set A phyloglm ##
tree <- read.tree('../z_data/All_Cat1_and_Cv1_hAT-Tip100.nwk')

SetA <- read.table('../z_data/FeatureSetA_Prewired_Tip100_Features.txt', header = T)

SetA <- SetA %>% 
  column_to_rownames('Tip100')

SetA <- SetA %>%
  mutate(type_bin = as.numeric(factor(type)) - 1) %>% 
  mutate(type_bin = 1 - type_bin)

SetA$species <- if_else(grepl('Cv1', rownames(SetA)), 'Cv1', 'Cat1')

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, rownames(SetA)))
tree_pruned$edge.length[tree_pruned$edge.length == 0] <- 0.01 ## https://github.com/lamho86/phylolm/issues/46

model1_A <- phyloglm(type_bin ~ PRE_SUM + any_pioneer_bound, phy = tree_pruned, data = SetA, method = "logistic_IG10") # ignore warnings
summary(model1_A)

## Test Set A+B phyloglm ##
SetB <- read.table('../z_data/FeatureSetB_Rewired_Tip100_Features.txt', header = T)

SetB <- SetB %>% 
  column_to_rownames('Tip100')

SetB <- SetB %>%
  mutate(type_bin = as.numeric(factor(type)) - 1) %>% 
  mutate(type_bin = 1 - type_bin)

SetAB <- SetA %>%
  rownames_to_column(var = "Tip100") %>%  # make rownames a column
  left_join(
    SetB %>%
      rownames_to_column(var = "Tip100") %>%
      select(-type_bin, -type),
    by = "Tip100"
  ) %>%
  column_to_rownames(var = "Tip100")

model2_AB <- phyloglm(type_bin ~ PRE_SUM + any_pioneer_bound + RE_SUM + any_tissuespecific_bound + co_binding_count, phy = tree_pruned, data = SetAB, method = "logistic_IG10")
summary(model2_AB)

## Test Set A+B+C phyloglm ##
SetC <- read.table('../z_data/FeatureSetC_Unwired_Tip100_Features.txt', header = T)

SetABC <- SetAB %>%
  rownames_to_column(var = "Tip100") %>%
  left_join(SetC %>% select(-type), by = "Tip100") %>%
  column_to_rownames(var = "Tip100")

model3_ABC <- phyloglm(type_bin ~ PRE_SUM + any_pioneer_bound + RE_SUM + any_tissuespecific_bound + co_binding_count + pre_wired_losses + hub_lost, phy = tree_pruned, data = SetABC, method = "logistic_IG10")
summary(model3_ABC) # cobinding counts and pre-wired motif losses are the key features
exp(coef(model3_ABC)) # odds ratios


#### Evaluations ####
AIC(model1_A)
AIC(model2_AB)
AIC(model3_ABC)

# Log-likelihoods and number of parameters
# Number of regression coefficients
df1 <- length(coef(model1_A))
df2 <- length(coef(model2_AB))
df3 <- length(coef(model3_ABC))

# Numeric log-likelihoods
ll1_num <- logLik(model1_A)$logLik
ll2_num <- logLik(model2_AB)$logLik
ll3_num <- logLik(model3_ABC)$logLik

# Likelihood ratio tests
LR12 <- 2 * (ll2_num - ll1_num)
df12 <- df2 - df1
p12 <- pchisq(LR12, df12, lower.tail = FALSE)

LR23 <- 2 * (ll3_num - ll2_num)
df23 <- df3 - df2
p23 <- pchisq(LR23, df23, lower.tail = FALSE)

data.frame(
  comparison = c("model1 vs model2", "model2 vs model3"),
  LR_stat = c(LR12, LR23),
  df = c(df12, df23),
  p_value = c(p12, p23)
)


aic_values <- data.frame(
  model = c("P", "PR", "PRU"),
  AIC = c(AIC(model1_A), AIC(model2_AB), AIC(model3_ABC))
)
aic_values$deltaAIC <- aic_values$AIC - min(aic_values$AIC)

p1 <- ggplot(aic_values, aes(x = model, y = deltaAIC, fill = model)) +
  geom_col() +
  geom_text(aes(label = round(deltaAIC, 2)), vjust = -0.5) +
  ylab("deltaAIC (relative to best)") +
  theme_minimal() +
  ggtitle("Nested Phyloglm Model Comparison")

lrt_df <- data.frame(
  comparison = c("P vs PR", "PR vs PRU"),
  p_value = c(p12, p23)
)
lrt_df$negLogP <- -log10(lrt_df$p_value)

p2 <- ggplot(lrt_df, aes(x = comparison, y = negLogP, fill = comparison)) +
  geom_col() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  ggtitle("Likelihood Ratio Test for Nested Phyloglm Models")

cowplot::plot_grid(p1,p2)