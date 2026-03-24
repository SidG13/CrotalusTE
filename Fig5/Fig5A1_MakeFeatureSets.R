library(tidyverse)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig5')

main_df <- read.table('../z_data/NEW_SVSP_hATs_trinary_matrix_for_Kaas_modelling_Dec2025.txt', header = T)

#### Make Pre-wired Features ####
## Feature 1: States count of pre-wired sites
main_df <- main_df %>%
  select(-ends_with("_RE")) %>%
  mutate(
    PRE_SUM = rowSums(select(., ends_with("_PRE")), na.rm = TRUE)
  )

main_df <- main_df %>%
  mutate(type_bin = as.numeric(factor(type)) - 1) %>% 
  mutate(type_bin = 1 - type_bin)

cor.test(main_df$PRE_SUM, main_df$type_bin, method = "spearman") # higher pre-wired count is associated with CRE, good sign

## Feature 2: At least 1 pioneer bound, all pioneers are pre-wired
pioneers <- grep('SOX9|SOX10|Cebpa|MEIS1', colnames(main_df), value = T)

main_df <- main_df %>%
  mutate(
    any_pioneer_bound = as.integer(
      rowSums(
        select(
          .,
          matches(paste0("^(", paste(pioneers, collapse = "|"), ")"))
        ) == 2,
        na.rm = TRUE
      ) > 0
    )
  )

cor.test(main_df$any_pioneer_bound, main_df$type_bin, method = "pearson") # Tips with pioneer binding present are more likely to be CRE

pre_df <- main_df %>% 
  select(Tip100, type, PRE_SUM, any_pioneer_bound)
# write.table(pre_df, 'FeatureSetA_Prewired_Tip100_Features.txt', col.names = T, row.names = F, quote = F, sep = '\t')

#### Make Re-wired Features ####
rm(list=ls())

## Feature 1: States count of re-wired sites
main_df <- read.table('../z_data/NEW_SVSP_hATs_trinary_matrix_for_Kaas_modelling_Dec2025.txt', header = T)

main_df <- main_df %>%
  select(-ends_with("_PRE")) %>%
  mutate(
    RE_SUM = rowSums(select(., ends_with("_RE")), na.rm = TRUE)
  )

main_df <- main_df %>%
  mutate(type_bin = as.numeric(factor(type)) - 1) %>% 
  mutate(type_bin = 1 - type_bin)

cor.test(main_df$RE_SUM, main_df$type_bin, method = "spearman") # higher re-wired count is associated with CRE, good sign

## Feature 2: At least 1 tissue-specific bound
tissuespecific <- grep('EHF|ELF5|Creb3l2|CREB3L1', colnames(main_df), value = T)

main_df <- main_df %>%
  mutate(
    any_tissuespecific_bound = as.integer(
      rowSums(
        select(
          .,
          matches(paste0("^(", paste(tissuespecific, collapse = "|"), ")"))
        ) == 2,
        na.rm = TRUE
      ) > 0
    )
  )

cor.test(main_df$any_tissuespecific_bound, main_df$type_bin, method = "pearson") # Tips with tissuespecific binding present are more likely to be CRE

re_df1 <- main_df %>% 
  select(Tip100, type, RE_SUM, any_tissuespecific_bound)

## Feature 3: Number of significant pre-re cobinding interactions ##
main_df <- read.table('../z_data/NEW_SVSP_hATs_trinary_matrix_for_Kaas_modelling_Dec2025.txt', header = T)
cobinding_df <- data.frame(
  stringsAsFactors = FALSE,
       check.names = FALSE,
                       `Pre_wired` = c("SOX10.1","SOX10.1","SOX10.1","SOX9",
                                       "SOX9","SOX9","SOX10.1","SOX10.1",
                                       "SOX10.1","SOX10.1","SOX10.1","SOX10.1",
                                       "SOX10.1","SOX10.1","SOX10.1","SOX9",
                                       "SOX10.1","SOX10.1","SOX10.1",
                                       "SOX10.1","SOX10.1","Arid3a","Arid3a",
                                       "Arid3a","Arid3a","SOX10","SOX9","Arid3a",
                                       "Arid3a","Arid3a","Arid3a","Arid3a",
                                       "NFIX","JUN","JUN","JUN","Arid3a",
                                       "Arid3a","ATF4","ATF4","ATF4","SOX9",
                                       "Arid3a.2","Arid3a.2","Arid3a.2","SOX9",
                                       "TBX3.1","TBX3.1","TBX3.1","NFIX",
                                       "NFIX","NFIX","Arid3a","SOX10",
                                       "Arid3a.1","TBX3.1","TBX3.1","ATF4","JUN",
                                       "Arid3a.2","SOX10","Arid3a","SOX9",
                                       "TBX3.1","RORC","NFIX","NFIX",
                                       "Ddit3Cebpa.1","Ddit3Cebpa.1","Ddit3Cebpa.1",
                                       "Arid3a","Arid3a.2","Arid3a.1","TBX3.1",
                                       "NFIX","Arid3a","JUN"),
                        `Re_wired` = c("NFATC1.1","EHF","ELF5","NFATC1.1",
                                       "EHF","ELF5","Creb3l2","TBX3",
                                       "NFATC1","TBX3.3","MEIS1.2","Irx2",
                                       "Irx2.1","TBX3.1","MEIS1.1","TBX3","CREB3L1",
                                       "MEIS1.3","POU6F1","TBX3.2",
                                       "Creb3l2.1","MEIS1.1","NFATC1.1","EHF","ELF5",
                                       "NFATC1","CREB3L1","MEIS1.3","Irx2",
                                       "Irx2.1","TBX3.1","CREB3L1","MEIS1.1",
                                       "NFATC1.1","EHF","ELF5","Creb3l2",
                                       "TBX3","NFATC1.1","EHF","ELF5",
                                       "MEIS1.2","NFATC1.1","EHF","ELF5","TBX3.3",
                                       "NFATC1.1","EHF","ELF5","NFATC1.1",
                                       "EHF","ELF5","NFATC1","MEIS1.2",
                                       "Creb3l2","CREB3L1","NFATC1","TBX3",
                                       "Creb3l2","Creb3l2","TBX3.3","MEIS1.2",
                                       "POU6F1","Creb3l2","Creb3l2","Irx2",
                                       "Irx2.1","NFATC1.1","EHF","ELF5","TBX3.3",
                                       "TBX3","TBX3","TBX3","CREB3L1",
                                       "POU6F1","MEIS1.2")
                ) %>% 
  mutate(Pre_wired = paste0(Pre_wired, '_PRE')) %>% 
  mutate(Re_wired = paste0(Re_wired, '_RE'))

main_df$co_binding_count <- 0

for (i in 1:nrow(main_df)) {
  
  count <- 0

  for (j in 1:nrow(cobinding_df)) {
    pre_col <- cobinding_df$Pre_wired[j]
    re_col  <- cobinding_df$Re_wired[j]
    
    if (pre_col %in% colnames(main_df) && re_col %in% colnames(main_df)) {
      
      if (main_df[i, pre_col] == 2 && main_df[i, re_col] == 2) {
        count <- count + 1
      }
    }
  }

  main_df$co_binding_count[i] <- count
}

# View result
head(main_df)

main_df <- main_df %>%
  mutate(type_bin = as.numeric(factor(type)) - 1) %>% 
  mutate(type_bin = 1 - type_bin)

cor.test(main_df$co_binding_count, main_df$type_bin, method = "pearson") # Tips with more cobinders are more likely to be CRE

re_df2 <- main_df %>% 
  select(Tip100, co_binding_count)

re_df <- re_df1 %>% 
  left_join(re_df2, by = 'Tip100')

# write.table(re_df, 'FeatureSetB_Rewired_Tip100_Features.txt', col.names = T, row.names = F, quote = F, sep = '\t')

#### Make Un-wired Features ####
rm(list=ls())

## Feature 1: Count number of pre-wired motif losses (0s) ##
main_df <- read.table('../z_data/NEW_SVSP_hATs_trinary_matrix_for_Kaas_modelling_Dec2025.txt', header = T)

pre_cols <- grep("_PRE$", colnames(main_df), value = TRUE)

# Count number of 0s row-wise and store in a new column
main_df <- main_df %>%
  rowwise() %>%
  mutate(pre_wired_losses = sum(c_across(all_of(pre_cols)) == 0)) %>%
  ungroup()

main_df <- main_df %>%
  mutate(type_bin = as.numeric(factor(type)) - 1) %>% 
  mutate(type_bin = 1 - type_bin)

cor.test(main_df$pre_wired_losses, main_df$type_bin, method = "pearson") # Tips with more unwired TFBSs are less likely to be CRE

## Feature 2: Binary lost hubs ##
hub_cols <- grep("ATF4|JUN|Ddit3Cebpa", colnames(main_df), value = TRUE) # double check

# Count the number of 0s row-wise for these columns
main_df <- main_df %>%
  rowwise() %>%
  mutate(hub_lost = if_else(any(c_across(all_of(hub_cols)) == 0), 0, 1)) %>%
  ungroup()

cor.test(main_df$hub_lost, main_df$type_bin, method = "pearson") # Tips with more unwired hubs are less likely to be CRE

unwired_df <- main_df %>% 
  select(Tip100, type, pre_wired_losses, hub_lost)

# write.table(unwired_df, 'FeatureSetC_Unwired_Tip100_Features.txt', col.names = T, row.names = F, quote = F, sep = '\t')

