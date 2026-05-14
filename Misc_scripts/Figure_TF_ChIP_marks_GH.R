library(tidyverse)
library(rtracklayer)
library(RColorBrewer)
library(chromVAR)
library(ggcoverage)
library(scales)
library(cowplot)
library(gggenes)

# setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Figures/Figure_TF_ChIP')

cvv_tip100 <- read.table('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/Ultimate_hAT-Tip100_name_conversion.txt', header = T)
cvv_tip100 <- cvv_tip100 %>% 
  filter(species == 'croAtr2') %>% 
  select(chrom, start, end, short_name)

#### Calculate average functional signals around Tip100s ####
## Start EHF ##
track.folder = '/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Figures/Figure_TF_ChIP/TF_ChIP_bams'

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder), value = T)))
sample.meta <- sample.meta %>% 
  filter(SampleName == '01_0L7A_Q-02W8Arlington_CV130814_EHF_other_i06-27_T1_1.srt.nodup') %>% 
  mutate(Type = 'EHF')

window_size <- 1000  # ±2kb

all_signals <- lapply(1:nrow(cvv_tip100), function(i) {
  element_name <- cvv_tip100[i, 'short_name']
  chrom <- cvv_tip100[i, 'chrom']
  center <- round((cvv_tip100[i, 'start'] + cvv_tip100[i, 'end']) / 2)
  start_pos <- center - window_size
  end_pos <- center + window_size
  
  track.df <- tryCatch({
    LoadTrackFile(track.folder = track.folder, format = "bam", 
                  bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                  meta.info = sample.meta,
                  single.nuc = TRUE,
                  single.nuc.region = paste0(chrom, ":", start_pos, "-", end_pos),
                  norm.method = 'RPGC')
  }, error = function(e) {
    message("Skipping problematic file for element ", element_name, ": ", e$message)
    return(NULL) 
  })
  
  if (is.null(track.df)) return(NULL) 
  
  track_avg <- track.df %>%
    mutate(rel_position = start - center) %>%
    select(rel_position, score) %>%
    mutate(name = element_name)
  
  return(track_avg)
})

all_signals <- Filter(Negate(is.null), all_signals)

# Try a different grouping method, where we compare Cat1-hAT-Tip100s that are bound according to ATAC FPs 
binding_mat <- read.table('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/TOBIAS_out_files/2SPECIES_Joint_Tip100_TE_TFBS_Ciiider_scan_footprint_bound_trinMat_May2025.txt', header = T)
TFBS_translation <- read.table('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/TFBS_homology_best_translation_table.txt', header = T)
EHF_sites <- TFBS_translation %>% 
  filter(Transcription.Factor.Name ==  'EHF') %>% 
  pull(TFBS_homology_group)

binding_mat_EHF <- binding_mat[, EHF_sites] %>% 
  filter(grepl('Cat1', rownames(.)))

EHF_bound_TEs <- rownames(binding_mat_EHF)[rowSums(binding_mat_EHF == 2) > 0] # Need to have at least one row/column be 2, which is present-bound
EHF_notBoundNoMotif_TEs <- rownames(binding_mat_EHF)[rowSums(binding_mat_EHF != 0) == 0]
EHF_notBoundHasMotif_TEs <- setdiff(rownames(binding_mat_EHF), c(EHF_bound_TEs, EHF_notBoundNoMotif_TEs))

nrow(binding_mat_EHF) == length(c(EHF_bound_TEs, EHF_notBoundNoMotif_TEs, EHF_notBoundHasMotif_TEs)) # CHECK that each TE has been classified

#signal_combined2 <- bind_rows(all_signals)
#signal_combined2 <- signal_combined2 %>% 
#  mutate(Type = if_else(name %in% EHF_bound_TEs, 'EHF_bound', 'EHF_not_bound'))

signal_combined3 <- bind_rows(all_signals)
signal_combined3 <- signal_combined3 %>% 
  mutate(Type = case_when(name %in% EHF_bound_TEs ~ 'EHF_bound_TEs',
                          name %in% EHF_notBoundNoMotif_TEs ~ 'EHF_notBoundNoMotif_TEs',
                          name %in% EHF_notBoundHasMotif_TEs ~ 'EHF_notBoundHasMotif_TEs'))

#score_final2 <- signal_combined2 %>%
#  group_by(rel_position, Type) %>%
#  summarise(score = mean(score, na.rm = TRUE)) %>%
#  ungroup()

score_final3 <- signal_combined3 %>%
  group_by(rel_position, Type) %>%
  summarise(score = mean(score, na.rm = TRUE)) %>%
  filter(!is.na(Type)) %>% 
  ungroup()

# ggplot(score_final2, aes(x = rel_position, y = score, color = Type)) +
#   geom_line(data = filter(score_final2, Type == "EHF_bound"), size = 1, alpha = 0.3) +
#   geom_line(data = filter(score_final2, Type == "EHF_not_bound"), size = 1, alpha = 0.3) +
#   theme_minimal() +
#   scale_color_manual(values = c('EHF_not_bound' = brewer.pal(8, 'Dark2')[3], 'EHF_bound' = brewer.pal(8, 'Dark2')[5])) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "position relative to element center", y = "average EHF ChIP-seq signal") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   geom_smooth(data = score_final2, aes(x = rel_position, y = score, color = Type), 
#               method = "loess", size = 1, alpha = 0.5, linetype = "solid")

p_EHF <- ggplot(score_final3, aes(x = rel_position, y = score, color = Type)) +
  geom_line(data = filter(score_final3, Type == "EHF_bound_TEs"), size = 1, alpha = 0.3) +
  geom_line(data = filter(score_final3, Type == "EHF_notBoundNoMotif_TEs"), size = 1, alpha = 0.3) +
  geom_line(data = filter(score_final3, Type == "EHF_notBoundHasMotif_TEs"), size = 1, alpha = 0.3) +
  theme_minimal() +
  scale_color_manual(values = c('EHF_bound_TEs' = brewer.pal(8, 'Dark2')[3], 'EHF_notBoundNoMotif_TEs' = brewer.pal(8, 'Dark2')[5], 'EHF_notBoundHasMotif_TEs' = brewer.pal(8, 'Dark2')[6])) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "position relative to element center", y = "average EHF ChIP-seq signal") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_smooth(data = score_final3, aes(x = rel_position, y = score, color = Type), 
              method = "loess", size = 1, alpha = 0.5, linetype = "solid")

## End EHF ##