setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig1')

library(GenomicRanges)
library(tidyverse)
library(RColorBrewer)

## Viridis
viridis_cres <- read.table('../z_data/Cviridis_SVSP_cres.txt', header = T)
viridis_tips <- read.table('../z_data/Ultimate_hAT-Tip100_name_conversion.txt', header = T)
viridis_tips <- viridis_tips %>%
  filter(species == 'croVir3') %>% 
  select(chrom, start, end)

viridis_cres_gr <- GRanges(seqnames = viridis_cres$seqnames, 
                           ranges = IRanges(start = viridis_cres$start, end = viridis_cres$end),
                           type = viridis_cres$type)
viridis_tips_gr <- GRanges(seqnames = viridis_tips$chrom, 
                           ranges = IRanges(start = viridis_tips$start, end = viridis_tips$end))
overlaps <- findOverlaps(viridis_cres_gr, viridis_tips_gr)

overlap_df_cvv <- cbind(viridis_cres[queryHits(overlaps),], viridis_tips[subjectHits(overlaps),])
colnames(overlap_df_cvv) <- c("seqnames", "cre_start", "cre_end", "cre_name", "cre_type", "cre_shape", "chr", "tip_start", "tip_end")
d1 <- as.data.frame(table(overlap_df_cvv$cre_type)/table(viridis_cres$type))
d1$species <- 'croVir3'
d1[2,2] <- 0.9 # in reality, only 1 promoter, SVSP11's does NOT overlap a Tip100, so 9/10 promoters do


# Atrox
atrox_cres <- read.table('../z_data/Catrox_SVSP_cres.txt', header = T)
atrox_tips <- read.table('../z_data/Ultimate_hAT-Tip100_name_conversion.txt', header = T)
atrox_tips <- atrox_tips %>%
  filter(species == 'croAtr2') %>% 
  select(chrom, start, end)

atrox_cres_gr <- GRanges(seqnames = atrox_cres$seqnames, 
                         ranges = IRanges(start = atrox_cres$start, end = atrox_cres$end),
                         type = atrox_cres$type)
atrox_tips_gr <- GRanges(seqnames = atrox_tips$chrom, 
                         ranges = IRanges(start = atrox_tips$start, end = atrox_tips$end))
overlaps <- findOverlaps(atrox_cres_gr, atrox_tips_gr)

overlap_df_catr <- cbind(atrox_cres[queryHits(overlaps),], atrox_tips[subjectHits(overlaps),])
colnames(overlap_df_catr) <- c("seqnames", "cre_start", "cre_end", "cre_name", "cre_type", "cre_shape", "chr", "tip_start", "tip_end")
d2 <- as.data.frame(table(overlap_df_catr$cre_type)/table(atrox_cres$type))
d2$species <- 'croAtr2'

df_combined <- rbind(d1,d2)

ggplot(df_combined, aes(x = Var1, y = Freq, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("croVir3" = brewer.pal(8, "Dark2")[1], "croAtr2" = brewer.pal(8, "Dark2")[6])) +  
  labs(
    title = "fraction of SVSP CREs overlapping (Cv1/Cat1)-hAT-Tip100s"
  ) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
