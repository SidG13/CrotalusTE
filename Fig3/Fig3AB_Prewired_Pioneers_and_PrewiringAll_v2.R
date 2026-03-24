###    ###
# This script is used to make a simplified figure to show percentages of pre-wiring, based on BWP's schematic #
###    ###

library(tidyverse)
library(ggtree)
library(ape)
library(RColorBrewer)
library(ggnewscale)
library(purrr)
library(GenomicRanges)
library(phylolm)

#### Fig 3A ####
trin_mat <- read.table('../z_data/2SPECIES_Joint_Tip100_TE_TFBS_Ciiider_scan_footprint_bound_trinMat_May2025.txt', header = T)
pre_sites <- data.frame(
  stringsAsFactors = FALSE,
  TFBS.Name = c("Arid3a","MEIS1",
                "NFIX","Arid3a.1","SOX9","Arid3a.2","SOX10","NFATC1",
                "ATF4","Arid3a.3","Ddit3Cebpa","Arid3a.4",
                "SOX10.1","NFATC1.1","ZBTB26","NFATC1.2","Arid3a.5",
                "Arid3a.6","JUN","TBX3","Dlx4","Dlx4.1","TBX3.1",
                "Ddit3Cebpa.1","TFAP4","Ddit3Cebpa.2","Dlx4.2",
                "ZBTB26.1","SOX9.1","NFATC1.3"),
  TFBS.Alignment.Start = c(369L,348L,375L,156L,
                           46L,421L,35L,41L,312L,620L,130L,598L,384L,
                           234L,491L,537L,606L,613L,311L,480L,325L,601L,
                           116L,315L,582L,113L,319L,538L,547L,171L),
  TFBS.Alignment.End = c(376L,366L,387L,161L,
                         59L,436L,45L,50L,335L,637L,165L,606L,419L,
                         243L,518L,553L,614L,620L,333L,511L,344L,610L,
                         128L,336L,601L,127L,336L,575L,575L,184L),
  TFBS.Homology.Group = c("L934","L878","L959",
                          "L138","L1194","L1070","L882","L1032","L729",
                          "L1704","L116","L1610","L994","L422","L1303",
                          "L1446","L1636","L1669","L720","L1285","L816",
                          "L1620","L47","L755","L1564","L30","L785","L1456",
                          "L1478","L211")
) %>% 
  arrange(TFBS.Alignment.Start)

trin_mat_pre <- trin_mat[, pre_sites$TFBS.Homology.Group]

trin_mat_pre2 <- trin_mat_pre %>%
  rownames_to_column('Tip100') %>%
  mutate(type = if_else(grepl('SVSP', Tip100), 'SVSP', 'nonSVSP'))

tfbs_map <- setNames(pre_sites$TFBS.Name, pre_sites$TFBS.Homology.Group)

colnames(trin_mat_pre2) <- colnames(trin_mat_pre2) %>%
  map_chr(~ ifelse(.x %in% names(tfbs_map), tfbs_map[.x], .x))

trin_long <- trin_mat_pre2 %>%
  pivot_longer(cols = -c(Tip100, type), names_to = "TFBS", values_to = "Value") %>%
  mutate(Value = factor(Value, levels = c(0,1,2)))

tfbs_order <- pre_sites %>%
  arrange(TFBS.Alignment.Start) %>%
  pull(TFBS.Name)

trin_long <- trin_long %>%
  mutate(TFBS = factor(TFBS, levels = tfbs_order))

trin_frac <- trin_long %>%
  group_by(type, TFBS, Value) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(type, TFBS) %>%
  mutate(fraction = count / sum(count))

p1 <- ggplot(trin_frac, aes(x = TFBS, y = fraction, fill = Value)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ type, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = c("0" = "white", "1" = "lightgoldenrod2", "2" = "firebrick")) +
  labs(x = "TFBS", y = "Fraction of values", fill = "Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#### End Fig 3A ####

#### Fig 3B ####
bedfiles <- list.files("../z_data/peak_beds", pattern = "\\.bed$", full.names = TRUE) # Directory of peak .bed files, find on Zenodo

bed_df <- map_dfr(bedfiles, ~ read.table(.x, header = FALSE) %>%
                    mutate(filename = basename(.x)))

bed_df_f <- bed_df %>% 
  mutate(species = case_when(grepl('CV12|Viridis', filename) ~ 'croVir3',
                             TRUE ~ 'croAtr2')) %>% 
  mutate(tissue = case_when(grepl('_12_', filename) ~ 'pancreas',
                            grepl('_5_', filename) ~ 'skin',
                            grepl('_6_', filename) ~ 'skelmusc',
                            TRUE ~ 'venom')) %>% 
  select(-filename) %>% 
  distinct()

# Cv1 Tip and nonvenom peak overlaps
cv1_tips <- read.table('../z_data/Cv1-hAT-Tip100_coords.bed')

cv1_tips_gr <- GRanges(
  seqnames = cv1_tips$V1,
  ranges = IRanges(start = cv1_tips$V2, end = cv1_tips$V3),
  tip_name = cv1_tips$V4
)

# Convert tip_coords dataframe to GRanges object
cv1_peaks <- bed_df_f %>% filter(species == 'croVir3')
cv1_peaks_gr <- GRanges(
  seqnames = cv1_peaks$V1,
  ranges = IRanges(start = cv1_peaks$V2, end = cv1_peaks$V3),
  tissue = cv1_peaks$tissue
)

# Find overlaps
overlaps <- findOverlaps(cv1_tips_gr, cv1_peaks_gr)

# Extract the overlapping tip_coords
overlap_df <- cbind(cv1_tips[queryHits(overlaps), ], cv1_peaks[subjectHits(overlaps),])
colnames(overlap_df) <- c("chr", "TE_start", "TE_end", "TE_name", "seqnames", "peak_start", "peak_end", "species", "tissue")

cv1_tips <- cv1_tips %>%
  mutate(OpenChrPanc = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'pancreas'],
         OpenChrSkin = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'skin'],
         OpenChrSkMusc = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'skelmusc'],
         OpenChrVG = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'venom'])

# Cat1 Tip and nonvenom peak overlaps
cat1_tips <- read.table('../z_data/Cat1-hAT-Tip100.bed')

cat1_tips_gr <- GRanges(
  seqnames = cat1_tips$V1,
  ranges = IRanges(start = cat1_tips$V2, end = cat1_tips$V3),
  tip_name = cat1_tips$V4
)

# Convert tip_coords dataframe to GRanges object
cat1_peaks <- bed_df_f %>% filter(species == 'croAtr2')
cat1_peaks_gr <- GRanges(
  seqnames = cat1_peaks$V1,
  ranges = IRanges(start = cat1_peaks$V2, end = cat1_peaks$V3),
  tissue = cat1_peaks$tissue
)

# Find overlaps
overlaps <- findOverlaps(cat1_tips_gr, cat1_peaks_gr)

# Extract the overlapping tip_coords
overlap_df <- cbind(cat1_tips[queryHits(overlaps), ], cat1_peaks[subjectHits(overlaps),])
colnames(overlap_df) <- c("chr", "TE_start", "TE_end", "TE_name", "seqnames", "peak_start", "peak_end", "species", "tissue")

cat1_tips <- cat1_tips %>%
  mutate(OpenChrPanc = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'pancreas'],
         OpenChrSkin = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'skin'],
         OpenChrSkMusc = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'skelmusc'],
         OpenChrVG = V4 %in% overlap_df$TE_name[overlap_df$tissue == 'venom'])


# Combine tip datasets
cv1_tips <- cv1_tips %>% 
  mutate(V4 = case_when(grepl('NonVen', V4) ~ gsub('NonVen_Tip100', 'NonVen_Cv1-hAT-Tip100', V4),
                        grepl('SVSP', V4) ~ gsub('SVSP_Tip100', 'SVSP_Cv1-hAT-Tip100', V4))) %>% 
  select(-V1, -V2, -V3)

cat1_tips <- cat1_tips %>% 
  select(-V1, -V2, -V3)

all_tips <- rbind(cv1_tips, cat1_tips)
all_tips <- all_tips %>% 
  column_to_rownames('V4')

df <- all_tips %>%
  select(OpenChrVG) %>%
  rownames_to_column("Tip100") %>%
  mutate(
    CRE_status = ifelse(grepl("_E|_P", Tip100), "CRE", "non-CRE"),
    OpenChrVG = as.logical(OpenChrVG)   # make sure it's TRUE/FALSE
  )

name_conv <- read.table('../z_data/Ultimate_hAT-Tip100_name_conversion.txt', header = T)
rownames(all_tips) <- name_conv$short_name[match(rownames(all_tips), name_conv$id)]

TF_trinMat <- read.table('../z_data/2SPECIES_Joint_Tip100_TE_TFBS_Ciiider_scan_footprint_bound_trinMat_May2025.txt', header = T)
binding_table <- read.table('../z_data/2SPECIES_Joint_Tip100_TE_TFBS_Ciiider_scan_footprint_bound_May2025.txt', header = T)
TFBS_translation_table <- read.table('../z_data/TFBS_homology_best_translation_table.txt', header = T)

# We are interested now in the pioneer pre-wired sites (SOX9 and DDIT3::CEBPA and MEIS1) (check Prewired_TFBSs_loci.txt)
# SOX9 (L1194, L1478) and DDIT3::CEBPA (L116, L755, L30) and MEIS1 (L878) and SOX10 (L994)

pre_pioneers <- pre_sites %>% 
  mutate(pioneer = if_else(grepl('SOX9|SOX10|MEIS1|Ddit3Cebpa', TFBS.Name), 'pioneer', 'non-pioneer'))

pioneer_sites <- pre_pioneers %>% 
  filter(pioneer == 'pioneer') %>% 
  pull(TFBS.Homology.Group)

TF_trinMat <- TF_trinMat[, pioneer_sites]

TF_trinMat$AtLeastOneBoundPrePioneer <- 0
# Update only the pioneer rows
TF_trinMat$AtLeastOneBoundPrePioneer <- ifelse(
  apply(TF_trinMat[, pioneer_sites] == 2, 1, any),  # checks each row if any value == 2
  1,  # assign 1 if TRUE
  0   # assign 0 if FALSE
)

TF_trinMat$NoPioneerBound <- 0
TF_trinMat$NoPioneerBound <- apply(
  TF_trinMat[, pioneer_sites],
  1,
  function(x) as.integer(!any(x == 2))
)

all_tips_TF <- all_tips %>%
  rownames_to_column("tip_id") %>%
  left_join(
    TF_trinMat %>% rownames_to_column("tip_id"),
    by = "tip_id"
  ) %>%
  column_to_rownames("tip_id") %>% 
  mutate(across(everything(), ~replace_na(., 0)))

all_tips_pioneers <- all_tips_TF %>% 
  rownames_to_column('Tip100') %>%
  mutate(type = if_else(grepl('SVSP', Tip100), 'SVSP', 'nonSVSP'))

p2 <- all_tips_pioneers %>% 
  select(OpenChrPanc, OpenChrSkin, OpenChrSkMusc, OpenChrVG, AtLeastOneBoundPrePioneer, type) %>% 
  pivot_longer(
    cols = starts_with("OpenChr"),
    names_to = "Tissue",
    values_to = "Euchrom"
  ) %>% 
  group_by(type, Tissue, Euchrom) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type, Tissue) %>%
  mutate(prop = n / sum(n)) %>% 
  ggplot(.,
         aes(x = Tissue, y = prop, fill = Euchrom)) +
  geom_col(width = 0.7) +
  facet_wrap(~type, ncol = 1) +
  scale_fill_manual(
    values = c("TRUE" = "#1f78b4",   # Euchrom (blue)
               "FALSE" = "#f1c40f"), # non-Euchrom (yellow)
    labels = c("FALSE" = "Closed", "TRUE" = "Euchrom")
  ) +
  labs(
    y = "Fraction",
    x = NULL,
    fill = "Chromatin"
  ) +
  theme_bw()

# cowplot::plot_grid(p1, p2, nrow = 2)
p2
p1
