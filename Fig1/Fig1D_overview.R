# Adapted from scripts used in Perry et al. 2022

library(tidyverse)
library(rtracklayer)
library(RColorBrewer)
library(ggcoverage)
library(scales)
library(cowplot)
library(gggenes)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig1')


#### START Gene and Hi-C tracks for C viridis and C atrox ####

### START VIRIDIS ###
GEX <- read.table('../z_data/All_croVir_VG_featurecounts_Jan2025_VST_counts.txt')
GEX <- GEX %>% 
  rownames_to_column('gene') %>% 
  filter(str_detect(gene, 'SVSP')) %>% 
  select(c(gene, CV1297, CV1298, CV1299)) %>% # 3 samples corresponding to hi-depth ATAC-seq libraries
  mutate(meanVST = rowMeans(across(c(CV1297, CV1298, CV1299)))) %>% 
  mutate(gene = gsub('Venom_', '', gene)) %>% 
  select(gene, meanVST)

gtf <- rtracklayer::import('../z_data/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf') %>% # find on Zenodo file
  as.data.frame()

SVSP.info <- gtf %>% 
  filter(str_detect(gene_id, 'SVSP')) %>% 
  filter(str_detect(gene_id, 'fgenesh')) %>% 
  filter(type == 'gene') %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(Name = gsub('_', '', Name)) %>% 
  left_join(GEX, by = c('Name' = 'gene'))


###
### ******* This is where you set the x-axis limits. (These are also the limits of the above .paf alignment file)
###
SVSP.reg.start <- 8495000 # SVSP
SVSP.reg.end <- 8986585 # SVSP

SVSP.reg.length <- paste(c(round((SVSP.reg.end-SVSP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

Tip100_TEs <- read.table('../z_data/Ultimate_hAT-Tip100_name_conversion.txt', header = T)

Tip100_TEs <- Tip100_TEs %>% 
  filter(species == 'croVir3') %>% 
  select(id, chrom, start, end) %>% 
  filter(grepl('SVSP', id)) %>% 
  mutate(type = case_when(
    grepl("_E", id) ~ "E",
    grepl("_P", id) ~ "P",
    grepl("_O", id) ~ "O",
    TRUE ~ "Unknown")) %>% 
  mutate(TE_colour = case_when(
    grepl("_E", id) ~ brewer.pal(8, 'Dark2')[1],
    grepl("_P", id) ~ brewer.pal(8, 'Dark2')[2],
    grepl("_O", id) ~ brewer.pal(8, 'Dark2')[6],
    TRUE ~ "black")) %>% 
  mutate(index = row_number()) %>% 
  mutate(site = (start + end) / 2)

#### Plotting Gene Track ####

p_genes <- ggplot(SVSP.info, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = meanVST)) +
  geom_segment(aes(x=SVSP.reg.start, xend=SVSP.reg.end, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  #geom_segment(inherit.aes = F,data=Tip100_TEs %>% filter(type == 'E'),aes(x=SVSP.reg.start, xend=SVSP.reg.end,y=' Tip100_Enhancer',yend=' Tip100_Enhancer'),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=Tip100_TEs %>% filter(type == 'E'), aes(x=site, y=' Tip100_Enhancer'),size=4, color = brewer.pal(8, 'Dark2')[1], alpha = 0.4) +
  #geom_segment(inherit.aes = F,data=Tip100_TEs %>% filter(type == 'P'),aes(x=SVSP.reg.start,xend=SVSP.reg.end,y='  Tip100_Promoter',yend='  Tip100_Promoter'),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=Tip100_TEs %>% filter(type == 'P'), aes(x=site, y='  Tip100_Promoter'),size=4, color = brewer.pal(8, 'Dark2')[2], alpha = 0.4) +
  #geom_segment(inherit.aes = F,data=Tip100_TEs %>% filter(type == 'O'),aes(x=SVSP.reg.start,xend=SVSP.reg.end,y='   Tip100_Other',yend='   Tip100_Other'),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=Tip100_TEs %>% filter(type == 'O'), aes(x=site, y='   Tip100_Other'),size=4, color = brewer.pal(8, 'Dark2')[6], alpha = 0.4) +
  #geom_text(inherit.aes = F, data = Tip100_TEs %>% filter(type == 'P'), aes(label = index, x = site, y = '  Tip100_Promoter', color = 'red')) +
  #geom_text(inherit.aes = F, data = Tip100_TEs %>% filter(type == 'O'), aes(label = index, x = site, y = '   Tip100_Other', color = 'red')) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(limits=c(SVSP.reg.start, SVSP.reg.end),expand=c(0,0),oob = scales::oob_keep) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank())

# Load data and plot
track.folder = '../z_data/hidepth_scaled_bigwigs' # Folder contains {CV1297,CV1298,CV1299}*scaleFactorNorm.bw files, find on Zenodo

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bw', '' , grep('bw$', list.files(track.folder), value = T)))
sample.meta <- sample.meta %>% 
  mutate(Type = str_extract(SampleName, "^[^-]+"))

# load regions of interest from bam files
chrom = "scaffold-mi2"
start_pos = SVSP.reg.start
end_pos = SVSP.reg.end
track.df = LoadTrackFile(track.folder = track.folder, format = "bw", # change to bam
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = F, # change to T for BAM
                         region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change region to single.nuc.region for BAM,
)

# Average 3 samples in each start-end bin
score_average <- track.df %>%
  group_by(seqnames, start, end) %>%
  summarise(score = mean(score)) %>% 
  ungroup()
score_average$Type <- 'Average'

## ## Plotting Hi-C curves ## ##
hic_loops <- read.table('../z_data/Viridis_HiC_loops_fromBlair.bedpe', header = T)
hic_loops <- hic_loops %>% 
  filter(chr1 == 'scaffold-mi2' & chr2 == 'scaffold-mi2') %>% 
  #mutate(Type = 'Average') %>% # to plot overtop the average ATAC-seq signal
  rowwise() %>% 
  mutate(
    # Overlap1: Check for overlap between x1/x2 and Tip100_TEs$site
    overlap1 = Tip100_TEs$type[which(Tip100_TEs$site >= x1 & Tip100_TEs$site <= x2)[1]],
    index1 = Tip100_TEs$index[which(Tip100_TEs$site >= x1 & Tip100_TEs$site <= x2)[1]],
    
    # Overlap2: Check for overlap between y1/y2 and Tip100_TEs$site
    overlap2 = Tip100_TEs$type[which(Tip100_TEs$site >= y1 & Tip100_TEs$site <= y2)[1]],
    index2 = Tip100_TEs$index[which(Tip100_TEs$site >= y1 & Tip100_TEs$site <= y2)[1]]
  ) %>% 
  ungroup() %>% 
  #mutate(loop_alpha = ifelse((!is.na(overlap1) | !is.na(overlap2)), 1, 0.5)) %>% 
  #mutate(loop_strength = case_when(
  #  !is.na(overlap1) & !is.na(overlap2) ~ 1,  # Both overlap1 and overlap2 are not NA
  #  !is.na(overlap1) | !is.na(overlap2) ~ 0.6,  # Either overlap1 or overlap2 is not NA
  #  TRUE ~ 0.2                               # All other cases
  #)) %>% 
  mutate(loop_colour = case_when(
    overlap1 == 'E' | overlap2 == 'E' ~ brewer.pal(8, 'Dark2')[4],
    overlap2 == 'P' | overlap2 == 'P' ~ brewer.pal(8, 'Dark2')[4],
    overlap1 == 'O' | overlap2 == 'O' ~ brewer.pal(8, 'Dark2')[4],
    TRUE ~ 'grey60'
  )) %>% 
  mutate(centroid1_correct = (x1 + x2)/2,
         centroid2_correct = (y1 + y2)/2)

# SVSP enhancers and promoters
enh_pr <- read.table('../z_data/Cviridis_SVSP_cres.txt', header = T)
enh_pr$cre_shape <- as.character(enh_pr$cre_shape)

p_atac_array <- 
  ggplot() +
  theme_minimal() +
  geom_col(data = score_average, aes(x = start, y = score, fill = Type)) +
  #annotate(geom = "rect",
  #         xmin = peaks$start,
  #         xmax = peaks$end,
  #         ymin = -Inf,
  #         ymax = +Inf,
  #         alpha = 0.2) +
  geom_curve(data = hic_loops, aes(x = centroid1_correct, y = 0, xend = centroid2_correct, yend = 0), 
             curvature = 0.15, color = hic_loops$loop_colour) + # make a curve showing a Hi-C contact
  geom_rect(data = Tip100_TEs, aes(xmin = start, xmax = end, ymin = 2100, ymax=2200), color = Tip100_TEs$TE_colour, fill = Tip100_TEs$TE_colour, alpha = 0.6) +
  geom_point(data = enh_pr, aes(x = (start + end)/2, y = 2200, pch = cre_shape, fill = 'red'), color = 'red', size = 3, alpha = 0.4) +
  facet_wrap(~factor(Type, levels = c(sample.meta$Type, 'Average')), ncol = 1, strip.position = "right", scales = 'free_y') +
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        #strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) +
  guides(fill = 'none') +
  scale_fill_manual(values = c("Average" = "blue")) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq score")

p_viridis_genes <- plot_grid(p_genes, p_atac_array, nrow = 2, axis = 'lr', align = 'hv', rel_heights = c(1,0.8))

### END VIRIDIS ###

### START ATROX ###
GEX <- read.table('../z_data/All_croAtr2_DESeq2_VSTcounts_Nov2024.txt')
GEX <- GEX %>% 
  rownames_to_column('gene') %>% 
  filter(str_detect(gene, 'SVSP')) %>% 
  select(c(gene, CV1308, CV1311, CV1312)) %>% # 3 samples corresponding to hi-depth ATAC-seq libraries
  mutate(meanVST = rowMeans(across(c(CV1308, CV1311, CV1312)))) %>% 
  mutate(gene = gsub('Venom_', '', gene)) %>% 
  select(gene, meanVST) %>% 
  filter(gene != 'SVSP-29') # remove SVSP-29

gtf <- rtracklayer::import('../z_data/croAtr2_ragtag_scaffold_REHEADER_withToxins_noDashes_10.09.24_MT_manualedits_namesMod.sorted.gtf') %>% # Find on Zenodo
  as.data.frame()

SVSP.info <- gtf %>% 
  filter(str_detect(gene_id, 'SVSP')) %>% 
  filter(type == 'transcript') %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold','', seqnames)) %>% 
  mutate(gene_id = gsub('Venom_', '', gene_id)) %>% 
  filter(gene_id != 'SVSP-29') %>% 
  left_join(GEX, by = c('gene_id' = 'gene'))


###
### ******* This is where you set the x-axis limits. (These are also the limits of the above .paf alignment file)
###
SVSP.reg.start <- 11846033 # SVSP
SVSP.reg.end <- 12982468 # SVSP

SVSP.reg.length <- paste(c(round((SVSP.reg.end-SVSP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

Tip100_TEs <- read.table('../z_data/Ultimate_hAT-Tip100_name_conversion.txt', header = T)

Tip100_TEs <- Tip100_TEs %>% 
  filter(species == 'croAtr2') %>% 
  select(id, chrom, start, end) %>% 
  filter(grepl('SVSP', id)) %>% 
  mutate(type = case_when(
    grepl("_E", id) ~ "E",
    grepl("_P", id) ~ "P",
    grepl("_O", id) ~ "O",
    TRUE ~ "Unknown")) %>% 
  mutate(TE_colour = case_when(
    grepl("_E", id) ~ brewer.pal(8, 'Dark2')[1],
    grepl("_P", id) ~ brewer.pal(8, 'Dark2')[2],
    grepl("_O", id) ~ brewer.pal(8, 'Dark2')[6],
    TRUE ~ "black")) %>% 
  mutate(index = row_number()) %>% 
  mutate(site = (start + end) / 2)

#### Plotting Gene Track ####

p_genes <- ggplot(SVSP.info, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = meanVST)) +
  geom_segment(aes(x=SVSP.reg.start, xend=SVSP.reg.end, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  #geom_segment(inherit.aes = F,data=Tip100_TEs %>% filter(type == 'E'),aes(x=SVSP.reg.start, xend=SVSP.reg.end,y=' Tip100_Enhancer',yend=' Tip100_Enhancer'),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=Tip100_TEs %>% filter(type == 'E'), aes(x=site, y=' Tip100_Enhancer'),size=4, color = brewer.pal(8, 'Dark2')[1], alpha = 0.4) +
  #geom_segment(inherit.aes = F,data=Tip100_TEs %>% filter(type == 'P'),aes(x=SVSP.reg.start,xend=SVSP.reg.end,y='  Tip100_Promoter',yend='  Tip100_Promoter'),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=Tip100_TEs %>% filter(type == 'P'), aes(x=site, y='  Tip100_Promoter'),size=4, color = brewer.pal(8, 'Dark2')[2], alpha = 0.4) +
  #geom_segment(inherit.aes = F,data=Tip100_TEs %>% filter(type == 'O'),aes(x=SVSP.reg.start,xend=SVSP.reg.end,y='   Tip100_Other',yend='   Tip100_Other'),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=Tip100_TEs %>% filter(type == 'O'), aes(x=site, y='   Tip100_Other'),size=4, color = brewer.pal(8, 'Dark2')[6], alpha = 0.4) +
  #geom_text(inherit.aes = F, data = Tip100_TEs %>% filter(type == 'P'), aes(label = index, x = site, y = '  Tip100_Promoter', color = 'red')) +
  #geom_text(inherit.aes = F, data = Tip100_TEs %>% filter(type == 'O'), aes(label = index, x = site, y = '   Tip100_Other', color = 'red')) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(limits=c(SVSP.reg.start, SVSP.reg.end),expand=c(0,0),oob = scales::oob_keep) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank())

# Load data and plot
track.folder = '../z_data/atrox_scaled_bws' # Folder contains {CV1308,CV1311,CV1312}*scaleFactorNorm.bw files, find on Zenodo
# Use CV1308, CV1311 and CV1312

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bw', '' , grep('bw$', list.files(track.folder), value = T)))
sample.meta <- sample.meta %>% 
  mutate(Type = str_extract(SampleName, "^[^-]+"))

# load regions of interest from bam files
chrom = "scaffoldmi1"
start_pos = SVSP.reg.start
end_pos = SVSP.reg.end
track.df = LoadTrackFile(track.folder = track.folder, format = "bw", # change to bam
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = F, # change to T for BAM
                         region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change region to single.nuc.region for BAM,
)

score_average <- track.df %>%
  group_by(seqnames, start, end) %>%
  summarise(score = mean(score)) %>% 
  ungroup()
score_average$Type <- 'Average'

## Plotting Hi-C curves ##
hic_loops <- read.table('../z_data/Catrox_enriched_pixels_merged_loops.bedpe', header = T)
hic_loops <- hic_loops %>% 
  filter(chr1 == 'scaffoldmi1' & chr2 == 'scaffoldmi1') %>% 
  #mutate(Type = 'Average') %>% # to plot overtop the average ATAC-seq signal
  rowwise() %>% 
  mutate(
    # Overlap1: Check for overlap between x1/x2 and Tip100_TEs$site
    overlap1 = Tip100_TEs$type[which(Tip100_TEs$site >= x1 & Tip100_TEs$site <= x2)[1]],
    index1 = Tip100_TEs$index[which(Tip100_TEs$site >= x1 & Tip100_TEs$site <= x2)[1]],
    
    # Overlap2: Check for overlap between y1/y2 and Tip100_TEs$site
    overlap2 = Tip100_TEs$type[which(Tip100_TEs$site >= y1 & Tip100_TEs$site <= y2)[1]],
    index2 = Tip100_TEs$index[which(Tip100_TEs$site >= y1 & Tip100_TEs$site <= y2)[1]]
  ) %>% 
  ungroup() %>% 
  #mutate(loop_alpha = ifelse((!is.na(overlap1) | !is.na(overlap2)), 1, 0.5)) %>% 
  #mutate(loop_strength = case_when(
  #  !is.na(overlap1) & !is.na(overlap2) ~ 1,  # Both overlap1 and overlap2 are not NA
  #  !is.na(overlap1) | !is.na(overlap2) ~ 0.6,  # Either overlap1 or overlap2 is not NA
  #  TRUE ~ 0.2                               # All other cases
  #)) %>% 
  mutate(loop_colour = case_when(
    overlap1 == 'E' | overlap2 == 'E' ~ brewer.pal(8, 'Dark2')[4],
    overlap2 == 'P' | overlap2 == 'P' ~ brewer.pal(8, 'Dark2')[4],
    overlap1 == 'O' | overlap2 == 'O' ~ brewer.pal(8, 'Dark2')[4],
    TRUE ~ 'grey60'
  )) %>% 
  mutate(centroid1_correct = (x1 + x2)/2,
         centroid2_correct = (y1 + y2)/2)

# SVS Enhancers and promoters (Catrox)
enh_pr <- read.table('../z_data/Catrox_SVSP_cres.txt', header = T)
enh_pr$start <- as.numeric(enh_pr$start)
enh_pr$end <- as.numeric(enh_pr$end)
enh_pr$cre_shape <- as.character(enh_pr$cre_shape)

p_atac_array <- 
  ggplot() +
  theme_minimal() +
  geom_col(data = score_average, aes(x = start, y = score, fill = Type)) +
  #annotate(geom = "rect",
  #         xmin = peaks$start,
  #         xmax = peaks$end,
  #         ymin = -Inf,
  #         ymax = +Inf,
  #         alpha = 0.2) +
  geom_curve(data = hic_loops, aes(x = centroid1_correct, y = 0, xend = centroid2_correct, yend = 0), 
             curvature = 0.15, color = hic_loops$loop_colour) + # make a curve showing a Hi-C contact
  geom_rect(data = Tip100_TEs, aes(xmin = start, xmax = end, ymin = 210, ymax=220), color = Tip100_TEs$TE_colour, fill = Tip100_TEs$TE_colour, alpha = 0.6) +
  geom_point(data = enh_pr, aes(x = (start + end)/2, y = 220, pch = cre_shape, fill = 'red'), color = 'red', size = 3, alpha = 0.4) +
  facet_wrap(~factor(Type, levels = c(sample.meta$Type, 'Average')), ncol = 1, strip.position = "right", scales = 'free_y') +
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        #strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) +
  guides(fill = 'none') +
  scale_fill_manual(values = c("Average" = "blue")) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq score")

p_atrox_genes <- plot_grid(p_genes, p_atac_array, nrow = 2, axis = 'lr', align = 'hv', rel_heights = c(1,0.8))

### plotting ###
plot_grid(p_viridis_genes, p_atrox_genes, nrow = 2, axis = 'lr', align = 'hv')

