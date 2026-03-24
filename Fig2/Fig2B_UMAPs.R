library(Seurat)
library(tidyverse)
library(Signac)
library(SCENIC)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

setwd('/Volumes/SeagatePortableDrive/Crotalus_TE_analysis/z_Github_Submission_Feb2026/Fig2')

#### Renaming process ####
Cvv_3VG_Merged <- readRDS('../z_data/Cvv_3VG_Merged_postClustering_renamedGenes.rds') # Find on Zenodo

# Plot UMAP
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
p1 <- DimPlot(Cvv_3VG_Merged, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, label.box = T, label.color = 'white', repel = TRUE, cols = safe_colorblind_palette) + ggtitle("RNA")
p2 <- DimPlot(Cvv_3VG_Merged, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, label.box = T, label.color = 'white', repel = TRUE, cols = safe_colorblind_palette) + ggtitle("ATAC")
p3 <- DimPlot(Cvv_3VG_Merged, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, label.box = T, label.color = 'white', repel = TRUE, cols = safe_colorblind_palette) + ggtitle("WNN")
p_umap1 <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p1 <- DimPlot(Cvv_3VG_Merged, reduction = "umap.rna", group.by = "orig.ident", label = TRUE, label.size = 2.5, label.box = T, label.color = 'white', repel = TRUE, cols = brewer.pal(n = 3, name = "Dark2")) + ggtitle("RNA")
p2 <- DimPlot(Cvv_3VG_Merged, reduction = "umap.atac", group.by = "orig.ident", label = TRUE, label.size = 2.5, label.box = T, label.color = 'white', repel = TRUE, cols = brewer.pal(n = 3, name = "Dark2")) + ggtitle("ATAC")
p3 <- DimPlot(Cvv_3VG_Merged, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, label.box = T, label.color = 'white', repel = TRUE, cols = brewer.pal(n = 3, name = "Dark2")) + ggtitle("WNN")
p_umap2 <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

plot_grid(p_umap1, p_umap2, nrow = 2)

#### Use UMAP to visualize ####
gene_AUC <- read.table('../z_data/eRegulon_AUC_geneBased.txt', header = T) # Find on Zenodo
region_AUC <- read.table('../z_data/eRegulon_AUC_regionBased.txt', header = T) # Find on Zenodo
gene_AUC <- gene_AUC %>% 
  column_to_rownames('Cell')
rownames(gene_AUC) <- sub("___.*", "", rownames(gene_AUC))
region_AUC <- region_AUC %>% 
  column_to_rownames('Cell')
rownames(region_AUC) <- sub("___.*", "", rownames(region_AUC))

all(rownames(gene_AUC) %in% colnames(Cvv_3VG_Merged))

Cvv_3VG_Merged <- AddMetaData(object = Cvv_3VG_Merged, metadata = gene_AUC)
Cvv_3VG_Merged@meta.data[is.na(Cvv_3VG_Merged@meta.data)] <- 0
Cvv_3VG_Merged <- AddMetaData(object = Cvv_3VG_Merged, metadata = region_AUC)
Cvv_3VG_Merged@meta.data[is.na(Cvv_3VG_Merged@meta.data)] <- 0

p1 <- FeaturePlot(Cvv_3VG_Merged, features = "SVSP7", order = T, reduction = 'wnn.umap', min.cutoff = "q10", cols = c('ghostwhite', 'dodgerblue2'))
p2 <- FeaturePlot(Cvv_3VG_Merged, features = "GRHL1_._._.77g.", order = T, reduction = 'wnn.umap', min.cutoff = "q75", cols = c('ghostwhite', 'firebrick4'))
p3 <- FeaturePlot(Cvv_3VG_Merged, features = "SVSP10", order = T, reduction = 'wnn.umap', min.cutoff = "q10", cols = c('ghostwhite', 'dodgerblue2'))
p4 <- FeaturePlot(Cvv_3VG_Merged, features = "CREB3L2_._._.36g.", order = T, reduction = 'wnn.umap', min.cutoff = "q75", cols = c('ghostwhite', 'firebrick4'))
p5 <- FeaturePlot(Cvv_3VG_Merged, features = "SVSP8", order = T, reduction = 'wnn.umap', cols = c('ghostwhite', 'dodgerblue2'))
p6 <- FeaturePlot(Cvv_3VG_Merged, features = "EHF_._._.136g.", order = T, reduction = 'wnn.umap', min.cutoff = "q75", cols = c('ghostwhite', 'firebrick4'))

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)