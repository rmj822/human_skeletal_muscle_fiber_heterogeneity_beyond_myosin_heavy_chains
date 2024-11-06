################################################################################################################################################
#################################################       USED LINKS      ########################################################################
################################################################################################################################################

# https://github.com/hbctraining/scRNA-seq/tree/master/schedule
# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/09_merged_SC_marker_identification.md


################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(cowplot)
library(ggpubr)


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/1 Clustering/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata")

# Assign identity of resolution 0.7 to clusters ---------------------------------------------
Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.7"

# Load in gene annotation file --------------------------------------------
annotations <- read.csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

# Set to SCT assay for visualization, and log-transform data ------------------------------
DefaultAssay(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "SCT"


# GOALS
    # To determine the gene markers for each of the clusters
    # To identify cell types of each cluster using markers
    # To determine whether there's a need to re-cluster based on cell type markers, perhaps clusters need to be merged or split

# OPTIONS
 # 1: Identification of all markers for each cluster: this analysis compares each cluster against all others and outputs the genes that are differentially expressed/present.
        # Advantage: Useful for identifying unknown clusters and improving confidence in hypothesized cell types.
 # 2: Identification of conserved markers for each cluster: genes that are differentially expressed/present within each condition first, and then reports those genes that are conserved in the cluster across all conditions. These genes can help to figure out the identity for the cluster.
        # Advantage: Useful with more than one condition to identify cell type markers that are conserved across conditions.
 # 3: Marker identification between specific clusters: this analysis explores differentially expressed genes between specific clusters.
        # Advantage: Useful for determining differences in gene expression between clusters that appear to be representing the same celltype (i.e with markers that are similar) from the above analyses.


# Options:
  # logfc.threshold = minimum log2 fold change for avg expression in cluster compared to avg of all other clusters (default 0.25)
      # CON 1: could miss those cell markers that are expressed in a small fraction of cells within the cluster of interest, if the average log2FC doesn't meet the threshold
      # CON 2: could return a lot of metabolic/ribosomal genes due to slight differences in metabolic output by different cell types, which are not as useful to distinguish cell type identities
  # min.diff.pct = min % diff in % of cell expressing gene in cluster compared to % expression in all other clusters
      # CON 1: could miss those cell markers that are expressed in all cells, but are highly up-regulated in this specific cell type
  # min.pct = only test genes that are present in X % of cells in cluster (default 0.1)
      # CON 1: if set to a very high value could incur many false negatives due to the fact that not all genes are detected in all cells (even if it is expressed)
  # only.pos: only detect genes that are higher in cluster compared to other clusters

################################################################################################################################################
#########################################      IDENTIFICATION OF SLOW VS FAST TYPE  ###################################################
################################################################################################################################################


# FeaturePlot
FeaturePlot(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
            features = c("MYH7","MYH2", "TNNT1", "TNNT3", "TPM3", "TPM1"),
            order = TRUE,
            min.cutoff = 'q10',
            reduction = "umap",
            label = F,
            repel = TRUE)

`################################################################################################################################################
##################################      IDENTIFICATION OF MARKERS BETWEEN ALL CLUSTERS     ##################################################
################################################################################################################################################


# Compare cluster 0 and 1 (both fast clusters) ----------------------------
markers_all <- FindAllMarkers(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Determine top10 genes per cluster
top10 <- markers_all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- markers_all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Make heatmap for top 10 genes
Seurat::DoHeatmap(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, features = top5$gene, raster=F) +
  scale_fill_gradientn(colors = c("#134057", "white", "#BC4749"))

    # Interpretation:
          # Cluster 0, 1 and 3 are fast fibers and do seem different
          # Cluster 6 has both slow and fast fibers and seem really different than other clusters
          # Cluster 5 seems different from other slow clusters
          # Cluster 2 and 4 different enough???

################################################################################################################################################
##################################      IDENTIFICATION OF MARKERS BETWEEN CLUSTER 2 AND 4     ##################################################
################################################################################################################################################

# Compare cluster 0 and 1 (both fast clusters) ----------------------------
markers_2vs4 <- FindMarkers(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                            ident.1 = 2,
                            ident.2 = 4,
                            logfc.threshold = 0.25)

# Add gene annotations ----------------------------
markers_2vs4 <- markers_2vs4 %>%
  rownames_to_column(var = "gene") %>%
  left_join(annotations, by = c("gene" = "GENENAME"))

# Make column with pct diff
markers_2vs4$pct.diff <- markers_2vs4$pct.1 - markers_2vs4$pct.2

# Determine top50 different genes
markers_2vs4_top50 <- markers_2vs4 %>% dplyr::arrange(desc(abs(avg_log2FC))) %>% slice(1:50)

# Make heatmap for top 50 genes
Seurat::DoHeatmap(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, features = markers_2vs4_top50$gene, raster=F) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749"))

# Only small differences --> combine into one clusters


################################################################################################################################################
#######################################################     REASSIGN CLUSTERS     ##############################################################
################################################################################################################################################

# Assign identity of resolution 0.7 to clusters ---------------------------------------------
Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.7"

# Reclustering
filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest <- RenameIdents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                                                   "0" = "Fast1",
                                                                   "1" = "Fast2",
                                                                   "2" = "Slow1",
                                                                   "3" = "Fast3",
                                                                   "4" = "Slow1",
                                                                   "5" = "Slow2",
                                                                   "6" = "Intermediate1")



# Check UMAP to confirm correct cluster assignment
DimPlot(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = TRUE)
