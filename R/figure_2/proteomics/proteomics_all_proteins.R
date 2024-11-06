################################################################################################################################################
#################################################       USED LINKS      ########################################################################
################################################################################################################################################

# https://github.com/hbctraining/scRNA-seq/tree/master/schedule
# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/07_SC_clustering_cells_SCT.md


################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)


# Load filtered proteomics data ---------------------------------------------
data_proteomics <- read.csv("~/single_fiber_heterogeneity/data/data_pca_proteomics.csv") # 974 fibers for 1685 proteins
data_proteomics <- data_proteomics |> dplyr::rename("Protein" = 1) |> tibble::column_to_rownames("Protein")

# Proteome metadata
metadata_proteomics <- read.csv("~/single_fiber_heterogeneity/data/metadata_proteomics.csv")
metadata_proteomics <- metadata_proteomics %>% rename("fiberID" = X)
rownames(metadata_proteomics) <- metadata_proteomics$fiberID

# Check if order of samples is the same
all(colnames(data_proteomics) %in% rownames(metadata_proteomics))
all(colnames(data_proteomics) == rownames(metadata_proteomics))

################################################################################################################################################
################################################       CREATE SEURAT OBJECT      ###############################################################
################################################################################################################################################

seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                            meta.data = metadata_proteomics)


################################################################################################################################################
########################################################      PCA   ############################################################################
################################################################################################################################################

# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                              selection.method = "vst")

# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)

#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,  features = Seurat::VariableFeatures(object = seurat_proteome))

# Explore PCA plots
PCA_PC1_2 <- DimPlot(seurat_proteome, reduction = "pca", dims = c(1,2), group.by = "fiber_type") # No separation of fiber type by PC1, good separation by PC2
PCA_PC1_3 <- DimPlot(seurat_proteome, reduction = "pca", dims = c(1,3), group.by = "fiber_type") # No separation at all by PC3

ggsave(PCA_PC1_2, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/PCA_proteome_PC1_PC2.png",  width = 128, height = 90, units="mm")
ggsave(PCA_PC1_3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/PCA_proteome_PC1_PC3.png",  width = 128, height = 90, units="mm")


################################################################################################################################################
################################################      Identify significant PCS    ##############################################################
################################################################################################################################################

# Explore heatmap of PCs (try to find PC where heatmap starts to look fuzzy, not so distinct between groups) --------

PCA_heatmaps <- Seurat::DimHeatmap(seurat_proteome,
                   dims = 1:15,
                   cells = 500,  # number of cells with most negative or positive PCA scores to use for plotting
                   balanced = TRUE) # Still decent clustering up until PC 15

ggsave(PCA_heatmaps, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/PCA_heatmaps_proteome.png",  width = 128, height = 90, units="mm")


# Elbow plot: visualizes SD of each PC, search for where SD begins to plateau --------
elbow_plot <- Seurat::ElbowPlot(object = seurat_proteome,
                  ndims = 40) # Plateau only after about 30 PCs

ggsave(elbow_plot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/elbowplot_proteome.png",  width = 128, height = 90, units="mm")


# Visual inspection: 40 PCs

# Quantitative determination elbow plateau (two metrics, choose lowest value)--------------------------------

# Metric 1: New PC only contributes 5%, and all other cumulatively contribute 90%
pct <- seurat_proteome[["pca"]]@stdev / sum(seurat_proteome[["pca"]]@stdev) * 100 # Calculate percent of variation for each PC
cumu <- cumsum(pct) # Calculate cumulative percents with each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 # PC 41

# Metric 2: PC where percent change to next PC is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
co2 # PC 23

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Conclusion: use first 9 PCs to cluster (number less important with new versions of Seurat, biggest difference is time to compute with more PCs)
# Other option: just use first 40 clusters with sctransform normalization, since the more clusters, the more variation is accounted for (only downside: computing time)



################################################################################################################################################
####################################################      CLUSTERING      ######################################################################
################################################################################################################################################


# Graph-based clustering using K-nearest neighbor graph  ---------------------------------

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4))

################################################################################################################################################
################################################      DIMENSIONALITY REDUCTION   ##############################################################
################################################################################################################################################

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6)

# Run T-SNE ----------------------------------------------------------------
seurat_proteome <- RunTSNE(seurat_proteome, dims = 1:6)



################################################################################################################################################
#################################################     TEST DIFFERENT RESOLUTIONS  ##############################################################
################################################################################################################################################


# Resolution 0.2 -------------------------------

UMAP_res0.2 <- Seurat::DimPlot(seurat_proteome,
                reduction = "umap",
                label = FALSE,
                label.size = 6,
                group.by = "RNA_snn_res.0.2") +
    ggtitle("Resolution 0.2")

ggsave(UMAP_res0.2, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_res_0.2.png",  width = 128, height = 90, units="mm")

# Resolution 0.4 -------------------------------

UMAP_res0.4 <- Seurat::DimPlot(seurat_proteome,
                               reduction = "umap",
                               label = FALSE,
                               label.size = 6,
                               group.by = "RNA_snn_res.0.4") +
    ggtitle("Resolution 0.4")

ggsave(UMAP_res0.4, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_res_0.4.png",  width = 128, height = 90, units="mm")

# Resolution 0.6 -------------------------------

UMAP_res0.6 <- Seurat::DimPlot(seurat_proteome,
                               reduction = "umap",
                               label = FALSE,
                               label.size = 6,
                               group.by = "RNA_snn_res.0.6") +
    ggtitle("Resolution 0.6")

ggsave(UMAP_res0.6, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_res_0.6.png",  width = 128, height = 90, units="mm")

TSNE_res0.6 <- Seurat::DimPlot(seurat_proteome,
                               reduction = "tsne",
                               label = FALSE,
                               label.size = 6,
                               group.by = "RNA_snn_res.0.6") +
    ggtitle("Resolution 0.6")

ggsave(TSNE_res0.6, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/TSNE_res_0.6.png",  width = 128, height = 90, units="mm")

# Resolution 0.8 -------------------------------

UMAP_res0.8 <- Seurat::DimPlot(seurat_proteome,
                               reduction = "umap",
                               label = FALSE,
                               label.size = 6,
                               group.by = "RNA_snn_res.0.8") +
    ggtitle("Resolution 0.8")

ggsave(UMAP_res0.8, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_res_0.8.png",  width = 128, height = 90, units="mm")

TSNE_res0.8 <- Seurat::DimPlot(seurat_proteome,
                               reduction = "tsne",
                               label = FALSE,
                               label.size = 6,
                               group.by = "RNA_snn_res.0.8") +
    ggtitle("Resolution 0.8")

ggsave(TSNE_res0.8, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/TSNE_res_0.8.png",  width = 128, height = 90, units="mm")

# Resolution 1.0 -------------------------------

UMAP_res1.0 <- Seurat::DimPlot(seurat_proteome,
                               reduction = "umap",
                               label = FALSE,
                               label.size = 6,
                               group.by = "RNA_snn_res.1") +
    ggtitle("Resolution 1.0")

ggsave(UMAP_res1.0, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_res_1.0.png",  width = 128, height = 90, units="mm")


##### START WITH RESOLUTION 0.8


################################################################################################################################################
#########################################      CLUSTERING QUALITY CONTROL PER SUBJECT    #######################################################
################################################################################################################################################

UMAP_subject <- Seurat::DimPlot(seurat_proteome,
                label = TRUE,
                reduction = "umap",
                group.by = "subject")

ggsave(UMAP_subject, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_subject.png",  width = 128, height = 90, units="mm")

################################################################################################################################################
#########################################      CLUSTERING QUALITY CONTROL PER MYH TYPE    #####################################################
################################################################################################################################################

UMAP_fibertype <- Seurat::DimPlot(seurat_proteome,
                                label = TRUE,
                                reduction = "umap",
                                group.by = "fiber_type")

ggsave(UMAP_fibertype, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/UMAP_fibertype.png",  width = 128, height = 90, units="mm")


# Cluster 0: Fast
# Cluster 1: Slow
# Cluster 2: Fast
# Cluster 3: Fast
# Cluster 4: Slow
# Cluster 5: Slow
# Cluster 6: Fast
# Cluster 7: Fast
# Cluster 8: Fast


################################################################################################################################################
##################################      IDENTIFICATION OF MARKERS BETWEEN ALL CLUSTERS     ##################################################
################################################################################################################################################

# Assign identity of resolution 0.7 to clusters ---------------------------------------------
Idents(object = seurat_proteome) <- "RNA_snn_res.0.8"

# Compare all clusters ----------------------------
markers_all_initial <- FindAllMarkers(seurat_proteome, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

############# Conclusion: FindMarkers() does not work for proteomics data #############

proteome_clusters <- seurat_proteome@meta.data
write.csv(proteome_clusters, file = "~/single_fiber_heterogeneity/data/proteomics clustering/proteome_clusters_6PCs.csv")
