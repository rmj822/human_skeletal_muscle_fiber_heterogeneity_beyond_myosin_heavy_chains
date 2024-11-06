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


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("7 Targeted fiber type markers/Fibers at rest/Inflection/filtered_normalized_fibertype_seurat_wo_MSTRG_rest.Rdata")

# Set assay to SCT
DefaultAssay(filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT"

################################################################################################################################################
########################################################      PCA   ############################################################################
################################################################################################################################################

#  Run PCA------------------------------------------------
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- Seurat::RunPCA(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest)

# Export PC genes
PC_genes <- Loadings(filtered_normalized_fibertype_seurat_wo_MSTRG_rest[["pca"]])

PC_genes <- as.data.frame(PC_genes)

write.csv(PC_genes, file = "~/single_fiber_heterogeneity/data/transcriptomics_PC_loadings.csv")


################################################################################################################################################
################################################      Identify significant PCS    ##############################################################
################################################################################################################################################

# Explore heatmap of PCs (try to find PC where heatmap starts to look fuzzy, not so distinct between groups) --------
pdf("~/single_fiber_heterogeneity/doc/figures/figure_4/transcriptomics_all_genes/PC heatmaps.pdf")
Seurat::DimHeatmap(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
           dims = 1:12,
           cells = 500,  # number of cells with most negative or positive PCA scores to use for plotting
           balanced = TRUE)
dev.off()

       # Number of significant PCs??? Difficult to see, around 10?


# Via factoextra package
transcriptomics_scaledata <- GetAssayData(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest, assay = "SCT", slot = "scale.data")
pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

scree_plot <- factoextra::fviz_eig(pca_object)
ggsave(scree_plot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/transcriptomics_all_genes/PC_scree_plot_fviz.png")

PCvar <- factoextra::get_eigenvalue(pca_object)
PCvar$rank <- c(1:length(PCvar$variance.percent))

# Quantitative determination elbow plateau (two metrics, choose lowest value)--------------------------------

# Metric 1: New PC only contributes 5%, and all other cumulatively contribute 90%
pct <- PCvar$variance.percent
cumu <- PCvar$cumulative.variance.percent
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 # PC 677

# Metric 2: PC where percent change to next PC is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
co2 # PC 6

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

# Use 6 PCs for clustering of transcriptomics data

 ################################################################################################################################################
####################################################      CLUSTERING      ######################################################################
################################################################################################################################################


# Graph-based clustering using K-nearest neighbor graph  ---------------------------------

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- Seurat::FindNeighbors(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                                   dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- Seurat::FindClusters(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                                  resolution = c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4))

# Resolution can be found in metadata object
filtered_normalized_fibertype_seurat_wo_MSTRG_rest@meta.data %>%
  View()


################################################################################################################################################
################################################      DIMENSIONALITY REDUCTION   ##############################################################
################################################################################################################################################

# Run UMAP ----------------------------------------------------------------
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- Seurat::RunUMAP(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                                                            dims = 1:6)

# Run T-SNE ----------------------------------------------------------------
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- RunTSNE(filtered_normalized_fibertype_seurat_wo_MSTRG_rest, dims = 1:6)



################################################################################################################################################
#################################################     TEST DIFFERENT RESOLUTIONS  ##############################################################
################################################################################################################################################


# Determine best resolution for our dataset - Resolution 0.2 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.2"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.2/UMAP clustering resolution 0.2.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.2")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.2/TSNE clustering resolution 0.2.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.2")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.2/PCA clustering resolution 0.2.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.2")
dev.off()


# Determine best resolution for our dataset - Resolution 0.4 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.4"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.4/UMAP clustering resolution 0.4.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
        reduction = "umap",
        label = TRUE,
        label.size = 6) +
  ggtitle("Resolution 0.4")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.4/TSNE clustering resolution 0.4.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.4")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.4/PCA clustering resolution 0.4.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.4")
dev.off()

# Determine best resolution for our dataset - Resolution 0.5 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.5"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.5/UMAP clustering resolution 0.5.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.5")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.5/TSNE clustering resolution 0.5.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.5")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.5/PCA clustering resolution 0.5.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.5")
dev.off()

 # Determine best resolution for our dataset - Resolution 0.6 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.6"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.6/UMAP clustering resolution 0.6.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.6")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.6/TSNE clustering resolution 0.6.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.6")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.6/PCA clustering resolution 0.6.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.6")
dev.off()

# Determine best resolution for our dataset - Resolution 0.7 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.7"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.7/UMAP clustering resolution 0.7.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.7")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.7/TSNE clustering resolution 0.7.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.7")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.7/PCA clustering resolution 0.7.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.7")
dev.off()

# Determine best resolution for our dataset - Resolution 0.8 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.8"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.8/UMAP clustering resolution 0.8.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.8")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.8/TSNE clustering resolution 0.8.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.8")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 0.8/PCA clustering resolution 0.8.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 0.8")
dev.off()


# Determine best resolution for our dataset - Resolution 1.0 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.1"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.0/UMAP clustering resolution 1.0.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.0")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.0/TSNE clustering resolution 1.0.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.0")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.0/PCA clustering resolution 1.0.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.0")
dev.off()


# Determine best resolution for our dataset - Resolution 1.2 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.1.2"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.2/UMAP clustering resolution 1.2.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.2")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.2/TSNE clustering resolution 1.2.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.2")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.2/PCA clustering resolution 1.2.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.2")
dev.off()

# Determine best resolution for our dataset - Resolution 1.4 -------------------------------

# Assign identity of clusters
Idents(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest) <- "SCT_snn_res.1.4"

# Plot the UMAP
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.4/UMAP clustering resolution 1.4.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.4")
dev.off()

# Plot the TSNE
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.4/TSNE clustering resolution 1.4.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.4")
dev.off()

# Plot the PCA
pdf("8 Fiber heterogeneity (only rested samples)/1 Clustering/Resolution 1.4/PCA clustering resolution 1.4.pdf")
Seurat::DimPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest,
                reduction = "pca",
                label = TRUE,
                label.size = 6) +
  ggtitle("Resolution 1.4")
dev.off()


################################################################################################################################################
#####################################################      SAVE SEURAT OBJECT    ###############################################################
################################################################################################################################################

    # Save Seurat object, including the clustering with different resolutions
    # Next steps are performed for each resolution seperately

# Rename Seurat object
filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest <- filtered_normalized_fibertype_seurat_wo_MSTRG_rest

# Save Seurat object ------------------------------------------------------
save(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, file="8 Fiber heterogeneity (only rested samples)/1 Clustering/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata")
