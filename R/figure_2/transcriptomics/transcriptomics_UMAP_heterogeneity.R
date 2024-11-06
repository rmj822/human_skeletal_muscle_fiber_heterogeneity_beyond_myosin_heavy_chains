
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
library(VennDiagram)
library(clusterProfiler)


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

# Set identify of clusters
Idents(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"

# Load in gene annotation file --------------------------------------------
annotations <- read.csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

# Set to SCT assay for visualization, and log-transform data ------------------------------
DefaultAssay(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "SCT"


################################################################################################################################################
##############################################      MAKE DIM REDUCTION PLOTS     ###############################################################
################################################################################################################################################

scales::viridis_pal(option = "plasma")(length(unique(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data$final_cluster)))

# Export UMAP, TSNE AND PCA plots with correct clustering
TSNE <- DimPlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                reduction = "tsne",
                label = FALSE,
                pt.size = 1,
                cols = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF")) +
    guides(color = guide_legend(override.aes = list(size=1))) +
    theme_classic() +
    ggtitle("tSNE Transcriptomics") +
    xlab("tSNE 1") +
    ylab("tSNE 2") +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=4)
    )

ggsave(TSNE, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/tSNE_transcriptomics.png",  width = 90, height = 60, units="mm")

UMAP <- DimPlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                reduction = "umap",
                label = FALSE,
                pt.size = 1,
                cols = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                group.by = "final_cluster") +
    guides(color = guide_legend(override.aes = list(size=1))) +
    theme_classic() +
    ggtitle("UMAP Transcriptomics") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=4)
    )

ggsave(UMAP, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/UMAP_transcriptomics.png",  width = 90, height = 60, units="mm")

PCA <- DimPlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
               reduction = "pca",
               label = FALSE,
               pt.size = 0.3) +
    guides(color = guide_legend(override.aes = list(size=1))) +
    scale_color_manual("Cluster", values = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF")) +
    theme_classic() +
    ggtitle("PCA Transcriptomics - by cluster") +
    xlab("PC1 (11.1%)") +
    ylab("PC2 (3.5%)") +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=4),
        legend.position = "right"
    )

ggsave(PCA, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/PCA_transcriptomics.png",  width = 90, height = 60, units="mm")

################################################################################################################################################
##################################      IDENTIFICATION OF MARKERS BETWEEN NEW CLUSTERS     ##################################################
################################################################################################################################################


# Compare cluster 0 and 1 (both fast clusters) ----------------------------
markers_all <- FindAllMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75)
markers_all <- markers_all %>%
    left_join(annotations, by = c("gene" = "GENENAME"))
write_csv(markers_all, file = "~/single_fiber_heterogeneity/data/transcriptomics_clustering/transcriptomics_markers_all_6PCs.csv")

# Determine top10 genes per cluster
top20 <- markers_all %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10 <- markers_all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- markers_all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write_csv(top10, file = "~/single_fiber_heterogeneity/data/transcriptomics_clustering/transcriptomics_top10_markers_6PCs.csv")

# Make heatmap for top 10 genes
heatmap_top10 <- Seurat::DoHeatmap(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, features = top10$gene, raster=F,
                                   group.colors = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                                   assay = "SCT", slot = "scale.data", size=1, angle = 45, hjust = 0.1) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749")) +
    theme(
        text = element_text(face="bold", colour="black", size=4.5),
        axis.text.y = element_text(size = 3),
        legend.position = "top",
        legend.key.size = unit(4, 'mm')
    ) +
    guides(color="none")

heatmap_top10_nolegend <- Seurat::DoHeatmap(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, features = top10$gene, raster=F,
                                            group.colors = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                                            assay = "SCT", slot = "scale.data", size=1, angle = 45, hjust = 0.1) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749")) +
    theme(
        text = element_text(face="bold", colour="black", size=4.5),
        axis.text.y = element_text(size = 3),
        legend.position = "none",
        legend.key.size = unit(4, 'mm')
    ) +
    guides(color="none")

ggsave(heatmap_top10_nolegend, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/heatmap_transcriptomics_top10.png",  width = 90, height = 60, units="mm")

heatmap_legend <- ggpubr::get_legend(heatmap_top10)
heatmap_legend_ggplot <- ggpubr::as_ggplot(heatmap_legend)

ggsave(heatmap_legend_ggplot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/heatmap_legend_transcriptomics_top10.png",  width = 90, height = 60, units="mm")


# Make heatmap for top 5 genes
heatmap_top5 <- Seurat::DoHeatmap(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, features = top5$gene, raster=F,
                                  group.colors = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                                  assay = "SCT", slot = "scale.data",size=2, angle = 25, hjust = -0.05) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749")) +
    theme(
        text = element_text(face="bold", colour="black", size=4.5),
        axis.text.y = element_text(size = 4),
        legend.position = "top",
        legend.key.size = unit(4, 'mm')
    ) +
    guides(color="none")

heatmap_top5_nolegend <- Seurat::DoHeatmap(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, features = top5$gene, raster=F,
                                           group.colors = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                                           assay = "SCT", slot = "scale.data", size=2, angle = 25, hjust = -0.05) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749")) +
    theme(
        text = element_text(face="bold", colour="black", size=4.5),
        axis.text.y = element_text(size = 4),
        legend.position = "none",
        legend.key.size = unit(4, 'mm')
    ) +
    guides(color="none")

ggsave(heatmap_top5_nolegend, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/heatmap_transcriptomics_top5.png",  width = 90, height = 60, units="mm")

heatmap_legend_top5 <- ggpubr::get_legend(heatmap_top5)
heatmap_legend_top5_ggplot <- ggpubr::as_ggplot(heatmap_legend_top5)

ggsave(heatmap_legend_top5_ggplot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/heatmap_legend_transcriptomics_top5.png",  width = 90, height = 60, units="mm")


################################################################################################################################################
################################################      SEGREGATION PER nUMI / nGene     ###########################################################
################################################################################################################################################

featureplot_nUMI1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                 features = c("nUMI", "nGene"),
                                 order = TRUE,
                                 min.cutoff = 'q10',
                                 label = F,
                                 repel = TRUE,
                                 reduction = "umap",
                                 cols = c("white", "#354D3C"),
                                 ncol = 2,
                                 pt.size = 0.25,
                                 combine = FALSE)

featureplot_nUMI2 <- lapply(X = featureplot_nUMI1, FUN = function(x) x +
                                theme_classic() +
                                theme(
                                    text = element_text(face="bold", colour="black", size=4),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, size=8),
                                    legend.position = "right",
                                    legend.key.size = unit(4, 'mm')
                                ) +
                                xlab("tSNE 1") +
                                ylab("tSNE 2")
)

featureplot_nUMI3 <- CombinePlots(plots = featureplot_nUMI2)

ggsave(featureplot_nUMI3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_nUMI_transcriptomics.png", width = 135, height = 90, units="mm")


################################################################################################################################################
################################################      FEATURE PLOTS    ################################################################
################################################################################################################################################

# NMJ markers --------------------------------------------------------------
featureplot_NMJ1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                features = c("CHRNA1", "LRP4", "COLQ"),
                                order = TRUE,
                                min.cutoff = 'q10',
                                label = F,
                                repel = TRUE,
                                reduction = "umap",
                                ncol = 2,
                                pt.size = 0.25,
                                combine = FALSE)

featureplot_NMJ2 <- lapply(X = featureplot_NMJ1, FUN = function(x) x +
                               theme_classic() +
                               scale_colour_viridis() +
                               theme(
                                   text = element_text(face="bold", colour="black", size=4),
                                   axis.text = element_blank(),
                                   axis.ticks = element_blank(),
                                   plot.title = element_text(hjust = 0.5, size=8),
                                   legend.position = "right",
                                   legend.key.size = unit(4, 'mm')
                               ) +
                               xlab("UMAP1") +
                               ylab("UMAP2")
)

featureplot_NMJ3 <- CombinePlots(plots = featureplot_NMJ2)

ggsave(featureplot_NMJ3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_NMJ_transcriptomics.png", width = 135, height = 90, units="mm")

# FAP markers --------------------------------------------------------------
featureplot_fibro1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                  features = c("APOD", "DCN", "CFD", "CXCL14"),
                                  order = TRUE,
                                  min.cutoff = 'q10',
                                  label = F,
                                  repel = TRUE,
                                  reduction = "umap",
                                  ncol = 2,
                                  pt.size = 0.25,
                                  combine = FALSE)

featureplot_fibro2 <- lapply(X = featureplot_fibro1, FUN = function(x) x +
                                 theme_classic() +
                                 scale_colour_viridis() +
                                 theme(
                                     text = element_text(face="bold", colour="black", size=4),
                                     axis.text = element_blank(),
                                     axis.ticks = element_blank(),
                                     plot.title = element_text(hjust = 0.5, size=8),
                                     legend.position = "right",
                                     legend.key.size = unit(4, 'mm')
                                 ) +
                                 xlab("UMAP1") +
                                 ylab("UMAP2")
)

featureplot_fibro3 <- CombinePlots(plots = featureplot_fibro2)

ggsave(featureplot_fibro3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_fibro_transcriptomics.png", width = 135, height = 90, units="mm")

# MTJ markers --------------------------------------------------------------
featureplot_MTJ1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                features = c("SORBS2", "ANKRD1"),
                                order = TRUE,
                                min.cutoff = 'q10',
                                label = F,
                                repel = TRUE,
                                reduction = "umap",
                                ncol = 2,
                                pt.size = 0.25,
                                combine = FALSE)

featureplot_MTJ2 <- lapply(X = featureplot_MTJ1, FUN = function(x) x +
                               theme_classic() +
                               scale_colour_viridis() +
                               theme(
                                   text = element_text(face="bold", colour="black", size=4),
                                   axis.text = element_blank(),
                                   axis.ticks = element_blank(),
                                   plot.title = element_text(hjust = 0.5, size=8),
                                   legend.position = "right",
                                   legend.key.size = unit(4, 'mm')
                               ) +
                               xlab("UMAP1") +
                               ylab("UMAP2")
)

featureplot_MTJ3 <- CombinePlots(plots = featureplot_MTJ2)

ggsave(featureplot_MTJ3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_MTJ_transcriptomics.png", width = 135, height = 90, units="mm")

# develop markers --------------------------------------------------------------
featureplot_develop1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("ENAH", "NOS1", "MYOD1", "MYOG", "NRAP", "FHOD3", "FLNC"),
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_develop2 <- lapply(X = featureplot_develop1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_develop3 <- CombinePlots(plots = featureplot_develop2)

ggsave(featureplot_develop3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_develop_transcriptomics.png", width = 135, height = 90, units="mm")

# EC markers --------------------------------------------------------------
featureplot_EC1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("VWF", "FABP4", "CLDN5", "AQP1", "PECAM1", "A2M", "BTNL9"),
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_EC2 <- lapply(X = featureplot_EC1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_EC3 <- CombinePlots(plots = featureplot_EC2)

ggsave(featureplot_EC3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_EC_transcriptomics.png", width = 135, height = 90, units="mm")

# immune markers --------------------------------------------------------------
featureplot_immune1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("IL32"), # IL32 is the only one detected
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_immune2 <- lapply(X = featureplot_immune1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_immune3 <- CombinePlots(plots = featureplot_immune2)

ggsave(featureplot_immune3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_immune_transcriptomics.png", width = 135, height = 90, units="mm")

# MoMac markers --------------------------------------------------------------
featureplot_MoMac1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("HLA-DRA", "HLA-DPA1", "RNASE1", "CD74"),
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_MoMac2 <- lapply(X = featureplot_MoMac1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_MoMac3 <- CombinePlots(plots = featureplot_MoMac2)

ggsave(featureplot_MoMac3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_MoMac_transcriptomics.png", width = 135, height = 90, units="mm")

# SMC/Pericytes --------------------------------------------------------------
featureplot_SMC1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("MYL9", "TAGLN", "MYH11", "RGS5", "IGFBP5"),
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_SMC2 <- lapply(X = featureplot_SMC1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_SMC3 <- CombinePlots(plots = featureplot_SMC2)

ggsave(featureplot_SMC3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_SMC_transcriptomics.png", width = 135, height = 90, units="mm")

# MuSC markers --------------------------------------------------------------
featureplot_MuSC1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("APOC1", "APOE", "DLK1", "MYF5", "SPATS2L", "PAX7", "CHRDL2"), # Almost no gene detected
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_MuSC2 <- lapply(X = featureplot_MuSC1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_MuSC3 <- CombinePlots(plots = featureplot_MuSC2)

ggsave(featureplot_MuSC3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_MuSC_transcriptomics.png", width = 135, height = 90, units="mm")

# neutrophils markers --------------------------------------------------------------
featureplot_neutrophils1 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    features = c("NAMPT", "G0S2"),
                                    order = TRUE,
                                    min.cutoff = 'q10',
                                    label = F,
                                    repel = TRUE,
                                    reduction = "umap",
                                    ncol = 2,
                                    pt.size = 0.25,
                                    combine = FALSE)

featureplot_neutrophils2 <- lapply(X = featureplot_neutrophils1, FUN = function(x) x +
                                   theme_classic() +
                                   scale_colour_viridis() +
                                   theme(
                                       text = element_text(face="bold", colour="black", size=4),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),
                                       plot.title = element_text(hjust = 0.5, size=8),
                                       legend.position = "right",
                                       legend.key.size = unit(4, 'mm')
                                   ) +
                                   xlab("UMAP1") +
                                   ylab("UMAP2")
)

featureplot_neutrophils3 <- CombinePlots(plots = featureplot_neutrophils2)

ggsave(featureplot_neutrophils3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_neutrophils_transcriptomics.png", width = 135, height = 90, units="mm")

# Fast 2 and Slow 2 markers  --------------------------------------------------------------
featureplot_fast2_slow2a <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                        features = c("RIMS2", "KAZN", "FAM155A", "RP11-129M6.1"),
                                        order = TRUE,
                                        min.cutoff = 'q10',
                                        label = F,
                                        repel = TRUE,
                                        reduction = "umap",
                                        ncol = 2,
                                        pt.size = 0.25,
                                        combine = FALSE)

featureplot_fast2_slow2b <- lapply(X = featureplot_fast2_slow2a, FUN = function(x) x +
                                       theme_classic() +
                                       scale_colour_viridis(option = "plasma") +
                                       theme(
                                           text = element_text(face="bold", colour="black", size=4),
                                           axis.text = element_blank(),
                                           axis.ticks = element_blank(),
                                           plot.title = element_text(hjust = 0.5, size=8),
                                           legend.position = "right",
                                           legend.key.size = unit(4, 'mm')
                                       ) +
                                       xlab("UMAP1") +
                                       ylab("UMAP2")
)

featureplot_fast2_slow2c <- CombinePlots(plots = featureplot_fast2_slow2b)

ggsave(featureplot_fast2_slow2c, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_fast2_slow2_transcriptomics.png", width = 135, height = 90, units="mm")

# Fast 3 markers  --------------------------------------------------------------
featureplot_fast3a <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                  features = c("NDUFAB1", "RPS29", "CHCHD10", "RPS23"),
                                  order = TRUE,
                                  min.cutoff = 'q10',
                                  label = F,
                                  repel = TRUE,
                                  reduction = "umap",
                                  ncol = 2,
                                  pt.size = 0.25,
                                  combine = FALSE)

featureplot_fast3b <- lapply(X = featureplot_fast3a, FUN = function(x) x +
                                 theme_classic() +
                                 scale_colour_viridis(option = "plasma") +
                                 theme(
                                     text = element_text(face="bold", colour="black", size=4),
                                     axis.text = element_blank(),
                                     axis.ticks = element_blank(),
                                     plot.title = element_text(hjust = 0.5, size=8),
                                     legend.position = "right",
                                     legend.key.size = unit(4, 'mm')
                                 ) +
                                 xlab("UMAP1") +
                                 ylab("UMAP2")
)

featureplot_fast3c <- CombinePlots(plots = featureplot_fast3b)

ggsave(featureplot_fast3c, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_fast3_transcriptomics.png", width = 135, height = 90, units="mm")

# Feature plots for main figure 2 of manuscript - v1 --------------------------------------------------------------

# Without legends
feature_plot_TCF7L2 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("TCF7L2"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_SLIT3 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("SLIT3"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_MAML2 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                          features = c("MAML2"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_NDUFAB1 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("NDUFAB1"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_CHCHD10 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                            features = c("CHCHD10"),
                                            pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")


feature_plots_UMAP <- ggpubr::ggarrange(
    feature_plot_TCF7L2,
    feature_plot_SLIT3,
    feature_plot_MAML2,
    feature_plot_NDUFAB1,
    feature_plot_CHCHD10,
    ncol = 5,
    nrow = 1)

feature_plots_UMAP <- annotate_figure(feature_plots_UMAP, top = text_grob("Transcriptomics", color = "black", face = "bold", size = 8))

ggsave(feature_plots_UMAP,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_fig4_transcriptomics.png",
       width = 195,
       height = 40,
       units="mm")

# Get legends

# TCF7L2
feature_plot_TCF7L2_legend <- feature_plot_TCF7L2 +
    ggplot2::theme(
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_TCF7L2_legend <- ggpubr::get_legend(feature_plot_TCF7L2_legend)
feature_plot_TCF7L2_legend <- ggpubr::as_ggplot(feature_plot_TCF7L2_legend)

ggsave(feature_plot_TCF7L2_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_transcriptomics_TCF7L2_legend.png",
       width = 90,
       height = 40,
       units="mm")

# SLIT3
feature_plot_SLIT3_legend <- feature_plot_SLIT3 +
    ggplot2::theme(
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_SLIT3_legend <- ggpubr::get_legend(feature_plot_SLIT3_legend)
feature_plot_SLIT3_legend <- ggpubr::as_ggplot(feature_plot_SLIT3_legend)

ggsave(feature_plot_SLIT3_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_transcriptomics_SLIT3_legend.png",
       width = 90,
       height = 40,
       units="mm")

# MAML2
feature_plot_MAML2_legend <- feature_plot_MAML2 +
    ggplot2::theme(
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_MAML2_legend <- ggpubr::get_legend(feature_plot_MAML2_legend)
feature_plot_MAML2_legend <- ggpubr::as_ggplot(feature_plot_MAML2_legend)

ggsave(feature_plot_MAML2_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_transcriptomics_MAML2_legend.png",
       width = 90,
       height = 40,
       units="mm")

# NDUFAB1
feature_plot_NDUFAB1_legend <- feature_plot_NDUFAB1 +
    ggplot2::theme(
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_NDUFAB1_legend <- ggpubr::get_legend(feature_plot_NDUFAB1_legend)
feature_plot_NDUFAB1_legend <- ggpubr::as_ggplot(feature_plot_NDUFAB1_legend)

ggsave(feature_plot_NDUFAB1_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_transcriptomics_NDUFAB1_legend.png",
       width = 90,
       height = 40,
       units="mm")

# CHCHD10
feature_plot_CHCHD10_legend <- feature_plot_CHCHD10 +
    ggplot2::theme(
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_CHCHD10_legend <- ggpubr::get_legend(feature_plot_CHCHD10_legend)
feature_plot_CHCHD10_legend <- ggpubr::as_ggplot(feature_plot_CHCHD10_legend)

ggsave(feature_plot_CHCHD10_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_transcriptomics_CHCHD10_legend.png",
       width = 90,
       height = 40,
       units="mm")

# Feature plots for main figure 2 of manuscript - v2 --------------------------------------------------------------

feature_plot_SLIT3 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("SLIT3"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin = ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_CHCHD10 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("CHCHD10"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin= ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_CTDNEP1 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("CTDNEP1"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin= ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

final_plot <- ggpubr::ggarrange(feature_plot_CHCHD10,
                                feature_plot_SLIT3,
                                feature_plot_CTDNEP1,
                                ncol = 3,
                                nrow = 1) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

final_plot <- ggpubr::annotate_figure(final_plot, top = ggpubr::text_grob("Transcriptomics",
                                                            color = "black", face = "bold", size = 8))

ggsave(final_plot,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_transcriptomics_fig2_v2.png",
       width = 120,
       height = 35,
       units="mm")



################################################################################################################################################
################################################      HEATMAPS SLOW2 vs SLOW1 and FAST2 vs FAST1    ################################################################
################################################################################################################################################

seurat_averaged <- AverageExpression(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                     return.seurat = T, assays = "SCT")

# MYH7
feature_plot_MYH7_with_legend <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                                     features = c("MYH7"),
                                                     pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "top",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_legend_MYH7 <- ggpubr::get_legend(feature_plot_MYH7_with_legend)
feature_legend_MYH7 <- ggpubr::as_ggplot(feature_legend_MYH7)

ggsave(feature_legend_MYH7,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/Featureplot_MYH7_MYH2_MYH1_transcriptomics_legend_MYH7.png",
       width = 90,
       height = 40,
       units="mm")Idents(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"

# Slow1 vs Slow2 ----------------------------

# Find markers
markers_Slow1vsSlow2 <- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    ident.1 = "Slow1",
                                    ident.2 = "Slow2")

markers_Slow1vsSlow2 <- markers_Slow1vsSlow2 %>%
    rownames_to_column(var = "gene")

# Add gene annotations
markers_Slow1vsSlow2 <- left_join(markers_Slow1vsSlow2, annotations, by = c("gene" = "GENENAME"))

# Make column with pct diff
markers_Slow1vsSlow2$pct.diff <- markers_Slow1vsSlow2$pct.1 - markers_Slow1vsSlow2$pct.2

# Save markers
write_csv(markers_Slow1vsSlow2, file = "~/single_fiber_heterogeneity/data/transcriptomics_clustering/Slow1_vs_Slow2_markers.csv")

# Fast1 vs Fast1 ----------------------------

# Find markers
markers_Fast1vsFast2 <- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    ident.1 = "Fast1",
                                    ident.2 = "Fast2")

markers_Fast1vsFast2 <- markers_Fast1vsFast2 %>%
    rownames_to_column(var = "gene")

# Add gene annotations
markers_Fast1vsFast2 <- left_join(markers_Fast1vsFast2, annotations, by = c("gene" = "GENENAME"))

# Make column with pct diff
markers_Fast1vsFast2$pct.diff <- markers_Fast1vsFast2$pct.1 - markers_Fast1vsFast2$pct.2

# Save markers
write_csv(markers_Fast1vsFast2, file = "~/single_fiber_heterogeneity/data/transcriptomics_clustering/Fast1_vs_Fast2_markers.csv")

# Combine two marker results outputs ----------------------------

genes_slow1vsslow2 <- markers_Slow1vsSlow2$gene

markers_Fast1vsFast2_overlap <- markers_Fast1vsFast2 %>% dplyr::filter(gene %in% genes_slow1vsslow2)
markers_Slow1vsSlow2_overlap <- markers_Slow1vsSlow2 %>% dplyr::filter(gene %in% markers_Fast1vsFast2_overlap$gene)

all(markers_Fast1vsFast2_overlap$gene %in% markers_Slow1vsSlow2_overlap$gene)

markers_Fast1vsFast2_overlap <- markers_Fast1vsFast2_overlap[ order(match(markers_Fast1vsFast2_overlap$gene, markers_Slow1vsSlow2_overlap$gene)), ]

all(markers_Fast1vsFast2_overlap$gene %in% markers_Slow1vsSlow2_overlap$gene)
all(markers_Fast1vsFast2_overlap$gene == markers_Slow1vsSlow2_overlap$gene)

combined_markers <- tibble(
    gene = markers_Fast1vsFast2_overlap$gene,
    slow = markers_Slow1vsSlow2_overlap$avg_log2FC,
    fast = markers_Fast1vsFast2_overlap$avg_log2FC
)

combined_markers$mean_LFC <- ((combined_markers$slow + combined_markers$fast) / 2)

combined_markers <- combined_markers %>% arrange(desc(abs(mean_LFC)))

# Filter top 50 different genes based on average LFC
combined_markers_top50 <- combined_markers %>% top_n(n = 50, wt = abs(mean_LFC))


# Create heatmap ----------------------------------------------------------
heatmap_clusters <- Seurat::DoHeatmap(seurat_averaged,
                                      features = combined_markers_top50$gene,
                                      raster=F,
                                      draw.lines=F,
                                      group.colors = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                                      assay = "SCT",
                                      slot = "scale.data",
                                      size=2,
                                      angle = 0,
                                      hjust = 0.5
) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749"), na.value = "white") +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text.y = element_text(size = 5),
        legend.position = "top",
        legend.key.size = unit(4, 'mm')
    ) +
    guides(color="none")

heatmap_clusters_nolegend <- Seurat::DoHeatmap(seurat_averaged,
                                               features = combined_markers_top50$gene,
                                               raster=F,
                                               draw.lines=F,
                                               group.colors = c("#30123BFF", "#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"),
                                               assay = "SCT",
                                               slot = "scale.data",
                                               size=2,
                                               angle = 0,
                                               hjust = 0.5
) +
    scale_fill_gradientn(colors = c("#134057", "white", "#BC4749"), na.value = "white") +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text.y = element_text(size = 5),
        legend.position = "none",
        legend.key.size = unit(4, 'mm')
    ) +
    guides(color="none")

ggsave(heatmap_clusters_nolegend, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/heatmap_transcriptomics_clusters.png",  width = 128, height = 90, units="mm")

heatmap_clusters_legend <- ggpubr::get_legend(heatmap_clusters)
heatmap_clusters_legend_ggplot <- ggpubr::as_ggplot(heatmap_clusters_legend)

ggsave(heatmap_clusters_legend_ggplot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/heatmap_legend_clusters_transcriptomics.png",  width = 128, height = 90, units="mm")


################################################################################################################################################
##############################################     VENN PLOT OVERLAP CLUSTER MARKERS    ########################################################
################################################################################################################################################

markers_Fast1vsFast2_LFC1 <- markers_Fast1vsFast2 %>% dplyr::filter(avg_log2FC > 1 | avg_log2FC < -1)
markers_Slow1vsSlow2_LFC1 <- markers_Slow1vsSlow2 %>% dplyr::filter(avg_log2FC > 1 | avg_log2FC < -1)

Venn_clusters <- venn.diagram(
    # General
    filename=NULL,
    disable.logging=T,
    x = list(
        markers_Fast1vsFast2_LFC1 %>% as.data.frame() %>% dplyr::select(gene) %>% unlist(use.names=F),
        markers_Slow1vsSlow2_LFC1 %>% as.data.frame() %>% dplyr::select(gene) %>% unlist(use.names=F)
    ),
    category.names = c("Fast2" , "Slow2"),
    main = "Markers per cluster",
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.cex = 0.5,


    # Circles
    lwd = 2,
    col=c("#F49D6E", "#8FB339"),
    fill = c(alpha("#F49D6E",0.3), alpha('#8FB339',0.3)),

    # Numbers
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans",
    cat.distance = c(0.05, 0.02),

    # Names
    cat.cex = 0.75,
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(0,0),
    cat.dist = c(-0.09, -0.18),
    cat.col = c("#F49D6E", "#8FB339")
)

ggsave(Venn_clusters, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/venn_clusters.png", width = 40, height = 40, units="mm")

################################################################################################################################################
################################################      BIOTYPE OF MARKER GENES    ###############################################################
################################################################################################################################################

# Slow2 -----------------------------------------------

# Determine how many genes per type
slow2_type <- markers_Slow1vsSlow2_LFC1 %>%
    group_by(GENEBIOTYPE) %>%
    summarise(n = n())

# Filter
slow2_type <- slow2_type %>% dplyr::filter(GENEBIOTYPE == "antisense" | GENEBIOTYPE == "lincRNA" | GENEBIOTYPE == "protein_coding")

# Compute percentages
slow2_type$fraction = slow2_type$n / sum(slow2_type$n)

# Compute the cumulative percentages (top of each rectangle)
slow2_type$ymax = cumsum(slow2_type$fraction)

# Compute the bottom of each rectangle
slow2_type$ymin = c(0, head(slow2_type$ymax, n=-1))

# Compute label position
slow2_type$labelPosition.y <- c(0.05, 0.45, 0.72)

slow2_type$labelPosition.x <- c(1.5,1.9,1.1)

# Compute a good label
slow2_type$label <- paste0(slow2_type$GENEBIOTYPE, "\n n: ", slow2_type$n)

# Make the plot
biotype_slow2 <- ggplot(slow2_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=GENEBIOTYPE)) +
    geom_rect() +
    geom_label(aes(x= labelPosition.x, y=labelPosition.y, label=label), size=2, label.padding= unit(0.1, "lines")) +
    scale_fill_brewer(palette=4) +
    coord_polar(theta="y") +
    xlim(c(-0.1, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle("Slow2 markers (LFC > 1)") +
    theme(
        text = element_text(face="bold", size=8, colour="black"),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
    )

ggsave(biotype_slow2, file="~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/biotype_slow2.png", width = 40, height = 60, units="mm")


# fast2 -----------------------------------------------

# Determine how many genes per type
fast2_type <- markers_Fast1vsFast2_LFC1 %>%
    group_by(GENEBIOTYPE) %>%
    summarise(n = n())

# Filter
fast2_type <- fast2_type %>% dplyr::filter(GENEBIOTYPE == "antisense" | GENEBIOTYPE == "lincRNA" | GENEBIOTYPE == "protein_coding")

# Compute percentages
fast2_type$fraction = fast2_type$n / sum(fast2_type$n)

# Compute the cumulative percentages (top of each rectangle)
fast2_type$ymax = cumsum(fast2_type$fraction)

# Compute the bottom of each rectangle
fast2_type$ymin = c(0, head(fast2_type$ymax, n=-1))

# Compute label position
fast2_type$labelPosition.y <- c(0.01, 0.23, 0.5)

fast2_type$labelPosition.x <- c(1.8,1.7,1.4)

# Compute a good label
fast2_type$label <- paste0(fast2_type$GENEBIOTYPE, "\n n: ", fast2_type$n)

# Make the plot
biotype_fast2 <- ggplot(fast2_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=GENEBIOTYPE)) +
    geom_rect() +
    geom_label(aes(x= labelPosition.x, y=labelPosition.y, label=label), size=2, label.padding= unit(0.1, "lines")) +
    scale_fill_brewer(palette=4) +
    coord_polar(theta="y") +
    xlim(c(-0.1, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle("Fast2 markers (LFC > 1)") +
    theme(
        text = element_text(face="bold", size=8, colour="black"),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
    )

ggsave(biotype_fast2, file="~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/biotype_fast2.png", width = 40, height = 60, units="mm")

################################################################################################################################################
################################################      ENRICHMENT FOR CLUSTER MARKERS   #########################################################
################################################################################################################################################

# Slow
markers_slow1vsslow2 <- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    ident.1 = "Slow1",
                                    ident.2 = "Slow2")

markers_slow2 <- markers_slow1vsslow2 %>% dplyr::filter(avg_log2FC < -0.58)
markers_slow1 <- markers_slow1vsslow2 %>% dplyr::filter(avg_log2FC > 0.58)


# Fast
markers_fast1vsfast2<- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                   ident.1 = "Fast1",
                                   ident.2 = "Fast2")

markers_fast1vsfast3 <- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    ident.1 = "Fast1",
                                    ident.2 = "Fast3")

markers_fast2 <- markers_fast1vsfast2 %>% dplyr::filter(avg_log2FC < -0.58)
markers_fast3<- markers_fast1vsfast3 %>% dplyr::filter(avg_log2FC < -0.58)

# Intermediate
markers_all <- FindAllMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
markers_intermediate1 <- markers_all %>% dplyr::filter(cluster == "Intermediate1") %>% dplyr::filter(avg_log2FC > 0.58)

# Get background gene set
allgenes_SCT <- rownames(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest)

# Enrichment for Slow2
ORA_slow2 <- enrichGO(gene = rownames(markers_slow2),
                      universe      = allgenes_SCT,
                      keyType       = "SYMBOL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)

ORA_slow2_simplify <- simplify(ORA_slow2, cutoff=0.7, by="p.adjust", select_fun=min)

ORA_slow2_simplify_df <- as.data.frame(ORA_slow2_simplify)

ORA_slow2_simplify_df <- mutate(ORA_slow2_simplify_df, foldEnrich =
                                    (as.numeric(sub("/\\d+", "", ORA_slow2_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_slow2_simplify_df$GeneRatio))) /
                                    (as.numeric(sub("/\\d+", "", ORA_slow2_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_slow2_simplify_df$BgRatio)))
)

write_csv(ORA_slow2_simplify_df, file="~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_slow2.csv")

# Enrichment for fast2
ORA_fast2 <- enrichGO(gene = rownames(markers_fast2),
                      universe      = allgenes_SCT,
                      keyType       = "SYMBOL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)

ORA_fast2_simplify <- simplify(ORA_fast2, cutoff=0.7, by="p.adjust", select_fun=min)

ORA_fast2_simplify_df <- as.data.frame(ORA_fast2_simplify)

ORA_fast2_simplify_df <- mutate(ORA_fast2_simplify_df, foldEnrich =
                                    (as.numeric(sub("/\\d+", "", ORA_fast2_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_fast2_simplify_df$GeneRatio))) /
                                    (as.numeric(sub("/\\d+", "", ORA_fast2_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_fast2_simplify_df$BgRatio)))
)

write_csv(ORA_fast2_simplify_df, file="~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_fast2.csv")

# Enrichment for fast3
ORA_fast3 <- enrichGO(gene = rownames(markers_fast3),
                      universe      = allgenes_SCT,
                      keyType       = "SYMBOL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)

ORA_fast3_simplify <- simplify(ORA_fast3, cutoff=0.7, by="p.adjust", select_fun=min)

ORA_fast3_simplify_df <- as.data.frame(ORA_fast3_simplify)

ORA_fast3_simplify_df <- mutate(ORA_fast3_simplify_df, foldEnrich =
                                    (as.numeric(sub("/\\d+", "", ORA_fast3_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_fast3_simplify_df$GeneRatio))) /
                                    (as.numeric(sub("/\\d+", "", ORA_fast3_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_fast3_simplify_df$BgRatio)))
)

write_csv(ORA_fast3_simplify_df, file="~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_fast3.csv")

# Enrichment for intermediate1
ORA_intermediate1 <- enrichGO(gene = rownames(markers_intermediate1),
                              universe      = allgenes_SCT,
                              keyType       = "SYMBOL",
                              OrgDb         = org.Hs.eg.db,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.05,
                              readable      = FALSE)

ORA_intermediate1_simplify <- simplify(ORA_intermediate1, cutoff=0.7, by="p.adjust", select_fun=min)

ORA_intermediate1_simplify_df <- as.data.frame(ORA_intermediate1_simplify)

ORA_intermediate1_simplify_df <- mutate(ORA_intermediate1_simplify_df, foldEnrich =
                                            (as.numeric(sub("/\\d+", "", ORA_intermediate1_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_intermediate1_simplify_df$GeneRatio))) /
                                            (as.numeric(sub("/\\d+", "", ORA_intermediate1_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_intermediate1_simplify_df$BgRatio)))
)

write_csv(ORA_intermediate1_simplify_df, file="~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_intermediate1.csv")


# Prepare data for plotting
ORA_slow2 <- read.csv("~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_slow2.csv")
ORA_fast2 <- read.csv("~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_fast2.csv")
ORA_fast3 <- read.csv("~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_fast3.csv")
ORA_intermediate1 <- read.csv("~/single_fiber_heterogeneity/data/transcriptomics_enrichment_clusters/enrichment_intermediate1.csv")


# Select categories for slow2
slow2_categories <- ORA_slow2 %>%
    dplyr::filter(
        Description == "transmembrane signaling receptor activity" |
            Description == "DNA-binding transcription factor activity" |
            Description == "cell adhesion" |
            Description == "GTPase regulator activity" |
            Description == "cytokine production" |
            Description == "protein tyrosine kinase activator activity"
    )

# Select categories for fast2
fast2_categories <- ORA_fast2 %>%
    dplyr::filter(
        Description == "DNA-binding transcription factor activity" |
            Description == "transmembrane signaling receptor activity" |
            Description == "GTPase regulator activity" |
            Description == "glycosyltransferase activity" |
            Description == "protein kinase activity"
    )

# Select categories for fast3
fast3_categories <- ORA_fast3 %>%
    dplyr::filter(
        Description == "respiratory chain complex" |
            Description == "oxidoreductase complex" |
            Description == "polysomal ribosome" |
            Description == "cytoplasmic translation"
    )

# Select categories for fast3
intermediate1_categories <- ORA_intermediate1 %>%
    dplyr::filter(
        Description == "oxidoreductase complex" |
            Description == "respiratory chain complex" |
            Description == "cell-substrate junction" |
            Description == "cytoplasmic translation"
    )


# Add required columns for plotting
slow2_categories$cluster <- rep("Slow2", nrow(slow2_categories))
fast2_categories$cluster <- rep("Fast2", nrow(fast2_categories))
fast3_categories$cluster <- rep("Fast3", nrow(fast3_categories))
intermediate1_categories$cluster <- rep("Intermediate1", nrow(intermediate1_categories))


# Combine all into one dataframe
data_plot <- bind_rows(slow2_categories, fast2_categories, fast3_categories, intermediate1_categories)
data_plot$cluster = factor(data_plot$cluster, levels=c(
    "Intermediate1",
    "Fast3",
    "Fast2",
    "Slow2"
))
data_plot$graphdir <- c(1:2, 4:6, 3, 5, 4, 6, 7, 8, 9:12, 9, 10, 11, 12)

# Create dot plot

dotplot_enrich <- ggplot(data_plot, aes(fct_reorder(Description, graphdir), cluster)) +
    geom_point(aes(color= cluster, alpha=foldEnrich, size = -log10(p.adjust))) +
    scale_color_manual(values = c("#E1DD37FF", "#46F884FF", "#3E9BFEFF", "#7A0403FF")) +
    scale_alpha_continuous(range = c(0.2, 1)) +
    scale_size(range=c(1.5, 5)) +
    ggtitle("Enrichment by cluster") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=4),
        axis.text.y = element_text(size=5),
        axis.title = element_blank(),
        text = element_text(face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, face="bold", vjust=0, size=8),
        legend.position = "none",
        panel.grid.major.y = element_blank() ,
        panel.grid.major.x = element_line( linewidth=.1, color="grey" )
    )

ggsave(dotplot_enrich, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_UMAP/transcriptomics_enrichment_clusters.png", width = 110, height = 60, units="mm")

# Export data frame used for creating plot
write.csv(data_plot, file = "~/single_fiber_heterogeneity/data/transcriptomics_enrichment_by_cluster/enrichment_cluster.csv")
