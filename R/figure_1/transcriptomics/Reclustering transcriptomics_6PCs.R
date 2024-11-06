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
library(VennDiagram)
library(clusterProfiler)


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/1 Clustering/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata")

# Assign identity of resolution 0.7 to clusters ---------------------------------------------
Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "SCT_snn_res.0.7"

# Reclustering
filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest <- RenameIdents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                                                             "0" = "Fast1",
                                                                             "1" = "Fast3",
                                                                             "2" = "Slow1",
                                                                             "3" = "Fast2",
                                                                             "4" = "Slow1",
                                                                             "5" = "Slow2",
                                                                             "6" = "Intermediate1")

# Reclustering to get in logical order
filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest <- RenameIdents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                                                   "Fast1" = "Fast1",
                                                                   "Fast2" = "Fast2",
                                                                   "Fast3" = "Fast3",
                                                                   "Intermediate1" = "Intermediate1",
                                                                   "Slow1" = "Slow1",
                                                                   "Slow2" = "Slow2")

# Add final clustering to metadata
clustering <- Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest)
df_cluster <- data.frame(cluster = clustering)

all(rownames(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data) %in% rownames(df_cluster))
all(rownames(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data) == rownames(df_cluster))

filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$final_cluster <- df_cluster$cluster


# Set identify of clusters
Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "final_cluster"

# Load in gene annotation file --------------------------------------------
annotations <- read.csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

# Set to SCT assay for visualization, and log-transform data ------------------------------
DefaultAssay(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "SCT"



################################################################################################################################################
################################################      FEATURE PLOTS    ################################################################
################################################################################################################################################


# Fiber type only MYH7 and MYH2 --------------------------------------------------------------
featureplot_fibertype1 <- FeaturePlot(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                      features = c("MYH2", "MYH7"),
                                      order = TRUE,
                                      min.cutoff = 'q10',
                                      label = F,
                                      repel = TRUE,
                                      reduction = "umap",
                                      ncol = 2,
                                      pt.size = 0.25,
                                      combine = FALSE)

featureplot_fibertype2 <- lapply(X = featureplot_fibertype1, FUN = function(x) x +
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

featureplot_fibertype3 <- CombinePlots(plots = featureplot_fibertype2)

ggsave(featureplot_fibertype3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_all_genes/6_PCs/Featureplot_fibertype_transcriptomics.png", width = 135, height = 60, units="mm")




# Feature plot for developmental MYHs  --------------------------------------------------------------

        # MYH3, MYH8, MYL4 not expressed





################################################################################################################################################
###############################################      CLUSTERING QUALITY CONTROL     ############################################################
################################################################################################################################################

# By subject
seurat_by_subject <- filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest

seurat_by_subject@meta.data <- seurat_by_subject@meta.data %>%
    dplyr::mutate(
        subject = dplyr::case_when(
            subject == "1" ~ "T1",
            subject == "2" ~ "T2",
            subject == "3" ~ "T3",
            subject == "4" ~ "T4",
            subject == "5" ~ "T5",
            subject == "6" ~ "T6",
            subject == "7" ~ "T7",
            subject == "8" ~ "T8",
            subject == "9" ~ "T9",
            subject == "10" ~ "T10",
            subject == "11" ~ "T11",
            subject == "12" ~ "T12",
            subject == "13" ~ "T13",
            subject == "14" ~ "T14",
            TRUE ~ "NA"
    ))

seurat_by_subject@meta.data$subject <- factor(seurat_by_subject@meta.data$subject, levels = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14"))

UMAP_subject <- Seurat::DimPlot(seurat_by_subject,
                                label = FALSE,
                                reduction = "umap",
                                pt.size = 0.3,
                                group.by = "subject") +
    ggtitle("UMAP Transcriptomics - by subject") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    guides(color = guide_legend(override.aes = list(size=1), ncol=2)) +
    scale_color_manual(values = c("#30123BFF", "#424AB3FF", "#467EF4FF", "#31AFF5FF", "#18DAC7FF", "#38F491FF", "#83FF52FF", "#BDF534FF", "#E9D539FF", "#FEAA33FF", "#F8721CFF", "#E03F08FF", "#B61C02FF", "#7A0403FF")) +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.text=element_text(size=4)
    )

ggsave(UMAP_subject, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/UMAP_transcriptomics_subject.png", width = 90, height = 60, units="mm")


# By condition
UMAP_condition <- Seurat::DimPlot(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                  label = FALSE,
                                  reduction = "umap",
                                  pt.size = 0.3,
                                  group.by = "condition") +
    ggtitle("UMAP Transcriptomics - by day") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    guides(color = guide_legend(override.aes = list(size=1))) +
    scale_color_manual(labels = c("Day A", "Day B", "Day C"), values = c("#B61C02FF", "#31AFF5FF", "#30123BFF")) +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.text=element_text(size=4)
    )

ggsave(UMAP_condition, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/UMAP_transcriptomics_condition.png", width = 90, height = 60, units="mm")


# By fiber type
UMAP_fibertype <- Seurat::DimPlot(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                  label = FALSE,
                                  reduction = "umap",
                                  pt.size = 0.3,
                                  group.by = "fiber_type_MYH_hybrids") +
    ggtitle("UMAP Transcriptomics - by fiber type") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    guides(color = guide_legend(override.aes = list(size=1))) +
    scale_color_manual("Fiber type", labels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"), values = c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325", "#D2631C")) +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.text=element_text(size=4)
    )

ggsave(UMAP_fibertype, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/UMAP_transcriptomics_fibertype.png", width = 60, height = 60, units="mm")


################################################################################################################################################
################################################      SAVE NEW SEURAT OBJECT    ################################################################
################################################################################################################################################

filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest <- filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest
save(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, file="8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

