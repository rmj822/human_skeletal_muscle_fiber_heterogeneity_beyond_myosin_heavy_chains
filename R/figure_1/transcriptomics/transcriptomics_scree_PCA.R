
library(factoextra)
library(org.Hs.eg.db)

################################################################################################################################################
###########################################      GET NORMALIZED COUNTS      ####################################################################
################################################################################################################################################

# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object (in data-raw folder) ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

# Extract SCT scaled data
transcriptomics_scaledata <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "scale.data")

# Extract SCT normalized but not log-transformed counts
transcriptomics_counts <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")

# Set assay to SCT
DefaultAssay(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "SCT"


################################################################################################################################################
###########################################      PERFORM PCA WITH PRCOMP      ##################################################################
################################################################################################################################################

pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2)

################################################################################################################################################
################################################      ADD METADATA      ########################################################################
################################################################################################################################################

# Load metadata
metadata <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data

# Check if fibers in same order
all(rownames(data_pca) %in% rownames(metadata))
all(rownames(data_pca) == rownames(metadata))

# Add metadata
data_pca <- data_pca |>
    tibble::add_column(
        subject = metadata$subject,
        fiberID = rownames(metadata),
        fiber_type = metadata$fiber_type_MYH_hybrids
    )

# Relevel fiber types
data_pca$fiber_type <- factor(data_pca$fiber_type,
                              levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))


################################################################################################################################################
################################################      VARIANCE PER PC      #####################################################################
################################################################################################################################################

# Via factoextra package
factoextra::fviz_eig(pca_object)
PCvar <- factoextra::get_eigenvalue(pca_object)

# Nice Elbow plot for paper
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
# Use 6PCs for clustering

################################################################################################################################################
################################################      MAKE SCREE PLOT      #####################################################################
################################################################################################################################################

# Extract %s ----------------------------------------

variance <- factoextra::get_eig(pca_object) |>
    as.data.frame() |>
    dplyr::select(variance.percent) |>
    dplyr::mutate(rank = 1:length(variance.percent)) |>
    dplyr::slice_head(n = 50)

# Nice Elbow plot for paper

screeplot <- ggplot2::ggplot(variance,
                             ggplot2::aes(x = rank,
                                          y = variance.percent,
                                          color = rank > co2)) +
    ggplot2::annotate("rect",
                      xmin=-Inf,
                      xmax=6.5,
                      ymin=-Inf,
                      ymax=Inf,
                      alpha=0.2,
                      fill= "#4F7A5D") +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c("#4F7A5D", "#B7DFB3")) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Principal Component (1-50)") +
    ggplot2::ylab("Variance per PC (%)") +
    ggplot2::ggtitle("Scree plot transcriptomics") +
    ggplot2::theme(
        text = ggplot2::element_text(face = "bold",size = 8, colour = "black"),
        axis.title = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = "bold")
    ) +
    ggplot2::scale_x_continuous(
        breaks = c(0,
                   10,
                   20,
                   30,
                   40,
                   50),
        labels = c(0,
                   10,
                   20,
                   30,
                   40,
                   50)
    )

ggsave(screeplot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/scree_plot_transcriptomics.png", width = 90, height = 60, units="mm")

################################################################################################################################################
################################################      CREATE PCA PLOT      #####################################################################
################################################################################################################################################

# Create plot
pca_transcriptomics <- data_pca |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = fiber_type,
            names = fiberID
        )
    ) +
    ggplot2::geom_point(
        size = 1,
        alpha = 0.65,
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA Transcriptomics") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        )
    ) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           vjust = 1.5)
    ) +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        axis.text = ggplot2::element_text(size=8),
    ) +
    ggplot2::xlab("PC1 (11.1%)") +
    ggplot2::ylab("PC2 (3.5%)") +
    ggplot2::scale_color_manual(
        "Fiber type",
        values = c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325", "#D2631C")
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    ) +
    scale_x_continuous(trans = "reverse",
                       breaks=c(-50, 0, 50, 100),
                       labels=c("50", "0", "-50", "-100")
    )

pca_transcriptomics # Conclusion: exactly the same as seurat, only %var numbers are different



################################################################################################################################################
#########################################      CLUSTERING QUALITY CONTROL PER SUBJECT    #######################################################
################################################################################################################################################

scales::viridis_pal(option = "turbo")(length(unique(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data$subject)))

# PCA plot for paper

seurat_by_subject <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest

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

PCA_subject <- Seurat::DimPlot(seurat_by_subject,
                                label = FALSE,
                                reduction = "pca",
                                pt.size = 0.5,
                                group.by = "subject") +
    ggtitle("PCA Transcriptomics - by subject") +
    xlab("PC1 (11.1%)") +
    ylab("PC2 (3.5%)") +
    guides(color = guide_legend(override.aes = list(size=1), ncol=2)) +
    scale_color_manual(values = c("#30123BFF", "#424AB3FF", "#467EF4FF", "#31AFF5FF", "#18DAC7FF", "#38F491FF", "#83FF52FF", "#BDF534FF", "#E9D539FF", "#FEAA33FF", "#F8721CFF", "#E03F08FF", "#B61C02FF", "#7A0403FF")) +
    theme(
        text = element_text(face="bold", colour="black", size=8),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = "right",
        legend.text=element_text(size=4)
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    )

ggsave(PCA_subject,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/PCA_transcriptomics_subject.png",
       width = 90,
       height = 60,
       units="mm")

################################################################################################################################################
#########################################      CLUSTERING QUALITY CONTROL PER CONDITION    #####################################################
################################################################################################################################################

scales::viridis_pal(option = "turbo")(length(unique(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data$condition)))

# PCA plot for paper
PCA_condition <- Seurat::DimPlot(seurat_by_subject,
                               label = FALSE,
                               reduction = "pca",
                               pt.size = 0.5,
                               group.by = "condition") +
    ggtitle("PCA Transcriptomics - by day") +
    xlab("PC1 (11.1%)") +
    ylab("PC2 (3.5%)") +
    guides(color = guide_legend(override.aes = list(size=1), ncol=1)) +
    scale_color_manual(labels = c("Day A", "Day B", "Day C"), values = c("#B61C02FF", "#31AFF5FF", "#30123BFF")) +
    theme(
        text = element_text(face="bold", colour="black", size=8),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = "right",
        legend.text=element_text(size=4)
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    )

ggsave(PCA_condition,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/PCA_transcriptomics_condition.png",
       width = 90,
       height = 60,
       units="mm")

################################################################################################################################################
##############################################      CLUSTERING QUALITY CONTROL PER FIBER TYPE    ###############################################
################################################################################################################################################

# Plot the PCA - for paper

PCA_fibertype <- Seurat::DimPlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                 label = FALSE,
                                 reduction = "pca",
                                 cols = c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325", "#D2631C"),
                                 pt.size = 1,
                                 group.by = "fiber_type_MYH_hybrids") +
    guides(color = guide_legend(override.aes = list(size=1))) +
    theme_minimal() +
    ggtitle("PCA Transcriptomics") +
    xlab("PC1 (11.1%)") +
    ylab("PC2 (3.5%)") +
    theme(
        text = element_text(face="bold", colour="black", size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    )

ggsave(PCA_fibertype, filename = "8 Fiber heterogeneity (only rested samples)/1 Clustering/By fiber type slowfast/PCA fibertype.png", width = 128, height = 90, units="mm")


################################################################################################################################################
###########################################       PC DRIVERS PLOT     ####################################################################
################################################################################################################################################

# Load PC gene loading list
PC_genes <- read_csv("~/single_fiber_heterogeneity/data/transcriptomics_PC_loadings.csv")

# Rename first column
colnames(PC_genes)[1] = "Gene"

# Data wrangling
Top_PC1_pos <- PC_genes %>% arrange(desc(PC_1)) %>% dplyr::slice(1:8) %>% dplyr::select(Gene, PC_1) %>% arrange(PC_1)
Top_PC1_neg <- PC_genes %>% arrange(PC_1) %>% dplyr::slice(1:8) %>% dplyr::select(Gene, PC_1)

Top_PC2_pos <- PC_genes %>% arrange(desc(PC_2)) %>% dplyr::slice(1:8) %>% dplyr::select(Gene, PC_2)  %>% arrange(PC_2)
Top_PC2_neg <- PC_genes %>% arrange(PC_2) %>% dplyr::slice(1:8) %>% dplyr::select(Gene, PC_2)

PC_df <- data.frame(
    gene = c(Top_PC1_neg$Gene, Top_PC1_pos$Gene, Top_PC2_neg$Gene, Top_PC2_pos$Gene),
    PC_score = c(Top_PC1_neg$PC_1, Top_PC1_pos$PC_1, -Top_PC2_neg$PC_2, -Top_PC2_pos$PC_2),
    PC = c(rep("PC1", 16), (rep("PC2", 16))),
    No = c(30:15, 1:8, 16, 22, 9:14)
)

# Create plot
PC_drivers_plot <- ggplot(PC_df, aes(PC_score, fct_reorder(gene, No))) +
    geom_col(fill = "#618F70") +
    facet_grid(~ PC, scales="free") +
    theme_minimal() +
    geom_vline(xintercept = 0, colour="black") +
    ggtitle("PC drivers transcriptomics") +
    ylab("Genes") +
    ylab("PC score") +
    theme(
        text = ggplot2::element_text(face = "bold",size = 6, colour = "black"),
        plot.title = element_text(hjust = 0.5, face="bold", size = 8),
        legend.position = "none",
        strip.background =element_rect(fill= c("#DFEAE2")),
        strip.text.x = element_text(size = 6, face = "bold", hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y = element_blank()
    )

ggsave(PC_drivers_plot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/PC_drivers_transcriptomics.png", width = 60, height = 90, units="mm")
