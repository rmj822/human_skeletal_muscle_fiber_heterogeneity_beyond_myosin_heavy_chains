
library(Seurat)
library(tidyverse)
library(patchwork)


################################################################################################################################################
################################################       PROTEOMICS DATA WRANGLING       ####@####################################################
################################################################################################################################################

data_proteomics <-read.csv(here::here("data/figure_2/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "X") |>
    tibble::column_to_rownames("Gene.name")

metadata <-read.csv(here::here("data/figure_2/metadata_proteomics.csv"))|>
    dplyr::rename("Sample" = "X") |>
    tibble::column_to_rownames("Sample")

pca_object <- prcomp(t(data_proteomics), scale = T)

all(colnames(data_proteomics) %in% rownames(metadata))
all(colnames(data_proteomics) == rownames(metadata))

data_pca_P <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    tibble::add_column(
        fiber_type = metadata$fiber_type,
        fiberID = rownames(metadata),
        digestion_batch = as.factor(metadata$digestion_batch),
        subject = as.factor(metadata$subject),
        date = metadata$date_isolation,
        length = metadata$fiber_length
    )

data_pca_P$fiber_type <- factor(data_pca_P$fiber_type,
                              levels = c("Hybrid 1/2A",
                                         "Hybrid 2A/2X",
                                         "Type 2A",
                                         "Type 1"))


################################################################################################################################################
################################################       TRANSCRIPTOMICS DATA WRANGLING       ####@####################################################
################################################################################################################################################

load(here::here("data/figure_2/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata"))

factor(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids,
       levels = c("Hybrid 1/2A", "Type 1", "Type 2A", "Type 2X","Hybrid 2A/2X"))

filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids <-
    factor(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids,
           levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))

Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "fiber_type_MYH_hybrids"

my_cols <- c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325", "#D2631C")


# Extract SCT scaled data
transcriptomics_scaledata <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "scale.data")

# Extract SCT normalized but not log-transformed counts
transcriptomics_counts <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")

# Set assay to SCT
DefaultAssay(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "SCT"

pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2)

# Load metadata
metadata <- filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data

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
################################################      FIGURE 2SA      #####################################################################
################################################################################################################################################

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

Seurat::DimPlot(seurat_by_subject,
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

ggplot2::ggsave(here::here("doc/figures/figure_2_S2/figure_2_S2A.png"),
                units = "mm",
                height = 60,
                width = 90)

################################################################################################################################################
#########################################      FIGURE 2 S2B    #####################################################
################################################################################################################################################

Seurat::DimPlot(seurat_by_subject,
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

ggplot2::ggsave(here::here("doc/figures/figure_2_S2/figure_2_S2B.png"),
                units = "mm",
                height = 60,
                width = 90)

################################################################################################################################################
################################################       FIGURE 2 S2C       ####@####################################################
################################################################################################################################################

# PCA by subject proteomics -----------------------------------------------
data_pca_P <- data_pca_P %>%
    dplyr::mutate(subject = dplyr::case_when(
        subject == "FOR2" ~ "P1",
        subject == "FOR4" ~ "P2",
        subject == "FOR9" ~ "P3",
        subject == "FOR10" ~ "P4",
        subject == "FOR11" ~ "P5",
        TRUE ~ "unknown"
    ))

data_pca_P$subject <- factor(data_pca_P$subject, levels = c("P1", "P2", "P3", "P4", "P5"))

data_pca_P |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = subject,
            names = fiberID
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        alpha = 1
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("PCA Proteomics by subject") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        )
    ) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::theme(
        legend.position = "right",
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        axis.text = ggplot2::element_text(size=8),
        legend.text=element_text(size=4)
    ) +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    guides(color = guide_legend(override.aes = list(size=1), ncol=1)) +
    scale_color_manual(values = c("#786A7F",
                                  "#8FCAE8",
                                  "#C8D685",
                                  "#EAAE7C",
                                  "#9F6461"),
                       name = "")


ggplot2::ggsave(here::here("doc/figures/figure_2_S2/figure_2_S2C.png"),
                units = "mm",
                height = 60,
                width = 90)


################################################################################################################################################
###########################################      FIGURE 2 S2D   ####################################################################
################################################################################################################################################

# Load PC gene loading list
PC_genes <- read_csv(here::here("data/figure_2/transcriptomics_PC_loadings.csv"))

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
ggplot(PC_df, aes(PC_score, fct_reorder(gene, No))) +
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

ggplot2::ggsave(here::here("doc/figures/figure_2_S2/figure_2_S2D.png"),
                units = "mm",
                height = 60,
                width = 90)


################################################################################################################################################
################################################       FIGURE 2 S2E       ####@####################################################
################################################################################################################################################

# PCA_drivers -------------------------------------------------------------

PC_genes <- pca_object$rotation |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    dplyr::rename("PC_1" = PC1,
                  "PC_2" = PC2) |>
    tibble::rownames_to_column("Gene")

# Data wrangling
Top_PC1_pos <- PC_genes |>
    dplyr::arrange(desc(PC_1)) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_1) |>
    dplyr::arrange(PC_1)

Top_PC1_neg <- PC_genes  |>
    dplyr::arrange(PC_1) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_1)

Top_PC2_pos <- PC_genes |>
    dplyr::arrange(desc(PC_2)) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_2)  |>
    dplyr::arrange(PC_2)

Top_PC2_neg <- PC_genes |>
    dplyr::arrange(PC_2) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_2)

PC_df <- data.frame(
    gene = c(Top_PC1_neg$Gene,
             Top_PC1_pos$Gene,
             Top_PC2_neg$Gene,
             Top_PC2_pos$Gene),
    PC_score = c(Top_PC1_neg$PC_1,
                 Top_PC1_pos$PC_1,
                 Top_PC2_neg$PC_2,
                 Top_PC2_pos$PC_2),
    PC = c(rep("PC1", 16),
           (rep("PC2", 16))))

PC_df <- PC_df %>%
    dplyr::mutate(PC_score = -PC_score) %>%
    dplyr::mutate(order = c(c(16:1), c(32:17)))

# Create plot
ggplot2::ggplot(PC_df,
                ggplot2::aes(x = PC_score,
                            y = forcats::fct_reorder(gene,
                                                     order,
                                                     .desc = TRUE))) +
    ggplot2::geom_col(fill = "#045a8d") +
    ggplot2::facet_grid(~ PC, scales="free") +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = 0, colour="black") +
    ggplot2::ggtitle("PC drivers proteomics") +
    ggplot2::ylab("Genes") +
    ggplot2::ylab("PC score") +
    ggplot2::theme(
        text = ggplot2::element_text(face = "bold",size = 7, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face="bold"),
        legend.position = "none",
        strip.background = ggplot2::element_rect(fill= c("#9ecae1")),
        strip.text.x = ggplot2::element_text(size = 6, face = "bold", hjust=0.5),
        axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
        axis.title.y = ggplot2::element_blank()
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
    )

ggplot2::ggsave(here::here("doc/figures/figure_2_S2/figure_2_S2E.png"),
                units = "mm",
                height = 60,
                width = 90)
