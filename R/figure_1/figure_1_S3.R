################################################################################################################################################
########################################################      Panel A   ############################################################################
################################################################################################################################################

load(here::here("data/figure_2/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata"))

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

ggsave(screeplot, filename = "doc/figures/figure_1_S3/figure_1_S3A.png", width = 90, height = 60, units="mm")

################################################################################################################################################
########################################################      Panel B   ############################################################################
################################################################################################################################################

# Load proteomics data and create Seurat object ---------------------------

data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

pca_object <- prcomp(t(data_proteomics), scale = TRUE)


# Extract %s and create scree plot ----------------------------------------

variance <- factoextra::get_eig(pca_object) |>
    as.data.frame() |>
    dplyr::select(variance.percent) |>
    dplyr::mutate(rank = 1:length(variance.percent)) |>
    dplyr::slice_head(n = 50)

# Nice Elbow plot for paper

ElbowPlot <- ggplot2::ggplot(variance,
                             ggplot2::aes(x = rank,
                                          y = variance.percent,
                                          color = rank > 6)) +
    ggplot2::annotate("rect",
                      xmin=-Inf,
                      xmax=6.5,
                      ymin=-Inf,
                      ymax=Inf,
                      alpha=0.2,
                      fill= "#045a8d") +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c("#045a8d", "#9ecae1")) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Principal Component (1-50)") +
    ggplot2::ylab("Variance per PC (%)") +
    ggplot2::ggtitle("Scree plot proteomics") +
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

ElbowPlot

ggplot2::ggsave(here::here("doc/figures/figure_1_S3/figure_1_S3B.png"),
                units = "mm",
                height = 60,
                width = 90)



################################################################################################################################################
########################################################      Panel C   ############################################################################
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


ggsave(UMAP_subject, filename = "doc/figures/figure_1_S3/figure_1_S3C.png", width = 90, height = 60, units="mm")


################################################################################################################################################
########################################################      Panel D   ############################################################################
################################################################################################################################################

# Filtered and processed data
data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

metadata <- vroom::vroom(
    here::here("data/metadata_proteomics_fiber_type.csv")
) |>
    dplyr::select(!1)

seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                              meta.data = metadata)


# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                                selection.method = "vst")

# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)

#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,  features = Seurat::VariableFeatures(object = seurat_proteome))

# Graph-based clustering using K-nearest neighbor graph  ---------------------------------

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6, seed.use = 42)

my_cols <- viridisLite::turbo(n = 5)

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 5,
    pt.size = 0.5,
    cols = ggplot2::alpha(my_cols, 0.65),
    group.by = "subject") +
    ggplot2::ggtitle("Resolution 0.4") +
    # scale_color_viridis_d(option = "magma")
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle("UMAP Proteomics - by subject") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        axis.text = ggplot2::element_text(size=6),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.text= ggplot2::element_text(size=4),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position = "right"
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_1_S3/figure_1_S3D.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )


################################################################################################################################################
########################################################      Panel E   ############################################################################
################################################################################################################################################

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

ggsave(UMAP_condition, filename = "doc/figures/figure_1_S3/figure_1_S3E.png", width = 90, height = 60, units="mm")
