

# Load packages (CRAN and Bioconductor)
library(tidyverse)
library(ggrepel)
library(Seurat)
library(ggpubr)
library(viridis)

################################################################################################################################################
########################################################       FIGURE 4A     ###################################################################
################################################################################################################################################

# Load result files with all genes
DEG <- read.csv(here::here("data/figure_4/slow_vs_fast_transcriptomics.csv", header = T)) %>%
    drop_na(GENEID) %>%
    drop_na(padj)

# Create Volcano dataframe
volcanodata <- DEG %>%
    arrange(padj) %>%
    dplyr::mutate("significant" = dplyr::case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "enriched in slow - LFC",
        padj < 0.05 & log2FoldChange < -1 ~ "enriched in fast - LFC",
        padj < 0.05 & log2FoldChange > 0 ~ "enriched in slow",
        padj < 0.05  & log2FoldChange < 0 ~ "enriched in fast",
        TRUE ~ "not significant"
    ))  %>%
    dplyr::mutate("value" = -log10(padj)) %>%
    dplyr::mutate("names" = dplyr::case_when(
        GENE %in% c(
            "MYH2",
            "LDHA",
            "IRX3",
            "NANOS1",
            "FAM166B",
            "SH3RF2",
            "MYH7",
            "XPO4",
            "LDHB",
            "RP11-451G4.2",
            "USP54"
        ) ~ GENE,
        TRUE ~ ""
    ))

# Create plot
volcanodata |>
    ggplot2::ggplot(ggplot2::aes(
        x = log2FoldChange,
        y = -log10(pvalue),
        color = significant,
        label = names
    )) +
    ggplot2::geom_point(
        size = 0.5
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "#c6ecc8",
            "#5DC863FF",
            "#f0b3fe",
            "#440154FF",
            "grey"
        ),
        name = ""
    ) +
    ggrepel::geom_label_repel(
        data = volcanodata |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = log2FoldChange,
            y = -log10(pvalue),
            fill = significant,
            label = names
        ),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#e5f5e0",
        "#efedf5",
        "#efedf5"
    )) +
    ggplot2::ggtitle("Slow Vs Fast Transcriptomics") +
    ggplot2::theme_classic() +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 7,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::xlab("log2FC (Slow - Fast)") +
    ggplot2::ylab("-log10(P-value)") +
    scale_x_continuous(limits = c(-6, 6))


ggplot2::ggsave(here::here("doc/figures/figure_4/figure_4A.png"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
########################################################       FIGURE 4A - inlay     ###################################################################
################################################################################################################################################


# Load filtered Seurat object ---------------------------------------------
load(here::here("data/figure_4/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata"))

# Assign identity of resolution 0.7 to clusters ---------------------------------------------
Idents(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"

filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data$final_cluster

# Assign slow vs fast clustering
metadata <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data |>
    dplyr::mutate(fiber_type_seurat = dplyr::case_when(
        final_cluster == "Slow1" ~ "slow",
        final_cluster == "Slow2" ~ "slow",
        final_cluster == "Fast1" ~ "fast",
        final_cluster == "Fast2" ~ "fast",
        final_cluster == "Fast3" ~ "fast",
        final_cluster == "Intermediate1" ~ "hybrid",
        TRUE ~ "NA"
    ))

filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest[["umap"]]@cell.embeddings |>
    as.data.frame() |>
    #dplyr::mutate(UMAP_1 = UMAP_1 * -1) |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = metadata$fiber_type_seurat)
    ) +
    ggplot2::geom_point(size = 0.025,
                        alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(
        "#5DC863FF",
        "grey",
        "#440154FF"
    )) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        text = ggplot2::element_blank(),
        legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
    )

ggplot2::ggsave(here::here("doc/figures/figure_4/figure_4A_inlay.png"),
                units = "mm",
                height = 12.5,
                width = 12.5)

################################################################################################################################################
########################################################       FIGURE 4B     ###################################################################
################################################################################################################################################

library(eulerr)
library(tidyverse)
library(ggplot2)

# Loading seurat clusters and merging fast and slow -----------------------

metadata <-
    vroom::vroom(here::here("data/metadata_proteomics_seurat_clusters.csv")) |>
    dplyr::rename(
        "fiber_type_seurat" = "seurat_clusters"
    )

proteomics_data <-vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    log2() |>
    as.data.frame()

pseudobulk_maker <- function(.data, metadata, subject_id, grouping, colname) {
    selection_vector <- metadata |>
        dplyr::filter(fiber_type_seurat == grouping) |>
        dplyr::filter(subject == subject_id) |>
        dplyr::pull("fiberID")

    pseudobulk_median_calculation <- .data |>
        as.data.frame() |>
        dplyr::select(selection_vector) |>
        t() |>
        as.data.frame() |>
        dplyr::mutate(dplyr::across(.cols = everything(), median, na.rm = T)) |>
        unique() |>
        t()

    colnames(pseudobulk_median_calculation) <- colname
    return(pseudobulk_median_calculation)
}

data_pseudobulk <- data.frame(
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P1",
        grouping = "slow",
        colname = "P1_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P2",
        grouping = "slow",
        colname = "P2_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P3",
        grouping = "slow",
        colname = "P3_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P4",
        grouping = "slow",
        colname = "P4_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P5",
        grouping = "slow",
        colname = "P5_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P1",
        grouping = "fast",
        colname = "P1_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P2",
        grouping = "fast",
        colname = "P2_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P3",
        grouping = "fast",
        colname = "P3_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P4",
        grouping = "fast",
        colname = "P4_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P5",
        grouping = "fast",
        colname = "P5_fast"
    )
)

data_grouping_pseudobulk <- data.frame(
    "sample_ID" = c(
        "P1_slow",
        "P2_slow",
        "P3_slow",
        "P4_slow",
        "P5_slow",
        "P1_fast",
        "P2_fast",
        "P3_fast",
        "P4_fast",
        "P5_fast"
    ),
    "fiber_type" = c(
        "slow",
        "slow",
        "slow",
        "slow",
        "slow",
        "fast",
        "fast",
        "fast",
        "fast",
        "fast"
    )
)

data_limma <- data_pseudobulk |>
    limma::normalizeBetweenArrays(method = "quantile")

vector_factor_fiber_type <- factor(data_grouping_pseudobulk$fiber_type,
                                   levels = c(
                                       "slow",
                                       "fast"
                                   )
)

design_matrix <-
    model.matrix(
        ~ 0 + vector_factor_fiber_type,
        data_grouping_pseudobulk
    )

colnames(design_matrix) <- c(
    "slow",
    "fast"
)

fit <- limma::lmFit(
    data_limma,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "slow_vs_fast" = slow - fast,
    levels = design_matrix
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp <- limma::eBayes(tmp)

DE_results <- limma::topTable(tmp,
                              sort.by = "P",
                              n = Inf
)

data_volcano <- DE_results |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "enriched in slow",
        adj.P.Val < 0.05 & logFC < 0 ~ "enriched in fast",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes")

data_volcano <- limma::topTable(tmp,
                                sort.by = "P",
                                n = Inf
) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val < 0.05 & logFC > 1 ~ "highly upregulated in slow",
        adj.P.Val < 0.05 & logFC > 0 & logFC < 2 ~ "upregulated in slow",
        adj.P.Val < 0.05 & logFC < -1 ~ "highly upregulated in fast",
        adj.P.Val < 0.05 & logFC < 0 & logFC > -2 ~ "upregulated in fast",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::mutate(labelling = dplyr::case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "slow",
        adj.P.Val < 0.05 & logFC < 0 ~ "fast",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("value" = -log10(P.Value)) |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "MYH2",
            "USP28",
            "ZC3H8",
            "USP48",
            "TNNC2",
            "LDHA",
            "MYH7",
            "SUB1",
            "TNNT1",
            "DPYSL3",
            "TPM3",
            "LDHB"
        ) ~ Genes,
        TRUE ~ ""
    ))

data_volcano$signifficant <- factor(data_volcano$signifficant,
                                    levels = c("highly upregulated in slow",
                                               "upregulated in slow",
                                               "not signifficant",
                                               "upregulated in fast",
                                               "highly upregulated in fast"))

data_volcano |>
    ggplot2::ggplot(ggplot2::aes(
        x = logFC,
        y = -log10(P.Value),
        color = signifficant,
        label = names
    )) +
    ggplot2::geom_point(
        size = 0.5,
        alpha = ifelse(data_volcano$signifficant == "not signifficant", 1, 1)
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "#440154FF",
            "#f0b3fe",
            "grey",
            "#c6ecc8",
            "#5DC863FF"
        ),
        name = ""
    ) +
    ggrepel::geom_label_repel(
        data = data_volcano |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = labelling,
            label = names
        ),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#e5f5e0",
        "#efedf5"
    )) +
    ggplot2::ggtitle("Slow Vs Fast Proteomics") +
    ggplot2::theme_classic() +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 12,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::xlab("log2FC (Slow - Fast)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7),
        legend.position = "none"
    )

ggplot2::ggsave(
    here::here("doc/figures/figure_4/figure_4_B.png"),
    device = "png",
    height = 60,
    width = 60,
    units = "mm"
)

# write.csv(data_pseudobulk,
#           here::here("data/data_proteomics_pseudobulk.csv"))

# write.csv(DE_results,
#           file = here::here("data/DE_analysis_slow_vs_fast_proteomics.csv")
# )

# Inlay plot: -------------------------------------------------------------

data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

metadata <- metadata |>
    tibble::column_to_rownames("...1")

# Seurat workflow ---------------------------------------------------------

seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                              meta.data = metadata)

seurat_proteome[["RNA"]]$data <- seurat_proteome[["RNA"]]$counts

################################################################################################################################################
########################################################      PCA   ############################################################################
################################################################################################################################################

# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                                selection.method = "vst")

# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)

#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,
                                  features = Seurat::VariableFeatures(object = seurat_proteome))

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))

################################################################################################################################################
################################################      DIMENSIONALITY REDUCTION   ##############################################################
################################################################################################################################################

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6)

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 1,
    pt.size = 0.05,
    group.by = "fiber_type_seurat"
)

seurat_proteome[["umap"]]@cell.embeddings |>
    as.data.frame() |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = metadata$fiber_type_seurat)
    ) +
    ggplot2::geom_point(size = 0.025,
                        alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(
        "#5DC863FF",
        "#440154FF"
    )) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        text = ggplot2::element_blank(),
        legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
    )

ggplot2::ggsave(here::here("doc/figures/figure_4/figure_4B_inlay.png"),
                units = "mm",
                height = 12.5,
                width = 12.5)

################################################################################################################################################
#############################################       FIGURE 4C - TRANSCRIPTOMICS    #############################################################
################################################################################################################################################


feature_plot_IRX3 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("IRX3"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_USP54 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                          features = c("USP54"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_LDHA <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("LDHA"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_LDHB <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("LDHB"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

featureplots <- ggpubr::ggarrange(
    feature_plot_LDHA,
    feature_plot_LDHB,
    feature_plot_IRX3,
    feature_plot_USP54,
    ncol = 2,
    nrow = 2) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

annotate_figure(featureplots, top = text_grob("Transcriptomics",
                                                              color = "black", face = "bold", size = 7))

ggsave(here::here("doc/figures/figure_4/figure_4C_transcriptomics.png"),
       width = 60,
       height = 45,
       units="mm")


# Feature plots proteomics: -----------------------------------------------

feature_plot_USP28 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("USP28"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_USP28 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("USP28"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_DPYSL3 <- Seurat::FeaturePlot(seurat_proteome,
                                           features = c("DPYSL3"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

featureplots <- ggpubr::ggarrange(feature_plot_USP28,
                                  feature_plot_DPYSL3,
                                  ncol = 2,
                                  nrow = 1) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

featureplots <- annotate_figure(featureplots, top = text_grob("Proteomics",
                                                              color = "black", face = "bold", size = 7))

ggsave(featureplots,
       filename = "doc/figures/figure_4/figure_4C_proteomics.png",
       width = 60,
       height = 25,
       units="mm")

################################################################################################################################################
####################################       FIGURE 4D    ##########################################################
################################################################################################################################################

ORA_slow_overlap_simplify_df <- read.csv(here::here("data/figure_4/Common_slow.csv"))
ORA_fast_overlap_simplify_df <- read.csv(here::here("data/figure_4/Common_fast.csv"))
ORA_slow_transcriptome_simplify_df <- read.csv(here::here("data/figure_4/Transcriptome-only_slow.csv"))
ORA_fast_transcriptome_simplify_df <- read.csv(here::here("data/figure_4/Transcriptome-only_fast.csv"))
ORA_slow_proteome_simplify_df <- read.csv(here::here("data/figure_4/Proteome-only_slow.csv"))
ORA_fast_proteome_simplify_df <- read.csv(here::here("data/figure_4/Proteome-only_fast.csv"))


# Select categories for Slow Common
slow_overlap <- ORA_slow_overlap_simplify_df %>%
    dplyr::filter(
        Description == "fatty acid beta-oxidation" |
            Description == "mitochondrial matrix" |
            Description == "oxidoreductase activity" |
            Description == "regulation of cellular ketone metabolic process" |
            Description == "myofilament"

    )

# Select categories for Fast Common
fast_overlap <- ORA_fast_overlap_simplify_df %>%
    dplyr::filter(
        Description == "carbohydrate catabolic process" |
            Description == "muscle contraction" |
            Description == "calcium ion transport" |
            Description == "protein serine/threonine phosphatase complex"
    )


# Select categories for Slow Transcriptomics
slow_transcriptome <- ORA_slow_transcriptome_simplify_df %>%
    dplyr::filter(
        Description == "extracellular matrix" |
            Description == "receptor ligand activity" |
            Description == "amide binding"
    )

# Select categories for Slow Proteomics
slow_proteome <- ORA_slow_proteome_simplify_df %>%
    dplyr::filter(
        Description == "mitochondrial translation" |
            Description == "NADH dehydrogenase complex" |
            Description == "branched-chain amino acid metabolic process"
    )

# Select categories for Fast Transcriptomics
fast_transcriptome <- ORA_fast_transcriptome_simplify_df %>%
    dplyr::filter(
        Description == "cell fate commitment" |
            Description == "regulation of transporter activity" |
            Description == "ligand-activated transcription factor activity"
    )

# Select categories for Fast Proteomics
fast_proteome <- ORA_fast_proteome_simplify_df %>%
    dplyr::filter(
        Description == "P-body"
    )


# Add required columns for plotting
slow_overlap$comparison <- rep("Slow - overlap", nrow(slow_overlap))
slow_overlap$data <- rep("Overlap", nrow(slow_overlap))

fast_overlap$comparison <- rep("Fast - overlap", nrow(fast_overlap))
fast_overlap$data <- rep("Overlap", nrow(fast_overlap))

slow_transcriptome$comparison <- rep("Slow - Transcriptome", nrow(slow_transcriptome))
slow_transcriptome$data <- rep("Transcriptome", nrow(slow_transcriptome))

slow_proteome$comparison <- rep("Slow - Proteome", nrow(slow_proteome))
slow_proteome$data <- rep("Proteome", nrow(slow_proteome))

fast_transcriptome$comparison <- rep("Fast - Transcriptome", nrow(fast_transcriptome))
fast_transcriptome$data <- rep("Transcriptome", nrow(fast_transcriptome))

fast_proteome$comparison <- rep("Fast - Proteome", nrow(fast_proteome))
fast_proteome$data <- rep("Proteome", nrow(fast_proteome))

# Combine all into one dataframe
data_plot <- bind_rows(slow_overlap, slow_transcriptome, slow_proteome, fast_overlap, fast_transcriptome, fast_proteome)
data_plot$comparison = factor(data_plot$comparison, levels=c(
    "Slow - overlap",
    "Slow - Transcriptome",
    "Slow - Proteome",
    "Fast - overlap",
    "Fast - Transcriptome",
    "Fast - Proteome"
))
data_plot$data = factor(data_plot$data, levels=c("Overlap", "Transcriptome", "Proteome"))
data_plot$graphdir <- rep(19:1)


# Create dot plot -  long format for Nat Comms revision

ggplot(data_plot, aes(y = fct_reorder(Description, graphdir), x = comparison)) +
    geom_point(aes(color= comparison, alpha=foldEnrich, size = -log10(p.adjust))) +
    scale_color_manual(values = c("#440154FF", "#440154FF", "#440154FF", "#5DC863FF", "#5DC863FF", "#5DC863FF")) +
    scale_alpha_continuous(
        name = "foldEnrich",
        limits = c(2.05, 8.55),
        breaks = c(2.05, 8.55),
        range = c(0.2, 1)
    ) +
    scale_size_continuous(
        name = "-log10(p.adjust)",
        limits = c(1.32, 12.29),
        breaks = c(1.32, 12.29),
        range = c(1.5, 5)
    ) +
    #ggtitle("Enrichment by fiber type") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=4),
        axis.text.y = element_text(size=5),
        axis.title = element_blank(),
        text = element_text(face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, face="bold", vjust=0, size=8),
        legend.position = "none",
        panel.grid.major.y = element_blank() ,
        #panel.grid.major.x = element_line( linewidth=.1, color="grey" )
    )

ggsave(here::here("doc/figures/figure_4/figure_4D.png"),
       width = 75,
       height = 85,
       units="mm")


################################################################################################################################################
####################################       FIGURE 4E    ##########################################################
################################################################################################################################################

# Load prioritization list
prioritization_fast <- readxl::read_excel(here::here("data/figure_4/priortization_df_COARSE.xlsx")) %>%
    dplyr::filter(cluster == "fast")

prioritization_slow <- readxl::read_excel(here::here("data/figure_4/priortization_df_COARSE.xlsx")) %>%
    dplyr::filter(cluster == "slow")

# Load DE lists
DE_transcriptomics_fast <- read_csv(here::here("data/figure_5/slow_vs_fast.csv")) %>%
    dplyr::rename(TF = GENE) %>%
    dplyr::rename(log2FoldChange_pseudobulk = log2FoldChange) %>%
    dplyr::rename(padj_pseudobulk = padj) %>%
    dplyr::filter(log2FoldChange_pseudobulk < 0) %>%  # Filter for only genes enriched in fast fibers
    dplyr::select(TF, log2FoldChange_pseudobulk, padj_pseudobulk) %>%
    dplyr::filter(TF %in% prioritization_fast$TF) %>%
    mutate(log2FoldChange_pseudobulk = -log2FoldChange_pseudobulk) # 114 potential TFs for fast fibers

DE_transcriptomics_slow <- read_csv(here::here("data/figure_5/slow_vs_fast.csv")) %>%
    dplyr::rename(TF = GENE) %>%
    dplyr::rename(log2FoldChange_pseudobulk = log2FoldChange) %>%
    dplyr::rename(padj_pseudobulk = padj) %>%
    dplyr::filter(log2FoldChange_pseudobulk > 0) %>%  # Filter for only genes enriched in slow fibers
    dplyr::select(TF, log2FoldChange_pseudobulk, padj_pseudobulk) %>%
    dplyr::filter(TF %in% prioritization_slow$TF) # 66 potential TFs for fast fibers

DE_proteomics <- read_csv(here::here("data/figure_4/DE_analysis_slow_vs_fast_proteomics.csv")) %>%
    dplyr::rename(TF = '...1') %>%
    dplyr::rename(log2FoldChange_pseudobulk = logFC) %>%
    dplyr::select(TF, log2FoldChange_pseudobulk) %>%
    dplyr::filter(TF %in% prioritization_fast$TF) # Only 9 TFs with low LFC, no follow-up

# Add pseudobulk DE score
prioritization_fast <- prioritization_fast %>%
    dplyr::filter(TF %in% DE_transcriptomics_fast$TF) %>%
    left_join(DE_transcriptomics_fast, by = "TF")

prioritization_slow <- prioritization_slow %>%
    dplyr::filter(TF %in% DE_transcriptomics_slow$TF) %>%
    left_join(DE_transcriptomics_slow, by = "TF")

# Check correlation DE score seurat vs LFC pseudobulk
prioritization_fast %>%
    ggplot(aes(DE_score, log2FoldChange_pseudobulk)) +
    geom_point()

prioritization_slow %>%
    ggplot(aes(DE_score, log2FoldChange_pseudobulk)) +
    geom_point()

# Create df with scaled scores = calculate prioritization score
priortization_fast_scaledAll = prioritization_fast %>%
    mutate(DE_score_scaled = nichenetr::scaling_zscore(log2FoldChange_pseudobulk), RSS_score_scaled = nichenetr::scaling_zscore(RSS)) %>%
    mutate(prioritization_score = DE_score_scaled + RSS_score_scaled) %>%
    arrange(-prioritization_score)

priortization_slow_scaledAll = prioritization_slow %>%
    mutate(DE_score_scaled = nichenetr::scaling_zscore(log2FoldChange_pseudobulk), RSS_score_scaled = nichenetr::scaling_zscore(RSS)) %>%
    mutate(prioritization_score = DE_score_scaled + RSS_score_scaled) %>%
    arrange(-prioritization_score)

priortization_fast_scaledAll <- priortization_fast_scaledAll %>%
    dplyr::mutate("names" = dplyr::case_when(
        TF %in% c(
            "MAFA",
            "PITX1",
            "EGR1",
            "MAF",
            "MYF6"
        ) ~ TF,
        TRUE ~ ""
    ))

priortization_fast_scaledAll %>%
    ggplot(aes(RSS, log2FoldChange_pseudobulk, fill = prioritization_score, label = names)) +
    geom_point(colour = "black", shape = 21, size = 1) +
    scale_fill_viridis_c("Prioritization score", option = "magma") +
    theme_bw() +
    ggtitle("Prioritized TFs - Fast fibers") +
    ggrepel::geom_label_repel(
        data = priortization_fast_scaledAll %>% dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = RSS,
            y = log2FoldChange_pseudobulk,
            #fill = prioritization_score,
            label = names
        ),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ylab("L2FC vs slow fibers") +
    xlab("Regulon specificity score") +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
        text = element_text(size = 6, face = "bold"),
        legend.key.size = unit(2.5, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(1,0,5,0),
        legend.box.margin=margin(-10,-10,-10,-10)
    )

ggsave(here::here("doc/figures/figure_4/figure_4E_part1.png"),
       units = "mm",
       height = 41,
       width = 75)

priortization_slow_scaledAll <- priortization_slow_scaledAll %>%
    dplyr::mutate("names" = dplyr::case_when(
        TF %in% c(
            "EPAS1",
            "ZSCAN30",
            "FOSL2",
            "ETS1",
            "KLF15",
            "RFXANK",
            "ATF5",
            "HOXA7"
        ) ~ TF,
        TRUE ~ ""
    ))

priortization_slow_scaledAll %>%
    ggplot(aes(RSS, log2FoldChange_pseudobulk, fill = prioritization_score, label = names)) +
    geom_point(colour = "black", shape = 21, size = 1) +
    scale_fill_viridis_c("Prioritization score", option = "magma") +
    theme_bw() +
    ggtitle("Prioritized TFs - Slow fibers") +
    ggrepel::geom_label_repel(
        data = priortization_slow_scaledAll %>% dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = RSS,
            y = log2FoldChange_pseudobulk,
            #fill = prioritization_score,
            label = names
        ),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ylab("L2FC vs fast fibers") +
    xlab("Regulon specificity score") +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
        text = element_text(size = 6, face = "bold"),
        legend.key.size = unit(2.5, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(1,0,5,0),
        legend.box.margin=margin(-10,-10,-10,-10)
    )

ggsave(here::here("doc/figures/figure_4/figure_4E_part2.png"),
       units = "mm",
       height = 41,
       width = 75)

################################################################################################################################################
####################################       FIGURE 4F    ##########################################################
################################################################################################################################################

MAFA <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                            features = c("MAFA"),
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")


EPAS1 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                             features = c("EPAS1"),
                             pt.size = 0.1,
                             order = F) +
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
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")


featureplots_TF <- ggpubr::ggarrange(
    MAFA,
    EPAS1,
    ncol = 1,
    nrow = 2) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggpubr::annotate_figure(featureplots_TF, top = ggpubr::text_grob("TF expression",
                                                                                    color = "black", face = "bold", size = 7))

ggsave(here::here("doc/figures/figure_4/figure_4F.png"),
       width = 50,
       height = 70,
       units="mm")






