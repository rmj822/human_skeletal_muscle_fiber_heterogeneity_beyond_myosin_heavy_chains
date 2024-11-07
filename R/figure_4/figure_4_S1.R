################################################################################################################################################
########################################################       FIGURE 4 S1A     ###################################################################
################################################################################################################################################

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

# LDH Plots

feature_plot_LDHA <- Seurat::FeaturePlot(seurat_proteome,
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

feature_plot_LDHB <- Seurat::FeaturePlot(seurat_proteome,
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

featureplots_LDH <- ggpubr::ggarrange(feature_plot_LDHA,
                                      feature_plot_LDHB,
                                      ncol = 2,
                                      nrow = 1) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

featureplots_LDH <- annotate_figure(featureplots_LDH, top = text_grob("Proteomics",
                                                                      color = "black", face = "bold", size = 7))

ggsave(featureplots_LDH,
       filename = "doc/figures/figure_4_S1/figure_4_S1.png",
       width = 60,
       height = 25,
       units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1B     ###################################################################
################################################################################################################################################

DE_proteins <- vroom::vroom(
    here::here("data/figure_4/DE_analysis_slow_vs_fast_proteomics.csv")
) |>
    dplyr::rename("Genes" = 1)

DE_genes <- vroom::vroom(
    here::here("data/figure_4/slow_vs_fast_transcriptomics.csv")
)


# Extract features per category
shared_genes <- DE_proteins |>
    dplyr::select("Genes") |>
    dplyr::bind_rows(DE_genes |>
                         dplyr::select("GENE") |>
                         dplyr::rename("Genes" = "GENE")) |>
    dplyr::filter(duplicated(Genes)) |>
    dplyr::pull(Genes)


slow_proteins <- DE_proteins |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::filter(logFC > 0) |>
    dplyr::filter(Genes %in% shared_genes) |>
    dplyr::pull(unique(Genes))

fast_proteins <- DE_proteins |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::filter(logFC < 0) |>
    dplyr::filter(Genes %in% shared_genes) |>
    dplyr::pull(unique(Genes))

slow_genes <- DE_genes |>
    dplyr::filter(padj < 0.05) |>
    dplyr::filter(log2FoldChange > 0) |>
    dplyr::filter(GENE %in% shared_genes) |>
    dplyr::pull(unique(GENE))

fast_genes <- DE_genes |>
    dplyr::filter(padj < 0.05) |>
    dplyr::filter(log2FoldChange < 0) |>
    dplyr::filter(GENE %in% shared_genes) |>
    dplyr::pull(unique(GENE))

# Generate list of genes per category
euler_object <- list("Slow proteins" = slow_proteins,
                     "Fast proteins" = fast_proteins,
                     "Slow genes" = slow_genes,
                     "Fast genes" = fast_genes)

euler_object <- lapply(euler_object, as.character)

# Make plot
euler_diagram <- ggplotify::as.ggplot(
    plot(
        euler(euler_object),
        labels = list(cex = 0.75),
        quantities = list(
            font = 1,
            cex = 0.75
        ),
        fill = c(
            ggplot2::alpha("#440154FF", 0.85),
            ggplot2::alpha("#5DC863FF", 0.85),
            ggplot2::alpha("#440154FF", 0.5),
            ggplot2::alpha("#5DC863FF", 0.5)
        ),
        alpha = 0.65,
        edges = list(
            col = c(
                ggplot2::alpha("#440154FF", 1),
                ggplot2::alpha("#5DC863FF", 1),
                ggplot2::alpha("#440154FF", 1),
                ggplot2::alpha("#5DC863FF", 1)
            ),
            lex = 2
        ),
    )
)

euler_diagram

ggplot2::ggsave(
    here::here("doc/figures/figure_4_S1/figure_1_S1B.png"),
    device = "png",
    height = 60,
    width = 90,
    units = "mm"
)

################################################################################################################################################
########################################################       FIGURE 4 S1C     ###################################################################
################################################################################################################################################

data_DE_proteomics <- DE_proteins |>
    dplyr::rename(
        "Gene.name" = 1
    ) |>
    dplyr::select(
        c(Gene.name, logFC)
    )

data_DE_transcriptomics <- DE_genes |>
    dplyr::select(
        c(GENE, log2FoldChange)
    )

# filtering genes in both matrices so they overlap:
combined_data <- data_DE_transcriptomics |>
    dplyr::filter(GENE %in% data_DE_proteomics$Gene.name) |>
    dplyr::rename(
        "Gene.name" = "GENE"
    ) |>
    dplyr::inner_join(data_DE_proteomics) |>
    dplyr::rename(
        "lfc_transcriptomics" = "log2FoldChange",
        "lfc_proteomics" = "logFC"
    ) |>
    dplyr::filter(!is.na(lfc_proteomics))

pearson_correlation <- cor.test(
    x = combined_data$lfc_transcriptomics,
    y = combined_data$lfc_proteomics
)

r_pearson <- pearson_correlation$estimate

message <- paste("r =", as.character(round(r_pearson, digits = 3)))

data_plot <- combined_data |>
    dplyr::mutate(
        color = dplyr::case_when(
            lfc_transcriptomics > 1 & lfc_proteomics > 1 ~ "high in both fast",
            lfc_transcriptomics < -1 & lfc_proteomics < -1 ~ "high in both slow",
            TRUE ~ ""
        )
    ) |>
    dplyr::mutate(
        label = dplyr::case_when(
            Gene.name %in% c("MYH2",
                             "TNNI2",
                             "TNNT1",
                             "MYH7",
                             "USP28",
                             "USP48",
                             "AKAP13",
                             "GOLGA4",
                             "MYH7B",
                             "TPM1") ~ Gene.name,
            TRUE ~ ""
        )
    )

data_plot |>
    ggplot2::ggplot(ggplot2::aes(
        x = lfc_transcriptomics,
        y = lfc_proteomics,
        label = label,
        color = color
    )) +
    ggplot2::geom_point(alpha = 0.65, size = 0.5) +
    ggplot2::scale_color_manual(values = c("gray", "#440154", "#5DC863FF")) +
    ggrepel::geom_label_repel(
        data = data_plot |>
            dplyr::filter(!label == ""),
        mapping = ggplot2::aes(
            x = lfc_transcriptomics,
            y = lfc_proteomics,
            fill = color,
            label = label),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#cccccc",
        "#efedf5",
        "#e5f5e0"
    )) +
    ggplot2::annotate(
        "text",
        x = -4.5,
        y = 5,
        label = message,
        colour = "black",
        size = 2.5
    ) +
    ggplot2::annotate(
        "text",
        x = -4.5,
        y = 4,
        label = "p < 0.001",
        colour = "black",
        size = 2.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 8, face = "bold")
    ) +
    ggplot2::xlab("slow - fast logFC transcriptomics") +
    ggplot2::ylab("slow - fast logFC proteomics")

ggplot2::ggsave(
    filename = here::here(
        "doc/figures/figure_4_S1/figure_4_S1C.png"),
    units = "mm",
    height = 70,
    width = 70)

################################################################################################################################################
########################################################       FIGURE 4 S1D     ###################################################################
################################################################################################################################################

library(DESeq2)
library(ggplot2)
library(patchwork)
load(here::here("data/figure_4/transcriptomics_DESeq_object/dds.RData"))

DE_results <- DE_genes |>
    dplyr::select(!1) |>
    # tibble::column_to_rownames("GENE") |>
    dplyr::mutate(
        padj = as.numeric(padj),
        log2FoldChange = as.numeric(log2FoldChange)
    )

# Function for plotting ---------------------------------------------------

visualizer_3genes <- function(data1,
                              gene1,
                              position_p_val_1,
                              plot_lim_1,
                              data2,
                              gene2,
                              position_p_val_2,
                              plot_lim_2,
                              data3,
                              gene3,
                              position_p_val_3,
                              plot_lim_3,
                              grouping,
                              title) {

    plot1 <- data1 |>
        dplyr::mutate(
            fibertype = factor(fibertype, levels = c("slow", "fast"))
        ) |>
        ggplot(aes(x=fibertype, y=count)) +
        geom_boxplot(aes(fill = fibertype), alpha = 0.85, outlier.shape = 21) +
        ggplot2::geom_point(aes(x = fibertype, y = count, fill = fibertype), shape = 21, stroke = 0.5, color = "black") +
        ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = position_p_val_1, label = paste("P.adj = ",
                                                                                              # round(
                                                                                              DE_results |>
                                                                                                  dplyr::filter(GENE == gene1) |>
                                                                                                  dplyr::pull(padj)
                                                                                              # 3
                                                                                              # )
        ),
        label.size = 2) +
        scale_fill_manual(values=c("#440154FF", "#5DC863FF")) +
        ggtitle(gene1) +
        labs(
            y = "Normalized counts",
            x = ""
        ) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(size = 7),
            text = element_text(face="bold", colour="black", size=7),
            plot.title = element_text(size = 7, hjust = 0.5),
        ) +
        ylim(plot_lim_1)

    plot2 <- data2 |>
        dplyr::mutate(
            fibertype = factor(fibertype, levels = c("slow", "fast"))
        ) |>
        ggplot(aes(x=fibertype, y=count)) +
        geom_boxplot(aes(fill = fibertype), alpha = 0.85, outlier.shape = 21) +
        ggplot2::geom_point(aes(x = fibertype, y = count, fill = fibertype), shape = 21, stroke = 0.5, color = "black") +
        ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = position_p_val_2, label = paste("P.adj = ",
                                                                                              # round(
                                                                                              DE_results |>
                                                                                                  dplyr::filter(GENE == gene2) |>
                                                                                                  dplyr::pull(padj)
                                                                                              #     3
                                                                                              # )
        ),
        label.size = 2) +
        scale_fill_manual(values=c("#440154FF", "#5DC863FF")) +
        ggtitle(gene2) +
        labs(
            y = "Normalized counts",
            x = ""
        ) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(size = 7),
            text = element_text(face="bold", colour="black", size=7),
            plot.title = element_text(size = 7, hjust = 0.5),
        ) +
        ylim(plot_lim_2)

    plot3 <- data3 |>
        dplyr::mutate(
            fibertype = factor(fibertype, levels = c("slow", "fast"))
        ) |>
        ggplot(aes(x=fibertype, y=count)) +
        geom_boxplot(aes(fill = fibertype), alpha = 0.85, outlier.shape = 21) +

        ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = position_p_val_3, label = paste("P.adj = ",
                                                                                              # round(
                                                                                              DE_results |>
                                                                                                  dplyr::filter(GENE == gene3) |>
                                                                                                  dplyr::pull(padj)
                                                                                              #     3
                                                                                              # )
        ),
        label.size = 2) +
        ggplot2::geom_point(aes(x = fibertype, y = count, fill = fibertype), shape = 21, stroke = 0.5, color = "black") +
        scale_fill_manual(values=c("#440154FF", "#5DC863FF")) +
        ggtitle(gene3) +
        labs(
            y = "Normalized counts",
            x = ""
        ) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(size = 7),
            text = element_text(face="bold", colour="black", size=7),
            plot.title = element_text(size = 7, hjust = 0.5),
        ) +
        ylim(plot_lim_3)

    combined_plot <- patchwork::wrap_plots(plot1,
                                           plot2,
                                           plot3,
                                           nrow = 1, ncol = 3) +
        plot_annotation(
            title = title,
            theme = theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))
        )

    return(combined_plot)

}

# Get data for individual genes

PPP3CB <- plotCounts(dds, gene = "PPP3CB", intgroup="fibertype", returnData=TRUE)

PPP1R3D <- plotCounts(dds, gene = "PPP1R3D", intgroup="fibertype", returnData=TRUE)

PPP3CA <- plotCounts(dds, gene = "PPP3CA", intgroup="fibertype", returnData=TRUE)

PPP1R3A <- plotCounts(dds, gene = "PPP1R3A", intgroup="fibertype", returnData=TRUE)

# Get plot

# phosphatase_plot <-

visualizer_3genes(
    data1 = PPP3CB,
    gene1 = "PPP3CB",
    position_p_val_1 = 380,
    plot_lim_1 = c(100, 395),
    data2 = PPP1R3D,
    gene2 = "PPP1R3D",
    plot_lim_2 = c(0, 35),
    position_p_val_2 = 30,
    data3 = PPP3CA,
    gene3 = "PPP3CA",
    position_p_val_3 = 250,
    plot_lim_3 = c(50, 270),
    grouping = "fibertype",
    title = "Protein Serine/Threonine phosphatase Transcriptomics"
)

ggsave(
    # phosphatase_plot,
    filename = here::here("doc/figures/figure_4_S1/figure_4_S1D.pdf"),
    width = 90,
    height = 45,
    units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1G     ###################################################################
################################################################################################################################################

# Get data for individual genes

SREBF1 <- plotCounts(dds, gene = "SREBF1", intgroup="fibertype", returnData=TRUE)

RXRG <- plotCounts(dds, gene = "RXRG", intgroup="fibertype", returnData=TRUE)

RXRA <- plotCounts(dds, gene = "RXRA", intgroup="fibertype", returnData=TRUE)

RORA <- plotCounts(dds, gene = "RORA", intgroup="fibertype", returnData=TRUE)

# Get plot

# TF_plot <-

visualizer_3genes(
    data1 = SREBF1,
    gene1 = "SREBF1",
    position_p_val_1 = 40,
    plot_lim_1 = c(15, 43),
    data2 = RXRG,
    gene2 = "RXRG",
    position_p_val_2 = 90,
    plot_lim_2 = c(5, 100),
    data3 =  RORA,
    gene3 = "RORA",
    position_p_val_3 = 1000,
    plot_lim_3 = c(200, 1100),
    grouping = "fibertype",
    title = "Ligand-activated TF activity Transcriptomics"
)

# Save plot
ggsave(
    filename = here::here("doc/figures/figure_4_S1/figure_4_S1G.pdf"),
    width = 90,
    height = 45,
    units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1H     ###################################################################
################################################################################################################################################

# Get data for individual genes

SCCPDH <- plotCounts(dds, gene = "SCCPDH", intgroup="fibertype", returnData=TRUE)

BDH1 <- plotCounts(dds, gene = "BDH1", intgroup="fibertype", returnData=TRUE)

DCXR <- plotCounts(dds, gene = "DCXR", intgroup="fibertype", returnData=TRUE)

TXN2 <- plotCounts(dds, gene = "TXN2", intgroup="fibertype", returnData=TRUE)

# Get plot

# oxidoreductase_plot <-
visualizer_3genes(
    data1 = BDH1,
    gene1 = "BDH1",
    position_p_val_1 = 40,
    plot_lim_1 = c(10, 45),
    data2 = DCXR,
    gene2 = "DCXR",
    position_p_val_2 = 300,
    plot_lim_2 = c(50, 320),
    data3 = TXN2,
    gene3 = "TXN2",
    position_p_val_3 = 165,
    plot_lim_3 = c(60, 170),
    grouping = "fibertype",
    title = "Oxidoreductase activity Transcriptomics"
)

# Save plot
ggsave(
    filename = here::here("doc/figures/figure_4_S1/figure_4_S1H.pdf"),
    width = 90,
    height = 45,
    units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1J     ###################################################################
################################################################################################################################################

# Get data for individual genes

HSPG2 <- plotCounts(dds, gene = "HSPG2", intgroup="fibertype", returnData=TRUE)

CTSD <- plotCounts(dds, gene = "CTSD", intgroup="fibertype", returnData=TRUE)

ADAMTSL4 <- plotCounts(dds, gene = "ADAMTSL4", intgroup="fibertype", returnData=TRUE)

LAMC1 <- plotCounts(dds, gene = "LAMC1", intgroup="fibertype", returnData=TRUE)

# Get plot

# ecm_plot <-
visualizer_3genes(
    data1 = CTSD,
    gene1 = "CTSD",
    position_p_val_1 = 190,
    plot_lim_1 = c(70, 200),
    data2 = ADAMTSL4,
    gene2 = "ADAMTSL4",
    position_p_val_2 = 55,
    plot_lim_2 = c(5, 62),
    data3 = LAMC1,
    gene3 = "LAMC1",
    position_p_val_3 = 60,
    plot_lim_3 = c(20, 70),
    grouping = "fibertype",
    title = "Extracellular matrix Transcriptomics"
)

# Save plot
ggsave(
    filename = here::here("doc/figures/figure_4_S1/figure_4_S1J.pdf"),
    width = 90,
    height = 45,
    units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1I     ###################################################################
################################################################################################################################################

# Get data for individual genes

CPTP <- plotCounts(dds, gene = "CPTP", intgroup="fibertype", returnData=TRUE)

SCP2 <- plotCounts(dds, gene = "SCP2", intgroup="fibertype", returnData=TRUE)

PFDN2 <- plotCounts(dds, gene = "PFDN2", intgroup="fibertype", returnData=TRUE)

CRYAB <- plotCounts(dds, gene = "CRYAB", intgroup="fibertype", returnData=TRUE)

# Get plot

# amide_plot <-
visualizer_3genes(
    data1 = CPTP,
    gene1 = "CPTP",
    position_p_val_1 = 75,
    plot_lim_1 = c(20, 80),
    data2 = PFDN2,
    gene2 = "PFDN2",
    position_p_val_2 = 80,
    plot_lim_2 = c(35, 90),
    data3 = CRYAB,
    gene3 = "CRYAB",
    position_p_val_3 = 1950,
    plot_lim_3 = c(500, 2000),
    grouping = "fibertype",
    title = "Amide binding Transcriptomics"
)

# Save plot
ggsave(
    filename = here::here("doc/figures/figure_4_S1/figure_4_S1I.pdf"),
    width = 90,
    height = 45,
    units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1K     ###################################################################
################################################################################################################################################

# Get data for individual genes

FNDC5 <- plotCounts(dds, gene = "FNDC5", intgroup="fibertype", returnData=TRUE)

SPX <- plotCounts(dds, gene = "SPX", intgroup="fibertype", returnData=TRUE)

NENF <- plotCounts(dds, gene = "NENF", intgroup="fibertype", returnData=TRUE)

FLRT2 <- plotCounts(dds, gene = "FLRT2", intgroup="fibertype", returnData=TRUE)

# Make final plot

# ligand_plot <-
visualizer_3genes(
    data1 = FNDC5,
    gene1 = "FNDC5",
    position_p_val_1 = 450,
    plot_lim_1 = c(50, 500),
    data2 = SPX,
    gene2 = "SPX",
    position_p_val_2 = 23,
    plot_lim_2 = c(5, 27),
    data3 = NENF,
    gene3 = "NENF",
    position_p_val_3 = 120,
    plot_lim_3 = c(40, 130),
    grouping = "fibertype",
    title = "Receptor ligand activity Transcriptomics"
)

# Save plot
ggsave(
    filename = here::here("doc/figures/figure_4_S1/figure_4_S1K.pdf"),
    width = 90,
    height = 45,
    units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1F     ###################################################################
################################################################################################################################################

data_pseudobulk <- vroom::vroom(here::here("data/data_proteomics_pseudobulk.csv")) |>
    tibble::column_to_rownames("...1")

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

DE_results <- DE_proteins |>
    tibble::column_to_rownames("Genes")

YTHDF3 <- data_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "YTHDF3") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping_pseudobulk |>
                          dplyr::select(c(sample_ID, fiber_type)))

YTHDF3$fiber_type <- factor(YTHDF3$fiber_type, levels = c("slow", "fast"))

YTHDF3_plot <- ggplot2::ggplot(YTHDF3,
                               ggplot2::aes(x = fiber_type,
                                            y = YTHDF3)) +
    ggplot2::geom_boxplot(aes(fill = fiber_type),
                          alpha = 0.85) +
    ggplot2::geom_point(aes(x = fiber_type, y = YTHDF3, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 12.3, label = paste("adj.P.val = ", round(
        DE_results |>
            tibble::rownames_to_column("Gene.name") |>
            dplyr::filter(Gene.name == "YTHDF3") |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    # ggplot2::scale_color_manual(values=c("#440154FF",
    #                                     "#5DC863FF")) +
    ggplot2::ggtitle("YTHDF3") +
    ggplot2::labs(
        y = "Norm. LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
    ) +
    ggplot2::ylim(NA, 12.4)


# TRIM21 --------------------------------------------------------------------

TRIM21 <- data_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "TRIM21") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping_pseudobulk |>
                          dplyr::select(c(sample_ID, fiber_type)))

TRIM21$fiber_type <- factor(TRIM21$fiber_type, levels = c("slow", "fast"))

TRIM21_plot <- ggplot2::ggplot(TRIM21,
                               ggplot2::aes(x = fiber_type,
                                            y = TRIM21)) +
    ggplot2::geom_boxplot(aes(fill = fiber_type),
                          alpha = 0.85) +
    ggplot2::geom_point(aes(x = fiber_type, y = TRIM21, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 12.5, label = paste("adj.P.val = ", round(
        DE_results |>
            tibble::rownames_to_column("Gene.name") |>
            dplyr::filter(Gene.name == "TRIM21") |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("TRIM21") +
    ggplot2::labs(
        y = "Norm. LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
    ) +
    ggplot2::ylim(NA, 12.7)

# LSM2 --------------------------------------------------------------------

LSM2 <- data_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "LSM2") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping_pseudobulk |>
                          dplyr::select(c(sample_ID, fiber_type)))

LSM2$fiber_type <- factor(LSM2$fiber_type, levels = c("slow", "fast"))

LSM2_plot <- ggplot2::ggplot(LSM2,
                             ggplot2::aes(x = fiber_type,
                                          y = LSM2)) +
    ggplot2::geom_boxplot(aes(fill = fiber_type),
                          alpha = 0.85) +
    ggplot2::geom_point(aes(x = fiber_type, y = LSM2, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 11.85, label = paste("adj.P.val = ", round(
        DE_results |>
            tibble::rownames_to_column("Gene.name") |>
            dplyr::filter(Gene.name == "LSM2") |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +

    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("LSM2") +
    ggplot2::labs(
        y = "Norm. LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
    ) +
    ggplot2::ylim(NA, 11.9)



# Joint plot --------------------------------------------------------------

p_body_plot <- patchwork::wrap_plots(YTHDF3_plot,
                                     TRIM21_plot,
                                     LSM2_plot,
                                     nrow = 1, ncol = 3) +
    patchwork::plot_annotation(
        title = "P-body Proteomics",
        theme = theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))
    )

ggsave(p_body_plot,
       filename = here::here("doc/figures/figure_4_S1/figure_4_S1F.pdf"),
       width = 90,
       height = 45,
       units="mm")

################################################################################################################################################
########################################################       FIGURE 4 S1E     ###################################################################
################################################################################################################################################

# PPP3CB -------------------------------------------------------------------

PPP3CB <- data_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "PPP3CB") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping_pseudobulk |>
                          dplyr::select(c(sample_ID, fiber_type)))

PPP3CB$fiber_type <- factor(PPP3CB$fiber_type, levels = c("slow", "fast"))

PPP3CB_plot <- ggplot2::ggplot(PPP3CB,
                               ggplot2::aes(x = fiber_type,
                                            y = PPP3CB)) +
    ggplot2::geom_boxplot(aes(fill = fiber_type),
                          alpha = 0.85) +
    ggplot2::geom_point(aes(x = fiber_type, y = PPP3CB, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 14.2, label = paste("adj.P.val = ", round(
        DE_results |>
            tibble::rownames_to_column("Gene.name") |>
            dplyr::filter(Gene.name == "PPP3CB") |>
            dplyr::pull(adj.P.Val),
        5
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("PPP3CB") +
    ggplot2::labs(
        y = "Norm. LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
    ) +
    ggplot2::ylim(NA, 14.35)

# PPP1R3D --------------------------------------------------------------------

PPP1R3D <- data_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "PPP1R3D") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping_pseudobulk |>
                          dplyr::select(c(sample_ID, fiber_type)))

PPP1R3D$fiber_type <- factor(PPP1R3D$fiber_type, levels = c("slow", "fast"))

PPP1R3D_plot <- ggplot2::ggplot(PPP1R3D,
                                ggplot2::aes(x = fiber_type,
                                             y = PPP1R3D)) +
    ggplot2::geom_boxplot(aes(fill = fiber_type),
                          alpha = 0.85) +
    ggplot2::geom_point(aes(x = fiber_type, y = PPP1R3D, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 11.35, label = paste("adj.P.val = ", round(
        DE_results |>
            tibble::rownames_to_column("Gene.name") |>
            dplyr::filter(Gene.name == "PPP1R3D") |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("PPP1R3D") +
    ggplot2::labs(
        y = "Norm. LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5)
    ) +
    ggplot2::ylim(NA, 11.45)

# PPP1R3A --------------------------------------------------------------------

PPP1R3A <- data_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "PPP1R3A") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping_pseudobulk |>
                          dplyr::select(c(sample_ID, fiber_type)))

PPP1R3A$fiber_type <- factor(PPP1R3A$fiber_type, levels = c("slow", "fast"))

PPP1R3A_plot <- ggplot2::ggplot(PPP1R3A,
                                ggplot2::aes(x = fiber_type,
                                             y = PPP1R3A)) +
    ggplot2::geom_boxplot(aes(fill = fiber_type),
                          alpha = 0.85) +
    ggplot2::geom_point(aes(x = fiber_type, y = PPP1R3A, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 14.25, label = paste("adj.P.val = ", round(
        DE_results |>
            tibble::rownames_to_column("Gene.name") |>
            dplyr::filter(Gene.name == "PPP1R3A") |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("PPP1R3A") +
    ggplot2::labs(
        y = "Norm. LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
    ) +
    ggplot2::ylim(NA, 14.35)

# Joint plot --------------------------------------------------------------

phosphatase_plot <- patchwork::wrap_plots(PPP3CB_plot,
                                          PPP1R3D_plot,
                                          PPP1R3A_plot,
                                          nrow = 1, ncol = 3) +
    patchwork::plot_annotation(
        title = "Protein Serine/Threonine phosphatase Proteomics",
        theme = theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))
    )

ggsave(phosphatase_plot,
       filename = here::here("doc/figures/figure_4_S1/figure_4_S1E.pdf"),
       width = 90,
       height = 45,
       units="mm")
