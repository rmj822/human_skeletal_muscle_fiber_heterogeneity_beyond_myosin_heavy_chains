
# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(cowplot)
library(ggpubr)
library(VennDiagram)
library(clusterProfiler)

################################################################################################################################################
##############################################     LOAD AND PREP TRANSCRIPTOME DATA    #########################################################
################################################################################################################################################

# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

# Set identify of clusters
Idents(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"

# Set to SCT assay for visualization, and log-transform data ------------------------------
DefaultAssay(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "SCT"

################################################################################################################################################
##############################################     FEATURE PLOTS TRANSCRIPTOMICS WITHOUT LEGENDS   #############################################
################################################################################################################################################

# MYH7
feature_plot_transcriptomics_MYH7 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("MYH7"),
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

# MYH2
feature_plot_transcriptomics_MYH2 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("MYH2"),
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

# MYH1
feature_plot_transcriptomics_MYH1 <- Seurat::FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                         features = c("MYH1"),
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


################################################################################################################################################
###############################################     FEATURE PLOTS PROTEOMICS WITHOUT LEGENDS   #################################################
################################################################################################################################################

# MYH7
feature_plot_proteomics_MYH7 <- Seurat::FeaturePlot(object = seurat_proteome,
                                                         features = c("MYH7"),
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

# MYH2
feature_plot_proteomics_MYH2 <- Seurat::FeaturePlot(object = seurat_proteome,
                                                    features = c("MYH2"),
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

# MYH1
feature_plot_proteomics_MYH1 <- Seurat::FeaturePlot(object = seurat_proteome,
                                                    features = c("MYH1"),
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


################################################################################################################################################
#####################################################     COMBINE PLOTS TOGETHER   #############################################################
################################################################################################################################################

# Transcriptomics
feature_plot_transcriptomics <- ggpubr::ggarrange(
    feature_plot_transcriptomics_MYH7,
    feature_plot_transcriptomics_MYH2,
    feature_plot_transcriptomics_MYH1,
    ncol = 3,
    nrow = 1)

feature_plot_transcriptomics <- annotate_figure(feature_plot_transcriptomics,
                                                top = text_grob("Transcriptomics",
                                                                color = "black",
                                                                face = "bold",
                                                                size = 8)
                                                )

ggsave(feature_plot_transcriptomics,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_all_genes/6_PCs/Featureplot_MYH7_MYH2_MYH1_transcriptomics.png",
       width = 100,
       height = 35,
       units="mm")

# Proteomics
feature_plot_proteomics <- ggpubr::ggarrange(
    feature_plot_proteomics_MYH7,
    feature_plot_proteomics_MYH2,
    feature_plot_proteomics_MYH1,
    ncol = 3,
    nrow = 1)

feature_plot_proteomics <- annotate_figure(feature_plot_proteomics,
                                                top = text_grob("Proteomics",
                                                                color = "black",
                                                                face = "bold",
                                                                size = 8)
)

# Combine Transcriptomics and proteomics
feature_plot_combined <- ggpubr::ggarrange(
    feature_plot_transcriptomics,
    feature_plot_proteomics,
    ncol = 2,
    nrow = 1)


################################################################################################################################################
#####################################################     GET LEGENDS TRANSCRIPTOMICS   ########################################################
################################################################################################################################################

# MYH7
feature_plot_transcriptomics_MYH7_legend <- feature_plot_transcriptomics_MYH7 +
    ggplot2::theme(
        legend.position = "top",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_transcriptomics_MYH7_legend <- ggpubr::get_legend(feature_plot_transcriptomics_MYH7_legend)
feature_plot_transcriptomics_MYH7_legend <- ggpubr::as_ggplot(feature_plot_transcriptomics_MYH7_legend)

ggsave(feature_plot_transcriptomics_MYH7_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_all_genes/6_PCs/Featureplot_MYH7_MYH2_MYH1_transcriptomics_legend_MYH7.png",
       width = 90,
       height = 40,
       units="mm")

# MYH2
feature_plot_transcriptomics_MYH2_legend <- feature_plot_transcriptomics_MYH2 +
    ggplot2::theme(
        legend.position = "top",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_transcriptomics_MYH2_legend <- ggpubr::get_legend(feature_plot_transcriptomics_MYH2_legend)
feature_plot_transcriptomics_MYH2_legend <- ggpubr::as_ggplot(feature_plot_transcriptomics_MYH2_legend)

ggsave(feature_plot_transcriptomics_MYH2_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_all_genes/6_PCs/Featureplot_MYH7_MYH2_MYH1_transcriptomics_legend_MYH2.png",
       width = 90,
       height = 40,
       units="mm")

# MYH1
feature_plot_transcriptomics_MYH1_legend <- feature_plot_transcriptomics_MYH1 +
    ggplot2::theme(
        legend.position = "top",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

feature_plot_transcriptomics_MYH1_legend <- ggpubr::get_legend(feature_plot_transcriptomics_MYH1_legend)
feature_plot_transcriptomics_MYH1_legend <- ggpubr::as_ggplot(feature_plot_transcriptomics_MYH1_legend)

ggsave(feature_plot_transcriptomics_MYH1_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_2/transcriptomics_all_genes/6_PCs/Featureplot_MYH7_MYH2_MYH1_transcriptomics_legend_MYH1.png",
       width = 90,
       height = 40,
       units="mm")

################################################################################################################################################
##################################     TRANSCRIPTOMICS: BLENDED UMAP MYH2 - MYH1 (For transcriptomics)  ########################################
################################################################################################################################################

# Get counts data
data_transcriptomics <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                       assay = "SCT",
                       slot = "counts")

data_transcriptomics <- as.data.frame(data_transcriptomics)

# Extract MYH counts
MYH1 <- c("MYH1")
counts_MYH1 <- data_transcriptomics[MYH1, ]
counts_MYH1 <- as.data.frame(as.matrix(counts_MYH1))
counts_MYH1 <- t(counts_MYH1)

MYH2 <- c("MYH2")
counts_MYH2 <- data_transcriptomics[MYH2, ]
counts_MYH2 <- as.data.frame(as.matrix(counts_MYH2))
counts_MYH2 <- t(counts_MYH2)

MYH7 <- c("MYH7")
counts_MYH7 <- data_transcriptomics[MYH7, ]
counts_MYH7 <- as.data.frame(as.matrix(counts_MYH7))
counts_MYH7 <- as.data.frame(t(counts_MYH7))
counts_MYH7$Row.names <- rownames(counts_MYH7)

# Merge MYH dataframes
counts_MYH <- merge(counts_MYH1, counts_MYH2, by="row.names")
counts_MYH <- merge(counts_MYH, counts_MYH7, by="Row.names")

# Delete Row.names column
rownames(counts_MYH) <- counts_MYH$Row.names
counts_MYH <- counts_MYH %>% dplyr::select(MYH1, MYH2, MYH7)

# Calculate percentages for each MYH
counts_MYH$sum_MYH <- counts_MYH$MYH7 + counts_MYH$MYH2 + counts_MYH$MYH1

counts_MYH$MYH7_fraction <- counts_MYH$MYH7 / counts_MYH$sum_MYH * 100
counts_MYH$MYH2_fraction <- counts_MYH$MYH2 / counts_MYH$sum_MYH * 100
counts_MYH$MYH1_fraction <- counts_MYH$MYH1 / counts_MYH$sum_MYH * 100

# Calculate difference between MYH2 and MYH1
counts_MYH$MYH2_MYH1_diff <- counts_MYH$MYH2_fraction - counts_MYH$MYH1_fraction
counts_MYH$MYHs_fast <- counts_MYH$MYH2 + counts_MYH$MYH1

# Extract UMAP coordinates and add MYH diff colouring
df_UMAP_blended <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest[["umap"]]@cell.embeddings %>%
    as.data.frame() %>%
    rownames_to_column("row.names")

df_UMAP_blended <- counts_MYH %>%
    dplyr::select(MYH2_MYH1_diff, MYHs_fast) %>%
    rownames_to_column("row.names") %>%
    left_join(df_UMAP_blended, by = "row.names")

# UMAP coloured by MYH
UMAP_blended <- df_UMAP_blended %>%
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color =MYH2_MYH1_diff)
    ) +
    ggplot2::geom_point(aes(alpha = MYHs_fast),
                        size = 1) +
    scale_alpha(range = c(0.2, 1),
                guide = "none") +
    ggplot2::scale_color_gradient2(
        "MYH2 % \n     - \nMYH1 %",
        low = "#D2631C",
        mid = "lightgrey",
        high = "#23CE6B",
        midpoint = 0,
        limits = c(-100, 100)
    ) +
    ggplot2::ggtitle("UMAP transcriptomics - by fast MYH blended expression") +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, "mm"),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 7, face = "bold")
        # plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm")
    )

ggplot2::ggsave(
    UMAP_blended,
    filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/fast_MYHs_blended_transcriptomics.png",
    units = "mm",
    height = 60,
    width = 90
)

################################################################################################################################################
#######################     TRANSCRIPTOMICS: BLENDED UMAP MYH2 - MYH1 - MYH7 (For transcriptomics) - RGB colouring  ############################
################################################################################################################################################

perc_MYHs <- counts_MYH |>
    rownames_to_column("fiber_ID") %>%
    dplyr::select(c(fiber_ID, MYH7_fraction, MYH2_fraction, MYH1_fraction)) |>
    dplyr::mutate(
        dplyr::across(
            .cols = !fiber_ID,
            ~ .x/100
        )
    )

perc_list <- perc_MYHs |>
    dplyr::group_split(
        fiber_ID
    )

names(perc_list) <- perc_MYHs |>
    dplyr::arrange(fiber_ID) |>
    dplyr::pull(fiber_ID)

rgb_maker <- function(.data, .fiber_ID){
    tmp2 <- .data[[.fiber_ID]]

    tmp2 <- tmp2 |>
        dplyr::mutate(
            color = rgb(
                red = tmp2$MYH1_fraction,
                green = tmp2$MYH2_fraction,
                blue = tmp2$MYH7_fraction,
                alpha = 0.7
            ))

    return(tmp2)
}

color_list <- purrr::map(
    names(perc_list),
    ~ rgb_maker(
        .data = perc_list,
        .fiber_ID = .x
    )
) |>
    purrr::list_rbind() |>
    tibble::column_to_rownames("fiber_ID")

# Extract UMAP coordinates and add MYH diff colouring
df_UMAP_blended_rgb <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest[["umap"]]@cell.embeddings %>%
    as.data.frame() %>%
    rownames_to_column("row.names")

df_UMAP_blended_rgb <- color_list %>%
    dplyr::select(color) %>%
    rownames_to_column("row.names") %>%
    left_join(df_UMAP_blended_rgb, by = "row.names")

df_UMAP_blended_rgb$color


# UMAP coloured by MYH RGB
df_UMAP_blended_rgb |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = color)
    ) +
    ggplot2::geom_point(
        # ggplot2::aes(alpha = coloring$MYHs_fast),
        size = 1,
        alpha = 0.65) +
    ggplot2::scale_color_identity() +
    # ggplot2::scale_alpha(range = c(0.2, 1),
    #                      guide = "none") +
    # ggplot2::scale_color_gradient2(
    #     "MYH2 % \n     - \nMYH1 %",
    #     low = "#D2631C",
    #     mid = "lightgrey",
    #     high = "#23CE6B",
    #     midpoint = 0,
    #     limits = c(-100, 100),
    #     breaks = c(-100, -50, 0, 50, 100)
    # ) +
ggplot2::ggtitle("UMAP transcriptomics - by MYH blended expression") +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6, face = "bold"),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, "mm"),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 7, face = "bold")
        # plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm")
    )

ggplot2::ggsave(
    here::here("doc/figures/figure_1/all_MYHs_blended_transcriptomics_scaled.pdf"),
    units = "mm",
    height = 60,
    width = 80
)


