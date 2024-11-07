
################################################################################################################################################
########################################################       FIGURE 1 S4A     ###################################################################
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


# ACTN3 -------------------------------------------------------------------

feature_plot_ACTN3 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                          features = c("ACTN3"),
                                          reduction = "umap",
                                          pt.size = 0.3,
                                          label.size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c("Norm.\ncounts", option = "plasma") +
    ggplot2::ggtitle("UMAP transcriptomics\nby ACTN3 expression") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.1, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

# MYLK2 -------------------------------------------------------------------

feature_plot_MYLK2 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                          features = c("MYLK2"),
                                          reduction = "umap",
                                          pt.size = 0.3,
                                          label.size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c("Norm.\ncounts", option = "plasma") +
    ggplot2::ggtitle("UMAP transcriptomics\nby MYLK2 expression") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.1, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        plot.margin = grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

patchwork::wrap_plots(feature_plot_ACTN3, feature_plot_MYLK2)

ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4A.png"),
       width = 99,
       height = 45,
       units="mm")

