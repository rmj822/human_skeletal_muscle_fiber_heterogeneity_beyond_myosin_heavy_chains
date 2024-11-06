################################################################################################################################################
########################################################      Panel A   ############################################################################
################################################################################################################################################

library(vroom)
library(ggplot2)
library(here)
library(ggrepel)

data_file <- here("data/Ribo_muscle_enrich.txt")
data <- vroom(data_file)
labels <- c("RPL38", "RPS13", "RPL18A")

plot <- ggplot(data, aes(x = Enrich, y = Rank, size = Enrich, color = Enrich)) +
    geom_point() +
    theme_minimal() +
    scale_color_gradient(low = "red", high = "blue") +
    scale_y_reverse() +
    scale_size(range=c(0.1,1.5)) +
    ylab("Rank") +
    xlab("log2FC\n (Muscle - whole body)") +
    ggtitle("Muscle-specific ribosomal gene signature") +
    theme(legend.position = "none",
          axis.text = element_text(size = 8, face = "bold"),
          axis.title = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5))


ggsave("doc/figures/figure_3_S1/figure_3_S1A.png", plot = plot, height = 60, width = 90, units = "mm")
################################################################################################################################################
########################################################      Panel B   ############################################################################
################################################################################################################################################


load(here::here("data/figure_2/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata"))

transcriptomics_scaledata <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                          assay = "SCT",
                                          slot = "scale.data")

# Extract SCT normalized but not log-transformed counts
transcriptomics_counts <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                       assay = "SCT",
                                       slot = "counts")

pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2)

data_pca$fiberID <- rownames(data_pca)

ribosome_filtered <- transcriptomics_counts |>
    as.data.frame() %>%
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_ribosome) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(transcriptomics_counts, na.rm = T)

rel_abundance_ribosomes <- (colSums(ribosome_filtered,
                                    na.rm = T) / total_intensity) * 100

data_pca_ribosomes <- as.data.frame(rel_abundance_ribosomes) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca)

data_ribosomal_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2) |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ .x * 160)) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% vector_ribosome)

data_ribosomal_proteins <- data_ribosomal_proteins |>
    dplyr::mutate("x.centroid" = 0) |>
    dplyr::mutate("y.centroid" = 0) |>
    dplyr::mutate(label = dplyr::case_when(
        Genes %in% c(
            "RPL18",
            "RPS13",
            "RPS18",
            "RPL38",
            "RPL35",
            "RPL31"
        ) ~ Genes,
        TRUE ~ ""))


custom_ribo_plot_with_legend <- data_pca_ribosomes |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = rel_abundance_ribosomes
    )) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_viridis_c("Rel. \nAbundance", option = "plasma") +
    ggplot2::geom_segment(
        data = data_ribosomal_proteins,
        ggplot2::aes(x = x.centroid,
                     y = y.centroid,
                     xend = Dim.1,
                     yend = Dim.2),
        color = "black",
        arrow = grid::arrow(type = "closed",
                            length = unit(1, "mm"))
    ) +
    ggrepel::geom_label_repel(data = data_ribosomal_proteins,
                              ggplot2::aes(
                                  x = Dim.1,
                                  y = Dim.2,
                                  label = label),
                              size = 2,
                              color = "black",
                              label.padding = 0.1,
                              min.segment.length = 0.1,
                              segment.size = 0.2,
                              force = 20,
                              max.overlaps = Inf) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA colored by Cytosolic Ribosome GO term \nTranscriptomics") +
    ggplot2::xlab("PC1 (11.1%)") +
    ggplot2::ylab("PC2 (3.5%)") +
    ggplot2::theme(
        text = element_text(face="bold", colour="black", size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 3),
        legend.position = "right",
        legend.key.width = ggplot2::unit(2, "mm")
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    ) +
    scale_x_continuous(trans = "reverse",
                       breaks=c(-50, 0, 50, 100),
                       labels=c("50", "0", "-50", "-100")
    )

ggsave(custom_ribo_plot_with_legend,
       filename = here::here("doc/figures/figure_3_S1/figure_3_S1B.png"),
       width = 90,
       height = 60,
       units="mm")

################################################################################################################################################
########################################################      Panel C   ############################################################################
################################################################################################################################################

source(here::here("R/figure_3/figure_3.R"))
median_ribosomes <- data_ribosomes |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tidyr::pivot_longer(
        cols = !sample_id,
        names_to = "genes",
        values_to = "LFQ"
    ) |>
    dplyr::group_by(genes) |>
    dplyr::summarise(
        median_lfq = median(LFQ, na.rm = TRUE)
    )

samples <- colnames(data_ribosomes)

data_dif_median <- data_ribosomes |>
    tibble::rownames_to_column("genes") |>
    tidyr::pivot_longer(
        cols = !genes,
        values_to = "LFQ",
        names_to = "sample_id"
    ) |>
    dplyr::inner_join(median_ribosomes) |>
    dplyr::mutate(dif_from_median = (LFQ - median_lfq)) |>
    dplyr::select(genes, sample_id, dif_from_median)

data_dif_median |>
    dplyr::mutate(
        coloring = dplyr::case_when(
            genes %in% c("RPS13", "RPL18", "RPS18") ~ "clust_2",
            genes %in% c("RPL38", "RPL35", "RPL31") ~ "clust_1",
            TRUE ~ ""
        )
    ) |>
    dplyr::filter(genes %in% c("RPS13", "RPL18", "RPS18", "RPL38", "RPL35", "RPL31")) |>
    dplyr::mutate(
        genes = factor(genes, levels = c("RPS13", "RPL18", "RPS18", "RPL38", "RPL35", "RPL31"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = dif_from_median,
            y = genes,
            fill = coloring
        )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::geom_violin(alpha = 0.5) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(fill = coloring),
        shape = 21,
        size = 0.75,
        stroke = 0.15,
        color = "black") +
    ggplot2::scale_fill_manual(
        values = c("#008080",
                   "#FB8022FF")
    ) +
    ggplot2::xlab("Difference from the median (Log2 LFQ intensity)") +
    ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        legend.position = "none",
        text = ggplot2::element_text(size = 6, face = "bold")
    )

ggsave("doc/figures/figure_3_S1/figure_3_S1C.png", height = 60, width = 90, units = "mm")

################################################################################################################################################
########################################################      Panel D   ############################################################################
################################################################################################################################################

# Load data ---------------------------------------------------------------

ribosomal_clusters <- ribosomal_clusters |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::mutate("ribosomal_subunit" = dplyr::case_when(
        grepl("L", Gene.name) ~ "Large subunit",
        grepl("S", Gene.name) ~ "Small subunit",
        TRUE ~ "unknown"
    )) |>
    dplyr::filter(!ribosomal_subunit == "unknown")

ribosomal_clusters |>
    ggplot2::ggplot(
        ggplot2::aes(
            # x = ribosomal_clusters,
            x = ribosomal_clusters,
            fill = ribosomal_clusters
        )
    ) +
    ggplot2::geom_bar() +
    ggplot2::scale_fill_manual("Ribosomal clusters",
                               values = c(ggplot2::alpha("#008080", 0.75),
                                          ggplot2::alpha("#FB8022FF", 0.75),
                                          ggplot2::alpha("grey", 0.75))) +
    ggplot2::theme_light() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x.bottom = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(vjust = -0.35),
        axis.title.y = ggplot2::element_text(vjust = 0.35),
        text = ggplot2::element_text(size = 6),
        legend.key.size = ggplot2::unit(3, units = "mm"),
        legend.position = "bottom"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 18),
                                expand = c(0, 0)) +
    ggplot2::facet_grid(~ ribosomal_subunit) +
    ggplot2::theme(
        strip.text.x = ggplot2::element_text(size = 8, face = "bold")
    )


# Percentage plot ---------------------------------------------------------

proportion_ribosomal_subunits <- data.frame(
    ribosomal_clusters |>
        dplyr::filter(ribosomal_clusters == "ribosomal_cluster_1") |>
        dplyr::group_by(ribosomal_subunit) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        as.data.frame() |>
        dplyr::mutate(cnt = as.numeric(cnt)) |>
        dplyr::mutate("ribosomal_cluster_1" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(ribosomal_subunit)) |>
        tibble::remove_rownames() |>
        tibble::column_to_rownames("ribosomal_subunit") |>
        dplyr::select(ribosomal_cluster_1),
    ribosomal_clusters |>
        dplyr::filter(ribosomal_clusters == "ribosomal_cluster_2") |>
        dplyr::group_by(ribosomal_subunit) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        as.data.frame() |>
        dplyr::mutate(cnt = as.numeric(cnt)) |>
        dplyr::mutate("ribosomal_cluster_2" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(ribosomal_subunit)) |>
        tibble::remove_rownames() |>
        tibble::column_to_rownames("ribosomal_subunit") |>
        dplyr::select(ribosomal_cluster_2),
    ribosomal_clusters |>
        dplyr::filter(ribosomal_clusters == "ribosomal_cluster_3") |>
        dplyr::group_by(ribosomal_subunit) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        as.data.frame() |>
        dplyr::mutate(cnt = as.numeric(cnt)) |>
        dplyr::mutate("ribosomal_cluster_3" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(ribosomal_subunit)) |>
        tibble::remove_rownames() |>
        tibble::column_to_rownames("ribosomal_subunit") |>
        dplyr::select(ribosomal_cluster_3)
) |>
    tibble::rownames_to_column("ribosomal_subunit") |>
    tidyr::pivot_longer(
        cols = c(2,3,4),
        names_to = "ribosomal_clusters",
        values_to = "percentage"
    )

proportion_ribosomal_subunits |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = ribosomal_clusters,
            y = percentage,
            fill = ribosomal_clusters
        )
    ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual("Ribosomal clusters",
                               values = c(ggplot2::alpha("#008080", 0.75),
                                          ggplot2::alpha("#FB8022FF", 0.75),
                                          ggplot2::alpha("grey", 0.75))) +
    ggplot2::theme_light() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x.bottom = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(vjust = -0.35),
        axis.title.y = ggplot2::element_text(vjust = 0.35, face = "bold"),
        text = ggplot2::element_text(size = 6),
        legend.key.size = ggplot2::unit(3, units = "mm"),
        legend.position = "bottom"
    ) +
    ggplot2::labs(y = "Percentage") +
    ggplot2::scale_y_continuous(limits = c(0, 80),
                                expand = c(0, 0)) +
    ggplot2::facet_grid(~ ribosomal_subunit) +
    ggplot2::theme(
        strip.text.x = ggplot2::element_text(size = 8, face = "bold")
    )

ggplot2::ggsave(here::here("doc/figures/figure_3_S1/figure_3_S1D.png"),
                height = 60,
                width = 90,
                units = "mm")
