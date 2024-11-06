pc1_pos <- vroom::vroom(here::here("data/GSEA_PCA_proteomics/results_PC1_positive.csv")) |>
    dplyr::slice_head(n = 3)
pc1_neg <- vroom::vroom(here::here("data/GSEA_PCA_proteomics/results_PC1_negative.csv")) |>
    dplyr::slice_head(n = 3)

pc2_pos <- vroom::vroom(here::here("data/GSEA_PCA_proteomics/results_PC2_positive.csv")) |>
    dplyr::slice_head(n = 3)
pc2_neg <- vroom::vroom(here::here("data/GSEA_PCA_proteomics/results_PC2_negative.csv")) |>
    dplyr::slice_head(n = 4)

# Calculate foldEnrich for proteomics
pc1_pos <- dplyr::mutate(pc1_pos, foldEnrich =
                         (as.numeric(sub("/\\d+", "", pc1_pos$GeneRatio)) / as.numeric(sub(".*/", "", pc1_pos$GeneRatio))) /
                         (as.numeric(sub("/\\d+", "", pc1_pos$BgRatio)) / as.numeric(sub(".*/", "", pc1_pos$BgRatio))))

pc1_neg <- dplyr::mutate(pc1_neg, foldEnrich =
                             (as.numeric(sub("/\\d+", "", pc1_neg$GeneRatio)) / as.numeric(sub(".*/", "", pc1_neg$GeneRatio))) /
                             (as.numeric(sub("/\\d+", "", pc1_neg$BgRatio)) / as.numeric(sub(".*/", "", pc1_neg$BgRatio))))

pc2_pos <- dplyr::mutate(pc2_pos, foldEnrich =
                             (as.numeric(sub("/\\d+", "", pc2_pos$GeneRatio)) / as.numeric(sub(".*/", "", pc2_pos$GeneRatio))) /
                             (as.numeric(sub("/\\d+", "", pc2_pos$BgRatio)) / as.numeric(sub(".*/", "", pc2_pos$BgRatio))))

pc2_neg <- dplyr::mutate(pc2_neg, foldEnrich =
                             (as.numeric(sub("/\\d+", "", pc2_neg$GeneRatio)) / as.numeric(sub(".*/", "", pc2_neg$GeneRatio))) /
                             (as.numeric(sub("/\\d+", "", pc2_neg$BgRatio)) / as.numeric(sub(".*/", "", pc2_neg$BgRatio))))

# Annotate every dataframe before merging:

pc1_pos <- pc1_pos |>
    dplyr::mutate(
        direction = "PC1_pos"
    )

pc1_neg <- pc1_neg |>
    dplyr::mutate(
        direction = "PC1_neg"
    )

pc2_pos <- pc2_pos |>
    dplyr::mutate(
        direction = "PC2_pos"
    )

pc2_neg <- pc2_neg |>
    dplyr::mutate(
        direction = "PC2_neg"
    )

PC_plot <- dplyr::bind_rows(pc2_pos,
                            pc2_neg,
                            pc1_pos,
                            pc1_neg)

PC_plot$graph <- rep(1, nrow(PC_plot))

ggplot2::ggplot(
    PC_plot, ggplot2::aes(graph,
                           forcats::fct_reorder(Description, direction))) +
    ggplot2::geom_point(ggplot2::aes(color=foldEnrich, size = p.adjust)) +
    ggplot2::scale_color_gradient("Fold \nenrichment", low = "#6baed6",high = "#08519c") +
    ggplot2::scale_size(trans = 'reverse', range=c(2, 3)) +
    ggplot2::facet_grid(~direction) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(limits = c(0.9,1.1), expand = c(0, 0)) +
    ggplot2::theme(
        text = ggplot2::element_text(face = "bold",size = 6, colour = "black"),
        axis.title.y= ggplot2::element_blank(),
        axis.title.x= ggplot2::element_blank(),
        axis.text.x= ggplot2::element_blank(),
        axis.line.x =ggplot2::element_blank(),
        axis.ticks.x =ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face="bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, "mm"),
        strip.background = ggplot2::element_rect(fill="#bdd7e7")
    ) +
    ggplot2::coord_fixed(ratio = 0.2)

ggplot2::ggsave(here::here("doc/figures/figure_3/dot_plot_PC_proteomics.png"),
       width = 100,
       height = 120,
       units="mm")
