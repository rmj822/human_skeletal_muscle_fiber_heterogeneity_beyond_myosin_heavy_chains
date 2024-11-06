
# Load proteomics data and create Seurat object ---------------------------

data_proteomics <- vroom::vroom(here::here("data/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "...1") |>
    tibble::column_to_rownames("Gene.name")

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

# ggplot2::ggsave(here::here("doc/figures/figure_1/scree_plot_proteomics.png"),
#                 units = "mm",
#                 height = 60,
#                 width = 90)
