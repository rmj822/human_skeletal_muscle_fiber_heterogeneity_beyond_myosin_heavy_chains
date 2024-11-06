source(here::here("R/figure_2/MYH_curves.r"))

metadata <- vroom::vroom(
    "C:/Users/jns822/Desktop/Scripts/Heterofiber/data-raw/metadata.txt"
)

metadata <- metadata |>
    dplyr::filter(fiberID %in% data_fiber_type$fiber_ID,
                  !duplicated(fiberID))

metadata_fiber_type <- metadata |>
    dplyr::inner_join(data_fiber_type |>
                          dplyr::rename("fiberID" = fiber_ID) |>
                          dplyr::select(fiber_type, fiberID)) |>
    dplyr::select(subject,
                  fiberID,
                  fiber_type)

    data.frame(
        metadata_fiber_type |>
            dplyr::filter(subject == "FOR2") |>
            dplyr::group_by(fiber_type) |>
            dplyr::summarise(cnt = dplyr::n()) |>
            dplyr::mutate("FOR2" = round(cnt / sum(cnt), 3) * 100) |>
            tibble::column_to_rownames("fiber_type") |>
            dplyr::select(FOR2),
        metadata_fiber_type |>
            dplyr::filter(subject == "FOR4") |>
            dplyr::group_by(fiber_type) |>
            dplyr::summarise(cnt = dplyr::n()) |>
            dplyr::mutate("FOR4" = round(cnt / sum(cnt), 3) * 100) |>
            tibble::column_to_rownames("fiber_type") |>
            dplyr::select(FOR4),
        metadata_fiber_type |>
            dplyr::filter(subject == "FOR9") |>
            dplyr::group_by(fiber_type) |>
            dplyr::summarise(cnt = dplyr::n()) |>
            dplyr::mutate("FOR9" = round(cnt / sum(cnt), 3) * 100) |>
            tibble::column_to_rownames("fiber_type") |>
            dplyr::select(FOR9),
        metadata_fiber_type |>
            dplyr::filter(subject == "FOR10") |>
            dplyr::group_by(fiber_type) |>
            dplyr::summarise(cnt = dplyr::n()) |>
            dplyr::mutate("FOR10" = round(cnt / sum(cnt), 3) * 100) |>
            tibble::column_to_rownames("fiber_type") |>
            dplyr::select(FOR10),
        metadata_fiber_type |>
            dplyr::filter(subject == "FOR11") |>
            dplyr::group_by(fiber_type) |>
            dplyr::summarise(cnt = dplyr::n()) |>
            dplyr::mutate("FOR11" = round(cnt / sum(cnt), 3) * 100) |>
            tibble::column_to_rownames("fiber_type") |>
            dplyr::select(FOR11)
    ) |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("subject") |>
        tidyr::pivot_longer(
            cols = c("Hybrid 1/2A", "Hybrid 2A/2X", "Type 1", "Type 2A"),
            names_to = "fiber_type"
        ) |>
        ggplot2::ggplot(ggplot2::aes(
            x = subject,
            y = value,
            fill = fiber_type
        )) +
        ggplot2::geom_col(alpha = 0.85,
                          width = 0.5,
                          color = "black") +
        # ggplot2::facet_grid(~subject) +
        ggplot2::scale_fill_manual("Fiber types",
                                   values = c("#3B528BFF", "#FDE725FF", "#440154FF", "#5DC863FF")
        ) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle("Proteomics fiber type \nby participant") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.title = ggplot2::element_text(size = 8,
                                                            face = "bold"),
                       legend.key.size = ggplot2::unit(3, "mm"),
                       legend.text = ggplot2::element_text(size = 7)) +
        ggplot2::theme(
            axis.title.x = ggplot2::element_text(vjust = -0.35),
            axis.title.y = ggplot2::element_text(vjust = 0.35),
            text = ggplot2::element_text(size = 8),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
        ) +
        # ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        # ggplot2::theme(axis.ticks.x = ggplot2::element_blank()) +
        ggplot2::labs(x = "Subject", y = "Percentage")

    ggplot2::ggsave(here::here("doc/figures/figure_1/fiber_type_participant_proteomics.png"),
                    units = "mm",
                    width = 60,
                    height = 60)
