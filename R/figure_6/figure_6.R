################################################################################################################################################
#################################################     Panel B  ##############################################################
################################################################################################################################################

metadata_MD <- vroom::vroom(here::here("data/metadata_MD_w_fiber_type_w_anonim.csv"))

# overall % of fiber types ------------------------------------------------

metadata_MD |>
    ggplot2::ggplot(
        ggplot2::aes(x = fiber_type,
                     fill = fiber_type)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::scale_fill_manual(values = c(
        "#3B528BFF",
        "#fdc325",
        "#440154FF",
        "#5DC863FF")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~ condition + subject)



MD_ft <- data.frame(
    metadata_MD |>
        dplyr::filter(subject == "A1") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("A1" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(A1),
    metadata_MD |>
        dplyr::filter(subject == "A2") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("A2" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(A2),
    metadata_MD |>
        dplyr::filter(subject == "A3") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("A3" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(A3),
    metadata_MD |>
        dplyr::filter(subject == "C1") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("C1" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(C1),
    metadata_MD |>
        dplyr::filter(subject == "C2") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("C2" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(C2),
    metadata_MD |>
        dplyr::filter(subject == "C3") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("C3" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(C3),
    metadata_MD |>
        dplyr::filter(subject == "T1") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("T1" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(T1),
    metadata_MD |>
        dplyr::filter(subject == "T2") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        t() |>
        as.data.frame() |>
        tibble::add_column(V4 = c("Hybrid 2A/2X", 0)) |>
        t() |>
        as.data.frame() |>
        dplyr::mutate(cnt = as.numeric(cnt)) |>
        dplyr::mutate("T2" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::remove_rownames() |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(T2),
    metadata_MD |>
        dplyr::filter(subject == "T3") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        t() |>
        as.data.frame() |>
        tibble::add_column(V4 = c("Hybrid 2A/2X", 0)) |>
        t() |>
        as.data.frame() |>
        dplyr::mutate(cnt = as.numeric(cnt)) |>
        dplyr::mutate("T3" = round(cnt / sum(cnt), 3) * 100) |>
        dplyr::arrange(desc(fiber_type)) |>
        tibble::remove_rownames() |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(T3)
) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("subject") |>
    tidyr::pivot_longer(
        cols = c("Hybrid 1/2A", "Hybrid 2A/2X", "Type 1", "Type 2A"),
        names_to = "fiber_type"
    )


MD_ft$fiber_type <- factor(MD_ft$fiber_type, levels = c(
    "Type 1",
    "Hybrid 1/2A",
    "Type 2A",
    "Hybrid 2A/2X"
))

MD_ft |>
    ggplot2::ggplot(ggplot2::aes(
        x = subject,
        y = value,
        fill = fiber_type
    )) +
    ggplot2::geom_col(
        alpha = 0.85,
        colour = NA,
        size = 0
    ) +
    # ggplot2::facet_grid(~subject) +
    ggplot2::scale_fill_manual("Fiber \nTypes",
                               values = c("#440154FF", "#3B528BFF", "#5DC863FF", "#fdc325"),
                               labels = c("Type 1", "Hybrid\n 1/2A", "Type 2A", "Hybrid\n 2A/2X")
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Fiber type by participant") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(
        legend.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        ),
        legend.key.size = ggplot2::unit(3, "mm"),
        legend.text = ggplot2::element_text(size = 7),
        legend.position = "right"
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(
        expand = c(0, 0)
    ) +
    ggplot2::theme(
        axis.title.x = ggplot2::element_text(vjust = -0.35),
        axis.title.y = ggplot2::element_text(vjust = 0.35),
        text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7)) +
    # ggplot2::theme(axis.ticks.y = ggplot2::element_blank()) +
    # ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
    ggplot2::labs(x = "Subject", y = "Percentage") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        strip.text = ggplot2::element_text(colour = "white"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent', colour = "transparent"), #transparent legend panel
        legend.position = "right"
    ) +
    ggplot2::annotate(
        geom = "rect",
        xmin = 0.56,
        xmax = 3.45,
        ymin = 101,
        ymax = 110,
        fill = "#E48C2AFF",
        alpha = 0.5
    ) +
    ggplot2::annotate(
        geom = "rect",
        xmin = 3.55,
        xmax = 6.45,
        ymin = 101,
        ymax = 110,
        fill = "#969594",
        alpha = 0.5
    ) +
    ggplot2::annotate(
        geom = "rect",
        xmin = 6.55,
        xmax = 9.45,
        ymin = 101,
        ymax = 110,
        fill = "#d662c4",
        alpha = 0.5
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 2.005,
        y = 106.5,
        label = "ACTA1-NM",
        fontface = "bold",
        size = 2.5
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 5,
        y = 106.5,
        label = "Control",
        fontface = "bold",
        size = 2.5
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 8,
        y = 106.5,
        label = "TNNT1-NM",
        fontface = "bold",
        size = 2.5
    ) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,2,0), "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
                   legend.box.margin= ggplot2::margin(t = -10,b = -10,2,-5),
                   legend.title = ggplot2::element_text(size = 7, face = "bold"),
                   legend.text = ggplot2::element_text(size = 5),
                   legend.key.size = ggplot2::unit(2, "mm"),
                   legend.spacing.x = ggplot2::unit(2, "mm"))

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6B.png"),
                units = "mm",
                height = 60,
                width = 70)
