library(ggplot2)
library(dplyr)
# Loading transcriptomics and making transcriptomics plot -----------------
counts_ft <- readr::read_csv(file = here::here("data-raw/counts_ft.csv"))

data_plot_fibertype_participant <- counts_ft %>%
    dplyr::mutate(
        fiber_type_MYH = dplyr::case_when(
            fiber_type_MYH == "Type 1" ~ "Type 1",
            fiber_type_MYH == "Type 2A" ~ "Type 2A",
            fiber_type_MYH == "Type 2X" ~ "Type 2X",
            fiber_type_MYH == "Hybrid 1/2A" ~ "Hybrid 1/2A",
            fiber_type_MYH == "Hybrid 1/2X" ~ "Type 1",
            fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction >= MYH1_fraction ~ "Hybrid 1/2A",
            fiber_type_MYH == "Hybrid 2A/2X" ~ "Hybrid 2A/2X",
            fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction < MYH1_fraction ~ "Hybrid 2A/2X",
            TRUE ~ "NA"
        ))

data_plot_fibertype_participant$fiber_type_MYH <- factor(data_plot_fibertype_participant$fiber_type_MYH, levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))

transcriptomics_plot <- data_plot_fibertype_participant %>%
    ggplot(aes(
        x = subject,
        fill = fiber_type_MYH
    )) +
    geom_bar(na.rm = TRUE, position="fill", alpha = 0.85, colour=NA, size=0) +
    scale_fill_manual("Fiber type", values = c("#440154FF", "#3B528BFF", "#5DC863FF", "#fdc325", "#D2631C")) +
    labs(
        x = "Subject",
        y = "Percentage"
    ) +
    theme_classic() +
    ggtitle("MYH-based fiber typing by participant (transcriptomics)\n") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = c("0", "25", "50", "75", "100")) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), labels = c("0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14")) +
    theme(
        text = element_text(face="bold", colour="black", size=8),
        strip.text = element_text(colour = "white"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent', colour = "transparent") #transparent legend panel
        # legend.position = "none",
    )

# Loading proteomics and making plot --------------------------------------

metadata_fiber_type <- readr::read_csv(file=here::here("data/metadata_proteomics.csv"))

proteomics_ft <- data.frame(
    metadata_fiber_type |>
        dplyr::filter(subject == "P1") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("P1" = round(cnt / sum(cnt), 3) * 100) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(P1),
    metadata_fiber_type |>
        dplyr::filter(subject == "P2") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("P2" = round(cnt / sum(cnt), 3) * 100) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(P2),
    metadata_fiber_type |>
        dplyr::filter(subject == "P3") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("P3" = round(cnt / sum(cnt), 3) * 100) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(P3),
    metadata_fiber_type |>
        dplyr::filter(subject == "P4") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("P4" = round(cnt / sum(cnt), 3) * 100) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(P4),
    metadata_fiber_type |>
        dplyr::filter(subject == "P5") |>
        dplyr::group_by(fiber_type) |>
        dplyr::summarise(cnt = dplyr::n()) |>
        dplyr::mutate("P5" = round(cnt / sum(cnt), 3) * 100) |>
        tibble::column_to_rownames("fiber_type") |>
        dplyr::select(P5)
) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("subject") |>
    tidyr::pivot_longer(
        cols = c("Hybrid 1/2A", "Hybrid 2A/2X", "Type 1", "Type 2A"),
        names_to = "fiber_type"
    )

proteomics_ft$fiber_type<- factor(proteomics_ft$fiber_type, levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X"))
proteomics_ft$subject <- factor(proteomics_ft$subject, levels = c("P1",
                                                                  "P2",
                                                                  "P3",
                                                                  "P4",
                                                                  "P5"))

proteomics_plot <- proteomics_ft |>
    ggplot2::ggplot(ggplot2::aes(x = subject,
                                 y = value,
                                 fill = fiber_type)) +
    ggplot2::geom_col(alpha = 0.85,
                      colour = NA,
                      size = 0) +
    # ggplot2::facet_grid(~subject) +
    ggplot2::scale_fill_manual("Fiber types",
                               values = c("#440154FF", "#3B528BFF", "#5DC863FF", "#fdc325")) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("MYH-based fiber typing \nby participant (proteomics)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(
        legend.title = ggplot2::element_text(size = 8,
                                             face = "bold"),
        legend.key.size = ggplot2::unit(3, "mm"),
        legend.text = ggplot2::element_text(size = 7),
        legend.position = "none"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(
        expand = c(0, 0),
        labels = c(
            "FOR2" = "P1",
            "FOR4" = "P2",
            "FOR9" = "P3",
            "FOR10" = "P4",
            "FOR11" = "P5"
        )
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
    theme(
        text = element_text(
            face = "bold",
            colour = "black",
            size = 8
        ),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill = "black"),
        panel.background = element_rect(fill = 'transparent'),
        #transparent panel bg
        plot.background = element_rect(fill = 'transparent', color = NA),
        #transparent plot bg
        panel.grid.major = element_blank(),
        #remove major gridlines
        panel.grid.minor = element_blank(),
        #remove minor gridlines
        legend.background = element_rect(fill = 'transparent'),
        #transparent legend bg
        legend.box.background = element_rect(fill = 'transparent') #transparent legend panel
    )

################################################################################################################################################
########################################################      Panel A   ############################################################################
################################################################################################################################################


combined_plot <- ggpubr::ggarrange(transcriptomics_plot,
                                   proteomics_plot,
                                   common.legend = TRUE,
                                   legend = "right",
                                   widths = c(2.5,1)
)

# ggsave(combined_plot,
#        filename = here::here("doc/figures/figure_1_S2/fiber_type_participant_combined.png"),
#        width = 195,
#        height = 60,
#        units = "mm")

# Exporting MYH abundances for dot blot-selected fibers -------------------

dot_blot_fibers <- c(
    "P3_fibNumber33",
    "P5_fibNumber175",
    "P3_fibNumber180",
    "P4_fibNumber162",
    "P1_fibNumber192",
    "P2_fibNumber36",
    "P5_fibNumber159",
    "P2_fibNumber90",
    "P2_fibNumber164",
    "P1_fibNumber177",
    "P5_fibNumber35",
    "P5_fibNumber187",
    "P5_fibNumber128",
    "P5_fibNumber78",
    "P4_fibNumber136",
    "P4_fibNumber159",
    "P4_fibNumber188",
    "P4_fibNumber60",
    "P4_fibNumber24",
    "P3_fibNumber108",
    "P4_fibNumber108",
    "P3_fibNumber110",
    "P4_fibNumber182",
    "P2_fibNumber29",
    "P4_fibNumber152",
    "P3_fibNumber177",
    "P2_fibNumber194",
    "P4_fibNumber199",
    "P3_fibNumber88",
    "P3_fibNumber193",
    "P4_fibNumber57",
    "P3_fibNumber100",
    "P1_fibNumber100",
    "P3_fibNumber185",
    "P3_fibNumber80",
    "P1_fibNumber171",
    "P1_fibNumber198",
    "P5_fibNumber199",
    "P4_fibNumber181",
    "P5_fibNumber119"
)

data_dot_blot <- data_fiber_type  |>
    dplyr::filter(fiber_ID %in% dot_blot_fibers)  |>
    dplyr::arrange(match(fiber_ID, dot_blot_fibers)) |>
    dplyr::select(fiber_ID, MYH7, MYH2, MYH1, fiber_type) |>
    dplyr::mutate("sample_order_left_to_right" = rep(c(1:10), 4))

readr::write_csv(data_dot_blot,
                 here::here("data/dot_blot/MS_fiber_type_dot_blot_fibers.csv"))
