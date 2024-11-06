################################################################################################################################################
#################################################     Panel A  ##############################################################
################################################################################################################################################
muscle_disease_data <- vroom::vroom(
    here::here("data/data_muscle_disease.csv")
) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene_name")

metadata_muscle_disease <- vroom::vroom(
    here::here("data/metadata_muscle_disease.csv")
) |>
    dplyr::select(!1) |>
    dplyr::mutate(subject = dplyr::case_when(
        subject == "T1" ~ "TM1",
        subject == "T2" ~ "TM2",
        subject == "T3" ~ "TM3",
        subject == "A1" ~ "AM1",
        subject == "A2" ~ "AM2",
        subject == "A3" ~ "AM3",
        subject == "C1" ~ "C1",
        subject == "C2" ~ "C2",
        subject == "C3" ~ "C3",
        TRUE ~ "error"
    ))

muscle_disease_number_quantified_proteins <- colSums(!is.na(muscle_disease_data)) |>
    as.data.frame() |>
    tibble::rownames_to_column("SampleID")

colnames(muscle_disease_number_quantified_proteins) <- c(
    "fiber_ID",
    "number_of_proteins"
)

muscle_disease_number_quantified_proteins <- metadata_muscle_disease |>
    dplyr::select(fiber_ID, condition, subject) |>
    dplyr::inner_join(muscle_disease_number_quantified_proteins) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::rename("study" = "condition")

muscle_disease_number_quantified_proteins |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = subject,
            y = number_of_proteins,
            fill = study
        )
    ) +
    ggplot2::geom_violin(alpha = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual("condition", values = c("#E48C2AFF",
                                                       "#969594",
                                                       "#d662c4")) +
    ggplot2::ylim(1000, 2500) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Nemaline Myopathy",
                  y = "Number of quantified proteins") +
    # ggplot2::scale_x_discrete(breaks = c("heterofiber"),
    #                           labels = c("")) +
    ggplot2::theme(text = ggplot2::element_text(size = 9.5)) +
    ggplot2::theme(
        legend.key.size = ggplot2::unit(4, "mm")
    )

# ggplot2::ggsave(here::here("doc/figures/figure_6_S1/figure_6_S1_A.png"),
#                 height = 60,
#                 width = 128,
#                 units = "mm")

################################################################################################################################################
#################################################     Panel B  ##############################################################
################################################################################################################################################

data_controls <- muscle_disease_data |>
    dplyr::select(dplyr::starts_with("c"))

sum_of_intensities <- colSums(data_controls,
                              na.rm = TRUE)

rel_abundance <- data_controls |>
    t() |>
    as.data.frame()

rel_abundance <- rel_abundance/sum_of_intensities * 100

mean_rel_abundance <- rel_abundance |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        mean,
        na.rm = TRUE
    )) |>
    dplyr::slice_head(n = 1) |>
    t() |>
    as.data.frame() |>
    log10()


colnames(mean_rel_abundance) <- "log10_rel_abundance"

Gene_names <- rownames(mean_rel_abundance)

mean_expression_all <- mean_rel_abundance |>
    tibble::rownames_to_column("gene")

mean_expression_all <- mean_expression_all |>
    dplyr::arrange(desc(log10_rel_abundance)) |>
    dplyr::mutate(
        order = 1:nrow(data_controls)
    )

library(ggplot2)
library(tidyverse)
library(ggrepel)

ggplot2::ggplot() +

    # Add all genes
    ggplot2::geom_point(data = mean_expression_all,
                        aes(x=order,
                            y=log10_rel_abundance),
                        colour = "#756bb1",
                        size = 0.25,
                        alpha = 0.5) +

    # Add horizontal lines to indicate % instead of log scale
    geom_hline(yintercept = log(20,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(5,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(1,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(0.1,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(0.01,10), linetype="dashed", linewidth=0.2) +

    # Add text for % expression
    annotate("text", x=2900, y=log(28,10), label= "20%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(7,10), label= "5%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(1.4,10), label= "1%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(0.14,10), label= "0.1%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(0.014,10), label= "0.01%", colour="black", fontface=2, size=2) +

    # Add custom labels:
    geom_label_repel(data = mean_expression_all %>% dplyr::filter(gene %in% c("MYH7",
                                                                              "MYH2",
                                                                              "ACTA1",
                                                                              "TNNT1",
                                                                              "TNNT3")),
                     mapping = aes(x = order,
                                   y = log10_rel_abundance,
                                   label = gene),
                     size = 1.8, label.padding=0.1, max.overlaps = Inf, min.segment.length=0.1, segment.size=0.2, force = 10) +
    # Change design
    ylab("% total intensities, 10log") +
    xlab("Protein rank") +
    theme_classic() +
    theme(
        text = element_text(face="bold", colour="black", size = 6),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        legend.position = "none",
    )

ggplot2::ggsave(here::here("doc/figures/figure_6_S1/figure_6_S1_B.png"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
#################################################     Panel C  ##############################################################
################################################################################################################################################

muscle_disease_data_controls <- vroom::vroom(
    here::here("data/data_muscle_disease.csv")
) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene_name") |>
    dplyr::select(dplyr::starts_with("c")) |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        log2
    ))

proteomics_data <- vroom::vroom(
    here::here("data/data_proteomics_filtered.csv"),
    col_select = !c(1)
) |>
    as.data.frame() |>
    dplyr::filter(
        !Gene.name == "",
        !duplicated(Gene.name)
    ) |>
    tibble::column_to_rownames("Gene.name") |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        log2
    ))

avg_intensities_MD <- muscle_disease_data_controls |>
    t() |>
    as.data.frame() |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        mean, na.rm = TRUE
    )) |>
    dplyr::slice_head(n = 1) |>
    t() |>
    as.data.frame() |>
    dplyr::rename("avg_intensities_MD" = 1) |>
    tibble::rownames_to_column("Genes")

avg_intensities_heterofiber <- proteomics_data |>
    t() |>
    as.data.frame() |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        mean, na.rm = TRUE
    )) |>
    dplyr::slice_head(n = 1) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% avg_intensities_MD$Genes) |>
    dplyr::rename(
        "avg_intensities_heterofiber" = 2
    )

avg_intensities_MD <- avg_intensities_MD |>
    dplyr::filter(Genes %in% avg_intensities_heterofiber$Genes)

joined_intensities <- avg_intensities_heterofiber |>
    dplyr::inner_join(avg_intensities_MD) |>
    dplyr::filter(!is.na(avg_intensities_MD)) |>
    dplyr::filter(!is.na(avg_intensities_heterofiber))

pearson_correlation <- cor.test(x = joined_intensities$avg_intensities_heterofiber,
                                y = joined_intensities$avg_intensities_MD)

R_pearson <- round(pearson_correlation$estimate, 3)

P_pearson <- formatC(pearson_correlation$p.value, format = "e", digits = 2)

joined_intensities |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = avg_intensities_heterofiber,
            y = avg_intensities_MD,
            names = Genes
        )
    ) +
    ggplot2::geom_point(alpha = 0.65,
                        size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::annotate(
        geom = "text",
        x = 10,
        y = 20,
        label = paste("r =", R_pearson),
        colour = "black",
        size = 2.5
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 10,
        y = 19,
        label = "p < 0.001",
        colour = "black",
        size = 2.5
    ) +
    ggplot2::xlab("Avg LFQ intensities 1000 fibers (Log2)") +
    ggplot2::ylab("Avg LFQ intensities muscle disease (Log2)") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7)
    )

ggplot2::ggsave(here::here("doc/figures/figure_6_S1/figure_6_S1_C.png"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
#################################################     Panel D  ##############################################################
################################################################################################################################################


# Sorting dataset with the genes that are used for fiber typing,
# plotting intensities and rank:
data_fiber_typing <- muscle_disease_data |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2", "MYH1") |>
    dplyr::mutate("Fiber_ID" = colnames(muscle_disease_data)) |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
    dplyr::arrange(desc(MYH7)) |>
    dplyr::mutate("order_7" = seq_len(ncol(muscle_disease_data))) |>
    dplyr::arrange(desc(MYH2)) |>
    dplyr::mutate("order_2" = seq_len(ncol(muscle_disease_data)))

data_fiber_typing |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_7,
            y = MYH7
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal()


sum_of_intensities <- muscle_disease_data |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2", "MYH1") |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
    tidyr::replace_na(list(
        "MYH7" = 0,
        "MYH2" = 0,
        "MYH1" = 0
    )) |>
    rowSums()

perc_MYHs <- muscle_disease_data |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2", "MYH1") |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
    tidyr::replace_na(list(
        "MYH7" = 0,
        "MYH2" = 0,
        "MYH1" = 0
    )) |>
    tibble::add_column(sum_of_intensities) |>
    dplyr::mutate(across(
        .cols = c("MYH7", "MYH2", "MYH1"),
        ~ .x / sum_of_intensities * 100
    )) |>
    as.data.frame() |>
    dplyr::arrange(desc(MYH7)) |>
    dplyr::mutate("order_7" = seq_len(ncol(muscle_disease_data))) |>
    dplyr::arrange(desc(MYH2)) |>
    dplyr::mutate("order_2" = seq_len(ncol(muscle_disease_data))) |>
    dplyr::arrange(desc(MYH1)) |>
    dplyr::mutate("order_1" = seq_len(ncol(muscle_disease_data))) |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::arrange(desc(fiber_ID)) |>
    dplyr::select(c(
        "MYH7",
        "MYH2",
        "MYH1",
        "order_7",
        "order_2",
        "order_1",
        "fiber_ID"
    )) |>
    tidyr::pivot_longer(
        cols = c("MYH7", "MYH2", "MYH1"),
        names_to = "MYHs",
        values_to = "values"
    )


# Individual curve for intensity of MYH7:
perc_MYHs |>
    dplyr::filter(MYHs == "MYH7") |>
    ggplot2::ggplot(
        ggplot2::aes(order_7,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::theme_minimal()


# MYH7 curve --------------------------------------------------------------

# Finding bottom knee of MYH 7 curve,
# that's the one we are using as a threshold.
# The function detects top knee so we do the inverse of the curve:

MYH_7_curve <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH7) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    t()

MYH_7_curve <- DropletUtils::barcodeRanks(MYH_7_curve, lower = 10)

bottom_knee_MYH7 <- 100 - MYH_7_curve@metadata$knee

perc_MYHs |>
    dplyr::filter(MYHs == "MYH7") |>
    ggplot2::ggplot(
        ggplot2::aes(order_7,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = "#54278f") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH7,
                        color = "#54278f",
                        size = 0.6,
                        alpha = 0.5) +
    ggplot2::xlab("Fiber rank") +
    ggplot2::ylab("MYH7 fraction (%)") +
    ggplot2::annotate(
        geom = "text",
        x = 65,
        y = 44,
        label = "MYH7 = 38%",
        size = 2,
        color = "#54278f",
        fontface = "bold"
    ) +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6.25,
                                     face = "bold")
    )

ggplot2::ggsave(here::here(
    "doc/figures/figure_6_S1/figure_6_S1_D_1.png"
),
units = "mm",
height = 40,
width = 40)

MYH_2_curve <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH2) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    t()

MYH_2_curve <- DropletUtils::barcodeRanks(MYH_2_curve, lower = 10)

bottom_knee_MYH2 <- 100 - MYH_2_curve@metadata$knee

perc_MYHs |>
    dplyr::filter(MYHs == "MYH2") |>
    ggplot2::ggplot(
        ggplot2::aes(order_2,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = "#54278f") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH2,
                        color = "#54278f",
                        size = 0.6,
                        alpha = 0.5) +
    ggplot2::xlab("Fiber rank") +
    ggplot2::ylab("MYH2 fraction (%)") +
    ggplot2::annotate(
        geom = "text",
        x = 65,
        y = 16,
        label = "MYH2 = 10%",
        size = 2,
        color = "#54278f",
        fontface = "bold"
    ) +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6.25,
                                     face = "bold")
    )

ggplot2::ggsave(here::here(
    "doc/figures/figure_6_S1/figure_6_S1_D_2.png"
),
units = "mm",
height = 40,
width = 40)


# MYH1 curve -------------------------------------------------------------------

MYH_1_curve <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH1) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    t()

MYH_1_curve <- DropletUtils::barcodeRanks(MYH_1_curve, lower = 10)

bottom_knee_MYH1 <- 100 - MYH_1_curve@metadata$knee

perc_MYHs |>
    dplyr::filter(MYHs == "MYH1") |>
    ggplot2::ggplot(
        ggplot2::aes(order_1,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = "#54278f") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH1,
                        color = "#54278f",
                        size = 0.6,
                        alpha = 0.5) +
    ggplot2::xlab("Fiber rank") +
    ggplot2::ylab("MYH1 fraction (%)") +
    ggplot2::annotate(
        geom = "text",
        x = 150,
        y = 22,
        label = "MYH1 = 17%",
        size = 2,
        color = "#54278f",
        fontface = "bold"
    ) +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6.25,
                                     face = "bold")
    )

ggplot2::ggsave(here::here(
    "doc/figures/figure_6_S1/figure_6_S1_D_3.png"
),
units = "mm",
height = 40,
width = 40)

################################################################################################################################################
#################################################     Panel E  ##############################################################
################################################################################################################################################

data_fiber_type <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    dplyr::mutate(
        fiber_type = dplyr::case_when(
            MYH7 >= bottom_knee_MYH7 &
                MYH2 <= bottom_knee_MYH2 ~ "Type 1",
            MYH7 >= bottom_knee_MYH7 &
                MYH2 >= bottom_knee_MYH2 ~ "Hybrid 1/2A",
            MYH1 >= bottom_knee_MYH1 &
                MYH2 >= bottom_knee_MYH2 ~ "Hybrid 2A/2X",
            MYH2 >= bottom_knee_MYH2 &
                MYH7 <= bottom_knee_MYH7 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::arrange(desc(fiber_ID))

# readr::write_csv(data_fiber_type,
#                  here::here("data/MD_perc_MYHs.csv"))

# Complete fiber type curves ----------------------------------------------

# Here I'm basing my final plot of fiber typing curves on Thibaux script:

fiber_types <- list(
    type_1 <- data_fiber_type |>
        dplyr::filter(fiber_type == "Type 1") |>
        dplyr::arrange(desc(MYH7)),
    Hybrid_1_2A <- data_fiber_type |>
        dplyr::filter(fiber_type == "Hybrid 1/2A") |>
        dplyr::arrange(desc(MYH7)),
    type_2A <- data_fiber_type |>
        dplyr::filter(fiber_type == "Type 2A") |>
        dplyr::arrange(MYH2),
    Hybrid_2A_2X <- data_fiber_type |>
        dplyr::filter(fiber_type == "Hybrid 2A/2X") |>
        dplyr::arrange(MYH1)
)
names(fiber_types) <- c(
    "type_1",
    "hybrid_1_2A",
    "type_2A",
    "hybrid_2A_2X"
)

order_matrix <- dplyr::bind_rows(
    fiber_types$type_1,
    fiber_types$hybrid_1_2A,
    fiber_types$type_2A,
    fiber_types$hybrid_2A_2X
) |>
    tibble::add_column("sample_MYH" = seq_len(nrow(data_fiber_type)))

rectangle_colors <- data.frame(
    start = c(0, 96, 159, 242),
    end = c(96, 159, 242, 272),
    fiber_type = c("Type I", "Hybrid I/IIA", "Type IIA", "Hybrid IIA/IIX")
)

order_matrix |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
        data = rectangle_colors,
        ggplot2::aes(
            xmin = start,
            xmax = end,
            ymin = -Inf,
            ymax = Inf,
            fill = fiber_type
        ),
        alpha = 0.15
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#3B528BFF",
        "#fdc325",
        "#440154FF",
        "#5DC863FF"
    )) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_point(
        ggplot2::aes(
            x = sample_MYH,
            y = MYH7
        ),
        colour = "#440154FF",
        size = 0.25
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
            x = sample_MYH,
            y = MYH2
        ),
        colour = "#5DC863FF",
        size = 0.25
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
            x = sample_MYH,
            y = MYH1
        ),
        colour = "#fdc325",
        size = 0.25
    ) +
    ggplot2::ylab("% MYH isoform expressed") +
    ggplot2::xlab("Sample (ranked by MYH expression)") +
    ggplot2::ggtitle("MYH-based fiber typing (proteomics)") +
    ggplot2::geom_vline(xintercept = 96, colour = "grey30", size = 0.25) +
    ggplot2::geom_vline(xintercept = 159, colour = "grey30", size = 0.25) +
    ggplot2::geom_vline(xintercept = 242, colour = "grey30", size = 0.25) +
    ggplot2::geom_hline(
        yintercept = 38,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = 10,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = 17,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
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
    ggplot2::annotate(
        "text",
        x = 48,
        y = 115,
        label = "Type 1 \n 96 (35.3%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 199,
        y = 115,
        label = "Type 2A \n 82 (30.1%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 256.5,
        y = 110,
        label = "Hybrid \n2A/2X \n31 (11.4%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 127,
        y = 110,
        label = "Hybrid \n1/2A \n63 (23.3%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 15,
        y = 95,
        label = "MYH7",
        colour = "#440154FF",
        fontface = 2,
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 232,
        y = 80,
        label = "MYH2",
        colour = "#5DC863FF",
        fontface = 2,
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 265,
        y = 67,
        label = "MYH1",
        colour = "#fdc325",
        fontface = 2,
        size = 2
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 272), clip = "off") +
    ggplot2::annotate(
        "text",
        x = 20,
        y = 43,
        label = "MYH7 threshold",
        colour = "black",
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 20,
        y = 6.5,
        label = "MYH2 threshold",
        colour = "black",
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 20,
        y = 22,
        label = "MYH1 threshold",
        colour = "black",
        size = 2
    ) +
    ggplot2::theme(text = ggplot2::element_text(size = 7)) +
    ggplot2::scale_x_continuous(limits = c(0, 272),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0,125),
                                expand = c(0, 0),
                                breaks = c(0, 25, 50, 75, 100),
                                labels = c("0", "25", "50", "75", "100")) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 7))

ggplot2::ggsave(
    here::here("doc/figures/figure_6_S1/figure_6_S1_E.png"),
    device = "png",
    width = 128,
    height = 60,
    units = "mm"
)

metadata_fiber_type <- data_fiber_type |>
    dplyr::select(fiber_ID, fiber_type) |>
    dplyr::inner_join(metadata_muscle_disease)

readr::write_csv(metadata_fiber_type,
                 here::here("data/metadata_MD_w_fiber_type_w_anonim.csv"))
