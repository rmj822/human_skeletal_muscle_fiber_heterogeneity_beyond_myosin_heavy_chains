
data_proteomics <- vroom::vroom(here::here("data/proteomics_fiber_type_markers.csv")) |>
  dplyr::select(!...1)


# Calculating % of MYH using only 7 and 2 ---------------------------------

proteomics_MYH_filtered <- vroom::vroom(
    here::here("data/proteomics_MYH_filtered.csv"),
    col_select = !c(1)
) |>
    as.data.frame() |>
    tibble::column_to_rownames("Gene.name")

sum_of_intensities <- proteomics_MYH_filtered |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2") |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2"), as.numeric)) |>
    tidyr::replace_na(list(
        "MYH7" = 0,
        "MYH2" = 0
    )) |>
    rowSums()


perc_MYHs <- proteomics_MYH_filtered |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2") |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2"), as.numeric)) |>
    tidyr::replace_na(list(
        "MYH7" = 0,
        "MYH2" = 0
    )) |>
    tibble::add_column(sum_of_intensities) |>
    dplyr::mutate(across(
        .cols = c("MYH7", "MYH2"),
        ~ .x / sum_of_intensities * 100
    )) |>
    as.data.frame() |>
    dplyr::arrange(desc(MYH7)) |>
    dplyr::mutate("order_7" = seq_len(ncol(proteomics_MYH_filtered))) |>
    dplyr::arrange(desc(MYH2)) |>
    dplyr::mutate("order_2" = seq_len(ncol(proteomics_MYH_filtered))) |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::arrange(desc(fiber_ID)) |>
    dplyr::select(c(
        "MYH7",
        "MYH2",
    "order_7",
        "order_2",
        "fiber_ID"
    ))

# troponin T fiber typing -------------------------------------------------

data_troponins <- data_proteomics |>
  dplyr::filter(Gene.name %in% c(
    "TNNT1",
    "TNNT3"
  )) |>
  tibble::column_to_rownames("Gene.name") |>
  dplyr::mutate(dplyr::across(.cols = everything(), as.numeric)) |>
  t() |>
  as.data.frame() |>
    tidyr::replace_na(replace = list(
        TNNT1 = 0,
        TNNT3 = 0
    ))

sum_of_intensities <- data_troponins |>
  t() |>
  colSums()

perc_troponins <- data_troponins |>
  dplyr::mutate(across(
    .cols = c("TNNT1", "TNNT3"),
    ~ .x / sum_of_intensities * 100
  )) |>
  as.data.frame() |>
  dplyr::arrange(desc(TNNT1)) |>
  dplyr::mutate("order_TNNT1" = seq_len(nrow(perc_MYHs))) |>
  dplyr::arrange(desc(TNNT3)) |>
  dplyr::mutate("order_TNNT3" = seq_len(nrow(perc_MYHs))) |>
  tibble::rownames_to_column("fiber_ID")

# Curve fitting for TNNT1 with MYH7

perc_troponins |>
  dplyr::inner_join(perc_MYHs) |>
  dplyr::arrange(desc(order_7)) |>
    dplyr::select(MYH7, order_7, TNNT1, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH7, TNNT1),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = order_7,
      y = perc_of_isoform_expression,
      color = Genes
    )
  ) +
  ggplot2::geom_point(
    alpha = 0.5
  ) +
  ggplot2::scale_color_manual(values = c("#440154FF",
                                                    "darkblue")) +
  ggplot2::theme_minimal()

# ggplot2::ggsave(here::here("doc/figures/figure_2/plot_MYH7_TNNT1.png"))

# Curve fitting for TNNT3 with MYH2

perc_troponins |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_2)) |>
    dplyr::select(MYH2, order_2, TNNT3, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH2, TNNT3),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_2,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("#5DC863FF",
                                                      "darkred")) +
                                                          ggplot2::theme_minimal()

# ggplot2::ggsave(here::here("doc/figures/figure_2/plot_MYH2_TNNT3.png"))

# ATP2A fiber typing ------------------------------------------------------

data_ATP2A <- data_proteomics |>
    dplyr::filter(Gene.name %in% c(
        "ATP2A1",
        "ATP2A2"
    )) |>
    tibble::column_to_rownames("Gene.name") |>
    dplyr::mutate(dplyr::across(.cols = everything(), as.numeric)) |>
    t() |>
    as.data.frame() |>
    tidyr::replace_na(replace = list(
        ATP2A1 = 0,
        ATP2A2 = 0
    ))

sum_of_intensities <- data_ATP2A |>
    t() |>
    colSums()

perc_ATP2A <- data_ATP2A |>
    dplyr::mutate(across(
        .cols = c("ATP2A1", "ATP2A2"),
        ~ .x / sum_of_intensities * 100
    )) |>
    as.data.frame() |>
    dplyr::arrange(desc(ATP2A1)) |>
    dplyr::mutate("order_ATP2A1" = seq_len(nrow(perc_MYHs))) |>
    dplyr::arrange(desc(ATP2A2)) |>
    dplyr::mutate("order_ATP2A2" = seq_len(nrow(perc_MYHs))) |>
    tibble::rownames_to_column("fiber_ID")

# Curve fitting for ATP2A2 with MYH7

perc_ATP2A |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_7)) |>
    dplyr::select(MYH7, order_7, ATP2A2, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH7, ATP2A2),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_7,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("darkblue", "#440154FF")) +
                                                          ggplot2::theme_minimal()

# ggplot2::ggsave(here::here("doc/figures/figure_2/plot_MYH7_ATP2A2.png"))

# Curve fitting for ATP2A1 with MYH2

perc_ATP2A |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_2)) |>
    dplyr::select(MYH2, order_2, ATP2A1, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH2, ATP2A1),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_2,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("darkred", "#5DC863FF")) +
                                                          ggplot2::theme_minimal()

# ggplot2::ggsave(here::here("doc/figures/figure_2/plot_MYH2_ATP2A1.png"))

# MYL fiber typing --------------------------------------------------------

data_MYL <- data_proteomics |>
    dplyr::filter(Gene.name %in% c(
        "MYL1",
        "MYL3"
    )) |>
    tibble::column_to_rownames("Gene.name") |>
    dplyr::mutate(dplyr::across(.cols = everything(), as.numeric)) |>
    t() |>
    as.data.frame() |>
    tidyr::replace_na(replace = list(
        MYL1 = 0,
        MYL3 = 0
    ))

sum_of_intensities <- data_MYL |>
    t() |>
    colSums()

perc_MYL <- data_MYL |>
    dplyr::mutate(across(
        .cols = c("MYL1", "MYL3"),
        ~ .x / sum_of_intensities * 100
    )) |>
    as.data.frame() |>
    dplyr::arrange(desc(MYL1)) |>
    dplyr::mutate("order_MYL1" = seq_len(nrow(perc_MYHs))) |>
    dplyr::arrange(desc(MYL3)) |>
    dplyr::mutate("order_MYL3" = seq_len(nrow(perc_MYHs))) |>
    tibble::rownames_to_column("fiber_ID")

# Curve fitting for MYL3 with MYH7

perc_MYL |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_7)) |>
    dplyr::select(MYH7, order_7, MYL3, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH7, MYL3),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_7,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("#440154FF",
                                                      "darkblue")) +
                                                          ggplot2::theme_minimal()

# ggplot2::ggsave(here::here("doc/figures/figure_2/plot_MYH7_MYL3.png"))

# Curve fitting for MYL1 with MYH2

perc_MYL |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_2)) |>
    dplyr::select(MYH2, order_2, MYL1, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH2, MYL1),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_2,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("#5DC863FF",
                                                      "darkred")) +
                                                          ggplot2::theme_minimal()

# ggplot2::ggsave(here::here("doc/figures/figure_2/plot_MYH2_MYL1.png"))


# Tropomyosins ------------------------------------------------------------

data_TPM <- data_proteomics |>
    dplyr::filter(Gene.name %in% c(
        "TPM1",
        "TPM3"
    )) |>
    tibble::column_to_rownames("Gene.name") |>
    dplyr::mutate(dplyr::across(.cols = everything(), as.numeric)) |>
    t() |>
    as.data.frame() |>
    tidyr::replace_na(replace = list(
        TPM1 = 0,
        TPM3 = 0
    ))

sum_of_intensities <- data_TPM |>
    t() |>
    colSums()

perc_TPM <- data_TPM |>
    dplyr::mutate(across(
        .cols = c("TPM1", "TPM3"),
        ~ .x / sum_of_intensities * 100
    )) |>
    as.data.frame() |>
    dplyr::arrange(desc(TPM1)) |>
    dplyr::mutate("order_TPM1" = seq_len(nrow(perc_MYHs))) |>
    dplyr::arrange(desc(TPM3)) |>
    dplyr::mutate("order_TPM3" = seq_len(nrow(perc_MYHs))) |>
    tibble::rownames_to_column("fiber_ID")

# Curve fitting for TPM3 with MYH7

perc_TPM |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_7)) |>
    dplyr::select(MYH7, order_7, TPM3, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH7, TPM3),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_7,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("#440154FF",
                                           "darkblue")) +
    ggplot2::theme_minimal()

# Curve fitting for TPM1 with MYH2

perc_TPM |>
    dplyr::inner_join(perc_MYHs) |>
    dplyr::arrange(desc(order_2)) |>
    dplyr::select(MYH2, order_2, TPM1, fiber_ID) |>
    tidyr::pivot_longer(cols = c(MYH2, TPM1),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_2,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(
        alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c("#5DC863FF",
                                           "darkred")) +
    ggplot2::theme_minimal()

# Defining thresholds for each curve --------------------------------------


# knees TNNTs -------------------------------------------------------------

TNNT1_curve <- perc_troponins |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(TNNT1) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        TNNT1 = 0
    )) |>
    t()

TNNT1_curve <- DropletUtils::barcodeRanks(TNNT1_curve, lower = 10)

bottom_knee_TNNT1 <- 100 - TNNT1_curve@metadata$knee

perc_troponins |>
    ggplot2::ggplot(
        ggplot2::aes(order_TNNT1,
                     TNNT1
        )
    ) +
    ggplot2::geom_point(color = "#440154FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_TNNT1, color = "red") +
    ggplot2::theme_minimal()

# TNNT3

TNNT3_curve <- perc_troponins |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(TNNT3) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        TNNT3 = 0
    )) |>
    t()

TNNT3_curve <- DropletUtils::barcodeRanks(TNNT3_curve, lower = 10)

bottom_knee_TNNT3 <- 100 - TNNT3_curve@metadata$knee

perc_troponins |>
    ggplot2::ggplot(
        ggplot2::aes(order_TNNT3,
                     TNNT3
        )
    ) +
    ggplot2::geom_point(color = "#5DC863FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_TNNT3, color = "red") +
    ggplot2::theme_minimal()


# knees ATP2As ------------------------------------------------------------

# ATP2A2

ATP2A2_curve <- perc_ATP2A |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(ATP2A2) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        ATP2A2 = 0
    )) |>
    t()

ATP2A2_curve <- DropletUtils::barcodeRanks(ATP2A2_curve, lower = 10)

bottom_knee_ATP2A2 <- 100 - ATP2A2_curve@metadata$knee

perc_ATP2A |>
    ggplot2::ggplot(
        ggplot2::aes(order_ATP2A2,
                     ATP2A2
        )
    ) +
    ggplot2::geom_point(color = "#440154FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_ATP2A2, color = "red") +
    ggplot2::theme_minimal()

# ATP2A1

ATP2A1_curve <- perc_ATP2A |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(ATP2A1) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        ATP2A1 = 0
    )) |>
    t()

ATP2A1_curve <- DropletUtils::barcodeRanks(ATP2A1_curve, lower = 10)

bottom_knee_ATP2A1 <- 100 - ATP2A1_curve@metadata$knee

perc_ATP2A |>
    ggplot2::ggplot(
        ggplot2::aes(order_ATP2A1,
                     ATP2A1
        )
    ) +
    ggplot2::geom_point(color = "#5DC863FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_ATP2A1, color = "red") +
    ggplot2::theme_minimal()


# knees MYLs --------------------------------------------------------------

# MYL3

MYL3_curve <- perc_MYL |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYL3) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        MYL3 = 0
    )) |>
    t()

MYL3_curve <- DropletUtils::barcodeRanks(MYL3_curve, lower = 10)

bottom_knee_MYL3 <- 100 - MYL3_curve@metadata$knee

perc_MYL |>
    ggplot2::ggplot(
        ggplot2::aes(order_MYL3,
                     MYL3
        )
    ) +
    ggplot2::geom_point(color = "#440154FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYL3, color = "red") +
    ggplot2::theme_minimal()

# MYL1

MYL1_curve <- perc_MYL |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYL1) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        MYL1 = 0
    )) |>
    t()

MYL1_curve <- DropletUtils::barcodeRanks(MYL1_curve, lower = 10)

bottom_knee_MYL1 <- 100 - MYL1_curve@metadata$knee

perc_MYL |>
    ggplot2::ggplot(
        ggplot2::aes(order_MYL1,
                     MYL1
        )
    ) +
    ggplot2::geom_point(color = "#5DC863FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYL1, color = "red") +
    ggplot2::theme_minimal()

# knees TPMs --------------------------------------------------------------

# TPM3

TPM3_curve <- perc_TPM |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(TPM3) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        TPM3 = 0
    )) |>
    t()

TPM3_curve <- DropletUtils::barcodeRanks(TPM3_curve, lower = 10)

bottom_knee_TPM3 <- 100 - TPM3_curve@metadata$knee

perc_TPM |>
    ggplot2::ggplot(
        ggplot2::aes(order_TPM3,
                     TPM3
        )
    ) +
    ggplot2::geom_point(color = "#440154FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_TPM3, color = "red") +
    ggplot2::theme_minimal()

# TPM1

TPM1_curve <- perc_TPM |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(TPM1) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        TPM1 = 0
    )) |>
    t()

TPM1_curve <- DropletUtils::barcodeRanks(TPM1_curve, lower = 10)

bottom_knee_TPM1 <- 100 - TPM1_curve@metadata$knee

perc_TPM |>
    ggplot2::ggplot(
        ggplot2::aes(order_TPM1,
                     TPM1
        )
    ) +
    ggplot2::geom_point(color = "#5DC863FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_TPM1, color = "red") +
    ggplot2::theme_minimal()

# knees MYHs --------------------------------------------------------------

# MYH7

MYH7_curve <- perc_MYHs |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH7) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
       MYH7 = 0
    )) |>
    t()

MYH7_curve <- DropletUtils::barcodeRanks(MYH7_curve, lower = 10)

bottom_knee_MYH7 <- 100 - MYH7_curve@metadata$knee

perc_MYHs |>
    ggplot2::ggplot(
        ggplot2::aes(order_7,
                     MYH7
        )
    ) +
    ggplot2::geom_point(color = "#440154FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH7, color = "red") +
    ggplot2::theme_minimal()

# MYH2

MYH2_curve <- perc_MYHs |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH2) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    tidyr::replace_na(replace = list(
        MYH2 = 0
    )) |>
    t()

MYH2_curve <- DropletUtils::barcodeRanks(MYH2_curve, lower = 10)

bottom_knee_MYH2 <- 100 - MYH2_curve@metadata$knee

perc_MYHs |>
    ggplot2::ggplot(
        ggplot2::aes(order_2,
                     MYH2
        )
    ) +
    ggplot2::geom_point(color = "#5DC863FF") +
    # ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH2, color = "red") +
    ggplot2::theme_minimal()


# Fiber_typing ------------------------------------------------------------


# fiber typing TNNTs ------------------------------------------------------

fiber_type_TNNTs <- perc_troponins |>
    dplyr::mutate(
        fiber_type_TNNTs = dplyr::case_when(
            TNNT1 >= bottom_knee_TNNT1 &
                TNNT3 <= bottom_knee_TNNT3 ~ "Type 1",
            TNNT3 >= bottom_knee_TNNT3 &
                TNNT1 <= bottom_knee_TNNT1 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::select(c(fiber_type_TNNTs, fiber_ID))


# fiber typing ATP2A ------------------------------------------------------

fiber_type_ATP2A <- perc_ATP2A |>
    dplyr::mutate(
        fiber_type_ATP2A = dplyr::case_when(
            ATP2A2 >= bottom_knee_ATP2A2 &
                ATP2A1 <= bottom_knee_ATP2A1 ~ "Type 1",
            ATP2A1 >= bottom_knee_ATP2A1 &
                ATP2A2 <= bottom_knee_ATP2A2 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::select(c(fiber_type_ATP2A, fiber_ID))


# fiber typing MYL --------------------------------------------------------

fiber_type_MYL <- perc_MYL |>
    dplyr::mutate(
        fiber_type_MYL = dplyr::case_when(
            MYL3 >= bottom_knee_MYL3 &
                MYL1 <= bottom_knee_MYL1 ~ "Type 1",
            MYL1 >= bottom_knee_MYL1 &
                MYL3 <= bottom_knee_MYL3 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::select(c(fiber_type_MYL, fiber_ID))

# fiber typing TPM --------------------------------------------------------

fiber_type_TPM <- perc_TPM |>
    dplyr::mutate(
        fiber_type_TPM = dplyr::case_when(
            TPM3 >= bottom_knee_TPM3 &
                TPM1 <= bottom_knee_TPM1 ~ "Type 1",
            TPM1 >= bottom_knee_TPM1 &
                TPM3 <= bottom_knee_TPM3 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::select(c(fiber_type_TPM, fiber_ID))


# fiber typing MYHs -------------------------------------------------------

fiber_type_MYHs <- perc_MYHs |>
    dplyr::mutate(
        fiber_type_MYH = dplyr::case_when(
            MYH7 >= bottom_knee_MYH7 &
                MYH2 <= bottom_knee_MYH2 ~ "Type 1",
            MYH2 >= bottom_knee_MYH2 &
                MYH7 <= bottom_knee_MYH7 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::select(c(fiber_type_MYH, fiber_ID))

all_fiber_types <- fiber_type_MYHs |>
    dplyr::inner_join(fiber_type_ATP2A) |>
    dplyr::inner_join(fiber_type_MYL) |>
    dplyr::inner_join(fiber_type_TNNTs) |>
    dplyr::inner_join(fiber_type_TPM) |>
    tibble::column_to_rownames("fiber_ID")

slow_fibers <- all_fiber_types |>
    dplyr::filter(fiber_type_MYH == "Type 1")

upset_slow <- list(
    MYH = slow_fibers$fiber_type_MYH,
    TNNT = slow_fibers$fiber_type_TNNTs,
    ATP2A = slow_fibers$fiber_type_ATP2A,
    MYL = slow_fibers$fiber_type_MYL,
    TPM = slow_fibers$fiber_type_TPM
)

ComplexHeatmap::list_to_matrix(upset_slow)

# UpSetR::upset(all_fiber_types)


# Upset plot --------------------------------------------------------------

data_upset_slow <- all_fiber_types  |>
    dplyr::mutate(type_MYH_Slow = ifelse(fiber_type_MYH == "Type 1", 1, 0)) |>
    dplyr::mutate(type_TNNT_Slow = ifelse(fiber_type_TNNTs == "Type 1", 1, 0)) |>
    dplyr::mutate(type_ATP2A_Slow = ifelse(fiber_type_ATP2A == "Type 1", 1, 0)) |>
    dplyr::mutate(type_MYL_Slow = ifelse(fiber_type_MYL == "Type 1", 1, 0)) |>
    dplyr::mutate(type_TPM_Slow = ifelse(fiber_type_TPM == "Type 1", 1, 0))

data_upset_fast <- all_fiber_types  |>
    dplyr::mutate(type_MYH_Fast = ifelse(fiber_type_MYH == "Type 2A", 1, 0)) |>
    dplyr::mutate(type_TNNT_Fast = ifelse(fiber_type_TNNTs == "Type 2A", 1, 0)) |>
    dplyr::mutate(type_ATP2A_Fast = ifelse(fiber_type_ATP2A == "Type 2A", 1, 0)) |>
    dplyr::mutate(type_MYL_Fast = ifelse(fiber_type_MYL == "Type 2A", 1, 0))

data_upset_hybrid <- all_fiber_types  |>
    dplyr::mutate(type_MYH_Hybrid = ifelse(fiber_type_MYH == "Hybrid 1/2A", 1, 0)) |>
    dplyr::mutate(type_TNNT_Hybrid = ifelse(fiber_type_TNNTs == "Hybrid 1/2A", 1, 0)) |>
    dplyr::mutate(type_ATP2A_Hybrid = ifelse(fiber_type_ATP2A == "Hybrid 1/2A", 1, 0)) |>
    dplyr::mutate(type_MYL_Hybrid = ifelse(fiber_type_MYL == "Hybrid 1/2A", 1, 0))

# slow upset:

ComplexUpset::upset(
    data_upset_slow,
    c("type_MYH_Slow", "type_MYL_Slow", "type_ATP2A_Slow", "type_TNNT_Slow", "type_TPM_Slow"),
    intersections=list(
        'type_MYH_Slow',
        'type_MYL_Slow',
        'type_ATP2A_Slow',
        "type_TNNT_Slow",
        "type_TPM_Slow",
        c('type_MYH_Slow', 'type_MYL_Slow'),
        c('type_MYH_Slow', 'type_ATP2A_Slow'),
        c('type_MYH_Slow', 'type_TNNT_Slow'),
        c('type_MYL_Slow', 'type_ATP2A_Slow'),
        c('type_MYL_Slow', 'type_TNNT_Slow'),
        c('type_ATP2A_Slow', 'type_TNNT_Slow'),
        c('type_MYH_Slow', 'type_MYL_Slow', 'type_ATP2A_Slow'),
        c('type_MYH_Slow', 'type_MYL_Slow', 'type_TNNT_Slow'),
        c('type_MYH_Slow', 'type_ATP2A_Slow', 'type_TNNT_Slow'),
        c('type_MYL_Slow', 'type_ATP2A_Slow', 'type_TNNT_Slow'),

        c('type_MYH_Slow', "type_MYL_Slow", "type_ATP2A_Slow", "type_TNNT_Slow", "type_TPM_Slow"),
        c('type_MYH_Slow', "type_MYL_Slow", "type_ATP2A_Slow", "type_TNNT_Slow")
    ),
    queries=list(
        ComplexUpset::upset_query(set='type_MYH_Slow', fill="#134057"),
        ComplexUpset::upset_query(set='type_MYL_Slow', fill='#8CB3E8'),
        ComplexUpset::upset_query(set='type_ATP2A_Slow', fill='#BC4749'),
        ComplexUpset::upset_query(set='type_TNNT_Slow', fill='#417B5A')
    ),
    ggplot2::labs(x = "Gene combinations"),
    base_annotations=list(
        'Intersection size'=(
            ComplexUpset::intersection_size(
                bar_number_threshold=1,  # show all numbers on top of bars
                width=0.5,   # reduce width of the bars
            )
            # add some space on the top of the bars
            + ggplot2::scale_y_continuous(expand= ggplot2::expansion(mult=c(0, 0.05)), limits = c(0, 900))
            + ggplot2::ylab('Intersection of slow fibers')
            + ggplot2::ggtitle('Slow-classified fibers')
            + ggplot2::theme(
                # hide grid lines
                panel.grid = ggplot2::element_blank(),
                # show axis lines
                axis.line=ggplot2::element_line(colour='black'),
                text = ggplot2::element_text(face="bold", size=15, colour="black"),
                strip.text = ggplot2::element_text(colour = "white"),
                strip.background = ggplot2::element_rect(fill="black"),
                legend.position = "none",
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
        )
    ),
    stripes=ComplexUpset::upset_stripes(
        geom=ggplot2::geom_segment(linewidth=12),  # make the stripes larger
        colors=c('grey95', 'grey95')
    ),
    # to prevent connectors from getting the colorured, use `fill` instead of `color`, together with `shape='circle filled'`
    matrix=ComplexUpset::intersection_matrix(
        geom=ggplot2::geom_point(
            shape='circle filled',
            size=3.5,
            stroke=0.45
        )
    ),
    set_sizes=(
        ComplexUpset::upset_set_size(geom=ggplot2::geom_bar(width=0.4))
        + ggplot2::theme(
            axis.line.x= ggplot2::element_line(colour='black'),
            axis.ticks.x= ggplot2::element_line(),
            panel.grid= ggplot2::element_blank()
        )
    ),
    sort_sets=FALSE,
    sort_intersections='descending'
) +
    ggplot2::theme(panel.grid= ggplot2::element_blank())

# ggplot2::ggsave(here::here("doc/figures/figure_2/upset_slow_proteomics.png"))

# fast upset

ComplexUpset::upset(
    data_upset_fast,
    c("type_MYH_Fast", "type_MYL_Fast", "type_ATP2A_Fast", "type_TNNT_Fast"),
    intersections=list(
        'type_MYH_Fast',
        'type_MYL_Fast',
        'type_ATP2A_Fast',
        "type_TNNT_Fast",
        c('type_MYH_Fast', 'type_MYL_Fast'),
        c('type_MYH_Fast', 'type_ATP2A_Fast'),
        c('type_MYH_Fast', 'type_TNNT_Fast'),
        c('type_MYL_Fast', 'type_ATP2A_Fast'),
        c('type_MYL_Fast', 'type_TNNT_Fast'),
        c('type_ATP2A_Fast', 'type_TNNT_Fast'),
        c('type_MYH_Fast', 'type_MYL_Fast', 'type_ATP2A_Fast'),
        c('type_MYH_Fast', 'type_MYL_Fast', 'type_TNNT_Fast'),
        c('type_MYH_Fast', 'type_ATP2A_Fast', 'type_TNNT_Fast'),
        c('type_MYL_Fast', 'type_ATP2A_Fast', 'type_TNNT_Fast'),

        c('type_MYH_Fast', "type_MYL_Fast", "type_ATP2A_Fast", "type_TNNT_Fast")
    ),
    queries=list(
        ComplexUpset::upset_query(set='type_MYH_Fast', fill="#134057"),
        ComplexUpset::upset_query(set='type_MYL_Fast', fill='#8CB3E8'),
        ComplexUpset::upset_query(set='type_ATP2A_Fast', fill='#BC4749'),
        ComplexUpset::upset_query(set='type_TNNT_Fast', fill='#417B5A')
    ),
    ggplot2::labs(x = "Gene combinations"),
    base_annotations=list(
        'Intersection size'=(
            ComplexUpset::intersection_size(
                bar_number_threshold=1,  # show all numbers on top of bars
                width=0.5,   # reduce width of the bars
            )
            # add some space on the top of the bars
            + ggplot2::scale_y_continuous(expand= ggplot2::expansion(mult=c(0, 0.05)), limits = c(0, 900))
            + ggplot2::ylab('Intersection of Fast fibers')
            + ggplot2::ggtitle('Fast-classified fibers')
            + ggplot2::theme(
                # hide grid lines
                panel.grid = ggplot2::element_blank(),
                # show axis lines
                axis.line=ggplot2::element_line(colour='black'),
                text = ggplot2::element_text(face="bold", size=15, colour="black"),
                strip.text = ggplot2::element_text(colour = "white"),
                strip.background = ggplot2::element_rect(fill="black"),
                legend.position = "none",
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
        )
    ),
    stripes=ComplexUpset::upset_stripes(
        geom=ggplot2::geom_segment(linewidth=12),  # make the stripes larger
        colors=c('grey95', 'grey95')
    ),
    # to prevent connectors from getting the colorured, use `fill` instead of `color`, together with `shape='circle filled'`
    matrix=ComplexUpset::intersection_matrix(
        geom=ggplot2::geom_point(
            shape='circle filled',
            size=3.5,
            stroke=0.45
        )
    ),
    set_sizes=(
        ComplexUpset::upset_set_size(geom=ggplot2::geom_bar(width=0.4))
        + ggplot2::theme(
            axis.line.x= ggplot2::element_line(colour='black'),
            axis.ticks.x= ggplot2::element_line(),
            panel.grid= ggplot2::element_blank()
        )
    ),
    sort_sets=FALSE,
    sort_intersections='descending'
) +
    ggplot2::theme(panel.grid= ggplot2::element_blank())

# ggplot2::ggsave(here::here("doc/figures/figure_2/upset_fast_proteomics.png"))

# Upset hybrids

ComplexUpset::upset(
    data_upset_hybrid,
    c("type_MYH_Hybrid", "type_MYL_Hybrid", "type_ATP2A_Hybrid", "type_TNNT_Hybrid"),
    intersections=list(
        'type_MYH_Hybrid',
        'type_MYL_Hybrid',
        'type_ATP2A_Hybrid',
        "type_TNNT_Hybrid",
        c('type_MYH_Hybrid', 'type_MYL_Hybrid'),
        c('type_MYH_Hybrid', 'type_ATP2A_Hybrid'),
        c('type_MYH_Hybrid', 'type_TNNT_Hybrid'),
        c('type_MYL_Hybrid', 'type_ATP2A_Hybrid'),
        c('type_MYL_Hybrid', 'type_TNNT_Hybrid'),
        c('type_ATP2A_Hybrid', 'type_TNNT_Hybrid'),
        c('type_MYH_Hybrid', 'type_MYL_Hybrid', 'type_ATP2A_Hybrid'),
        c('type_MYH_Hybrid', 'type_MYL_Hybrid', 'type_TNNT_Hybrid'),
        c('type_MYH_Hybrid', 'type_ATP2A_Hybrid', 'type_TNNT_Hybrid'),
        c('type_MYL_Hybrid', 'type_ATP2A_Hybrid', 'type_TNNT_Hybrid'),

        c('type_MYH_Hybrid', "type_MYL_Hybrid", "type_ATP2A_Hybrid", "type_TNNT_Hybrid")
    ),
    queries=list(
        ComplexUpset::upset_query(set='type_MYH_Hybrid', fill="#134057"),
        ComplexUpset::upset_query(set='type_MYL_Hybrid', fill='#8CB3E8'),
        ComplexUpset::upset_query(set='type_ATP2A_Hybrid', fill='#BC4749'),
        ComplexUpset::upset_query(set='type_TNNT_Hybrid', fill='#417B5A')
    ),
    ggplot2::labs(x = "Gene combinations"),
    base_annotations=list(
        'Intersection size'=(
            ComplexUpset::intersection_size(
                bar_number_threshold=1,  # show all numbers on top of bars
                width=0.5,   # reduce width of the bars
            )
            # add some space on the top of the bars
            + ggplot2::scale_y_continuous(expand= ggplot2::expansion(mult=c(0, 0.05)), limits = c(0, 900))
            + ggplot2::ylab('Intersection of Hybrid fibers')
            + ggplot2::ggtitle('Hybrid-classified fibers')
            + ggplot2::theme(
                # hide grid lines
                panel.grid = ggplot2::element_blank(),
                # show axis lines
                axis.line=ggplot2::element_line(colour='black'),
                text = ggplot2::element_text(face="bold", size=15, colour="black"),
                strip.text = ggplot2::element_text(colour = "white"),
                strip.background = ggplot2::element_rect(fill="black"),
                legend.position = "none",
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
        )
    ),
    stripes=ComplexUpset::upset_stripes(
        geom=ggplot2::geom_segment(linewidth=12),  # make the stripes larger
        colors=c('grey95', 'grey95')
    ),
    # to prevent connectors from getting the colorured, use `fill` instead of `color`, together with `shape='circle filled'`
    matrix=ComplexUpset::intersection_matrix(
        geom=ggplot2::geom_point(
            shape='circle filled',
            size=3.5,
            stroke=0.45
        )
    ),
    set_sizes=(
        ComplexUpset::upset_set_size(geom=ggplot2::geom_bar(width=0.4))
        + ggplot2::theme(
            axis.line.x= ggplot2::element_line(colour='black'),
            axis.ticks.x= ggplot2::element_line(),
            panel.grid= ggplot2::element_blank()
        )
    ),
    sort_sets=FALSE,
    sort_intersections='descending'
) +
    ggplot2::theme(panel.grid= ggplot2::element_blank())

# ggplot2::ggsave(here::here("doc/figures/figure_2/upset_hybrid_proteomics.png"))

# All curves together -----------------------------------------------------

all_curves <- perc_MYHs |>
    dplyr::inner_join(perc_ATP2A) |>
    dplyr::inner_join(perc_MYL) |>
    dplyr::inner_join(perc_TPM) |>
    dplyr::inner_join(perc_troponins)

cor_ATP2A2 <- cor.test(all_curves$MYH7, all_curves$ATP2A2,  method = "spearman")
r_ATP2A2 <- round(cor_ATP2A2$estimate, digits = 2)
p_ATP2A2 <- cor_ATP2A2$p.value

cor_MYL3 <- cor.test(all_curves$MYH7, all_curves$MYL3,  method = "spearman")
r_MYL3 <- round(cor_MYL3$estimate, digits = 2)
p_MYL3 <- cor_MYL3$p.value

cor_TPM3 <- cor.test(all_curves$MYH7, all_curves$TPM3,  method = "spearman")
r_TPM3 <- round(cor_TPM3$estimate, digits = 2)
p_TPM3 <- cor_TPM3$p.value

cor_TNNT1 <- cor.test(all_curves$MYH7, all_curves$TNNT1,  method = "spearman")
r_TNNT1 <- round(cor_TNNT1$estimate, digits = 2)
p_TNNT1 <- cor_TNNT1$p.value

all_curves |>
    dplyr::select(MYH7,
                  TPM3,
                  TNNT1,
                  ATP2A2,
                  MYL3,
                  order_7) |>
    tidyr::pivot_longer(cols = c(MYH7,
                                 TPM3,
                                 TNNT1,
                                 ATP2A2,
                                 MYL3),
                        values_to = "perc_of_isoform_expression",
                        names_to = "Genes") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_7,
            y = perc_of_isoform_expression,
            color = Genes
        )
    ) +
    ggplot2::geom_point(alpha = 0.5, size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c(
        "#BC4749",
        "#134057",
        "orange",
        "#417B5A",
        "#8CB3E8"
    )) +
    ggplot2::annotate("text", x=600, y=80, label= expression(underline(bold("Correlation vs MYH7:"))), colour="#134057", fontface=2, hjust=0, parse=T, size=2.5) +
    ggplot2::annotate("text", x=600, y=40, label= paste("ATP2A2: r =", r_ATP2A2, ", p < 0.001"), colour="#BC4749", fontface=2, hjust=0, size=2.5) +
    ggplot2::annotate("text", x=600, y=50, label= paste("TNNT1: r =", r_TNNT1, ", p < 0.001"), colour="#417B5A", fontface=2, hjust=0, size=2.5) +
    ggplot2::annotate("text", x=600, y=60, label= paste("TPM3; r =", r_TPM3, ", p < 0.001"), colour="#8CB3E8", fontface=2, hjust=0, size=2.5) +
    ggplot2::annotate("text", x=600, y=30, label= paste("MYL3: r =", r_MYL3, ", p < 0.001"), colour="orange", fontface=2, hjust=0, size=2.5) +
    ggplot2::ggtitle("Slow fiber markers proteomics") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(face="bold", size=8, colour="black"),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill="black"),
        plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::ylab("% slow isoform expressed") +
    ggplot2::xlab("Sample (ranked by MYH7 expression)")

# ggplot2::ggsave(here::here("doc/figures/figure_2/slow_fiber_markers_curves_spearman.png"),
#                 units = "mm",
#                 height = 60,
#                 width = 128)

# Simple upset proteomics -------------------------------------------------

data_upset_slow <- all_fiber_types  |>
    dplyr::mutate(type_MYH_Slow = ifelse(fiber_type_MYH == "Type 1", 1, 0)) |>
    dplyr::mutate(type_TNNT_Slow = ifelse(fiber_type_TNNTs == "Type 1", 1, 0)) |>
    dplyr::mutate(type_ATP2A_Slow = ifelse(fiber_type_ATP2A == "Type 1", 1, 0)) |>
    dplyr::mutate(type_MYL_Slow = ifelse(fiber_type_MYL == "Type 1", 1, 0)) |>
    dplyr::mutate(type_TPM_Slow = ifelse(fiber_type_TPM == "Type 1", 1, 0))

data_upset_slow |>
    dplyr::select(dplyr::starts_with("type")) |>
    rowSums() |>
    table()

number_slow <- data_upset_slow |>
    dplyr::select(dplyr::starts_with("type")) |>
    rowSums() |>
    as.data.frame() |>
    dplyr::rename("sums" = 1) |>
    dplyr::filter(!sums == 0) |>
    nrow()

upset_slow_df <- data.frame(
    Five = 231/number_slow*100,
    Four = 80/number_slow*100,
    Three = 14/number_slow*100,
    Two = 7/number_slow*100,
    One = 13/number_slow*100
)  |>
    tidyr::pivot_longer(everything())

upset_slow_df$graph <- c(1:5)

upset_slow_df  |>
    ggplot2::ggplot(ggplot2::aes(
        x = graph,
        y = value
    )) +
    ggplot2::geom_bar(ggplot2::aes(fill = name),
             position = "dodge", stat = "summary", fun = "mean",
             width = 0.5, size = 0.5, alpha = 0.85, color = "black") +
    # ggplot2::geom_text(ggplot2::aes(label=value), position = ggplot2::position_dodge(width=0.9), vjust=-0.3, size=2) +
    ggplot2::scale_fill_manual(values = c("#08519c",
                                          "#3182bd",
                                          "#eff3ff",
                                          "#6baed6",
                                          "#bdd7e7")) +
    ggplot2::labs(
        x = "Protein combination",
        y = "Number of slow fibers (%)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    ggplot2::scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("Five", "Four", "Three", "Two", "One")) +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill="black"),
        legend.position = "none",
    )

# ggplot2::ggsave(here::here("doc/figures/figure_1/simple_upset_proteomics.png"),
#                 units = "mm",
#                 height = 60,
#                 width = 60)
