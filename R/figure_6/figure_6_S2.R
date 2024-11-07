################################################################################################################################################
##################################       Panel A      #####################################################
################################################################################################################################################

# Load data
actin_v_control <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls.csv")) |>
    dplyr::mutate(sign_actin = dplyr::case_when(
        signifficant == "Upregulated (π < 0.05)" ~ "UP",
        signifficant == "Upregulated (FDR < 0.05)" ~ "UP",
        signifficant == "Downregulated (π < 0.05)" ~ "DOWN",
        signifficant == "Downregulated (FDR < 0.05)" ~ "DOWN",
        TRUE ~ "Not significant"
    )) %>%
    dplyr::select(Genes, logFC, sign_actin) %>%
    dplyr::rename(logFC_actin = logFC)

troponin_v_control <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls.csv"))  |>
    dplyr::mutate(sign_troponin = dplyr::case_when(
        signifficant == "Upregulated (π < 0.05)" ~ "UP",
        signifficant == "Upregulated (FDR < 0.05)" ~ "UP",
        signifficant == "Downregulated (π < 0.05)" ~ "DOWN",
        signifficant == "Downregulated (FDR < 0.05)" ~ "DOWN",
        TRUE ~ "Not significant"
    )) |>
    dplyr::select(Genes, logFC, sign_troponin) %>%
    dplyr::rename(logFC_troponin = logFC)


# Combine both in same df
combined_data <- left_join(actin_v_control, troponin_v_control, by = "Genes")
combined_data$logFC_mean <- (combined_data$logFC_actin + combined_data$logFC_troponin) / 2

# Choose which proteins to highlight
combined_data <- combined_data %>%
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "CNP",
            "MFAP4",
            "H3-3B",
            "NRDC",
            "UGDH",
            "MYH4",
            "SKIC8",
            "TNNT3",
            "TNNI2",
            "TNNC2",
            "THBS1",
            "MARS2",
            "TPM3",
            "MYL2",
            "TNNI1",
            "MYH7",
            "MYH2",
            "MYH1",
            "TNNT1"
        ) ~ Genes,
        TRUE ~ ""
    ))

# Make column for colouring
combined_data <- combined_data %>%
    dplyr::mutate("colour_plot" = dplyr::case_when(
        Genes %in% c("CNP") ~ "Both",
        Genes %in% c("MFAP4",
                     "H3-3B",
                     "NRDC") ~ "TNNT1 up",
        Genes %in% c("UGDH",
                     "SKIC8") ~ "TNNT1 down",
        Genes %in% c("THBS1") ~ "ACTA1 down",
        Genes %in% c("MARS2",
                     "TNNT1") ~ "ACTA1 up",
        Genes %in% c("MYH4",
                     "TNNT3",
                     "TNNI2",
                     "TNNC2",
                     "TPM3",
                     "MYL2",
                     "TNNI1",
                     "MYH7",
                     "MYH2",
                     "MYH1") ~ "Sarcomeric",
        TRUE ~ "Standard"
    ))

# Perform correlation analysis
cor <- cor.test(combined_data$logFC_actin, combined_data$logFC_troponin,  method = "spearman") # r = 0.736
round(cor$estimate, 3)
formatC(cor$p.value, format = "e", digits = 2)



# Create plot
logFC_comparison <- ggplot(combined_data,
                           aes(x = logFC_actin,
                               y = logFC_troponin,
                               label = names)) +
    geom_point(aes(colour = colour_plot), size = 0.5) +
    scale_color_manual(values = c("Standard" = "grey",
                                  "Both" = "#5DC863FF",
                                  "TNNT1 up" = "#d662c4",
                                  "TNNT1 down" = "#d662c4",
                                  "ACTA1 up" = "#E48C2AFF",
                                  "ACTA1 down" = "#E48C2AFF",
                                  "Sarcomeric" = "steelblue"
    )) +
    labs(
        title = "ACTA1-NM vs TNNT-NM",
        x = "log2FC (ACTA1-NM - control)",
        y = "log2FC (TNNT1-NM - control)",
    ) +
    theme_classic() +
    theme(
        text = element_text(face="bold", size=7, colour="black"),
        axis.text = element_text(colour="black", face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 8)
    ) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_label_repel(data = combined_data %>% dplyr::filter(!names == ""),
                     mapping = aes(x = logFC_actin,
                                   y = logFC_troponin,
                                   fill = colour_plot,
                                   label = names),
                     size = 1.8, colour = "black", label.padding = 0.1, min.segment.length = 0.1, segment.size = 0.2, max.overlaps = Inf) +
    scale_fill_manual(values = c("Standard" = "grey",
                                 "Both" = "#5DC863FF",
                                 "TNNT1 up" = "#d662c4",
                                 "TNNT1 down" = "#d662c4",
                                 "ACTA1 up" = "#E48C2AFF",
                                 "ACTA1 down" = "#E48C2AFF",
                                 "Sarcomeric" = "steelblue"
    )) +
    annotate("text", x=4, y=-2, label= "r = 0.736", colour="black", fontface=2, size=2) +
    annotate("text", x=4, y=-2.5, label= "p < 0.001", colour="black", fontface=2, size=2)




ggplot2::ggsave(logFC_comparison,
                file = "doc/figures/figure_6_S2/figure_6_S2A.pdf",
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
##################################       Panel B      #####################################################
################################################################################################################################################

data_pseudobulk <- vroom::vroom(here::here("data/data_MD_pseudobulk.csv")) |>
    tibble::column_to_rownames("Genes")

data_grouping <- data.frame(
    "sample_ID" = colnames(data_pseudobulk),
    "condition" = c(
        "actin",
        "actin",
        "actin",
        "control",
        "control",
        "control",
        "troponin",
        "troponin",
        "troponin"
    )
) |>
    dplyr::mutate(subject = gsub(
        pattern = ".*_",
        replacement = "",
        colnames(data_pseudobulk)
    )) |>
    dplyr::mutate(
        subject = dplyr::case_when(
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
        )
    )

data_pca <- data_pseudobulk |>
    limma::normalizeQuantiles() |>
    na.omit()

pca_object <- prcomp(t(data_pca), scale. = TRUE)

factoextra::fviz_eig(pca_object)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    dplyr::mutate(PC1 = PC1 * -1) |>
    # dplyr::mutate(PC2 = PC2 * -1) |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(data_grouping)

data_pca |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = condition,
        label = subject
    )) +
    ggplot2::geom_point(
        size = 1,
        alpha = 0.75,
    ) +
    ggplot2::theme_minimal() +
    # ggplot2::ggtitle("PCA Myopathies") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        ),
        legend.position = "none",
    ) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           vjust = 1.5)
    ) +
    ggplot2::theme(
        # legend.position = "none",
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        axis.text = ggplot2::element_text(size=8),
    ) +
    ggplot2::xlab("PC1 (35.4%)") +
    ggplot2::ylab("PC2 (22.5%)") +
    ggplot2::scale_color_manual(
        "condition:",
        values = c("#E48C2AFF",
                   "#969594",
                   "#d662c4")
    ) +
    ggrepel::geom_label_repel(data = data_pca,
                              mapping = ggplot2::aes(x = PC1,
                                                     y = PC2,
                                                     fill = condition),
                              color = "black",
                              size = 2,
                              label.padding = 0.1,
                              min.segment.length = 0.1,
                              segment.size = 0.2,
                              force = 10
    ) +
    ggplot2::scale_fill_manual(
        values = c("#f0bd87",
                   "#d5d4d4",
                   "#eab0e1")
    ) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))


ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2B.png"),
                units = "mm",
                height = 55,
                width = 55)

################################################################################################################################################
##################################       Panel C      #####################################################
################################################################################################################################################

# First Comparison actin vs control ---------------------------------------

up_actin_v_controls <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls.csv")) |>
    dplyr::mutate(up_actin_v_controls = dplyr::case_when(
        signifficant == "Upregulated (π < 0.05)" ~ 1,
        signifficant == "Upregulated (FDR < 0.05)" ~ 1,
        TRUE ~ 0
    )) |>
    dplyr::select(up_actin_v_controls, Genes)

down_actin_v_controls <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls.csv")) |>
    dplyr::mutate(down_actin_v_controls = dplyr::case_when(
        signifficant == "Downregulated (π < 0.05)" ~ 1,
        signifficant == "Downregulated (FDR < 0.05)" ~ 1,
        TRUE ~ 0
    )) |>
    dplyr::select(down_actin_v_controls, Genes)


# Second comparison Troponin vs control -----------------------------------

up_troponin_v_controls <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls.csv")) |>
    dplyr::mutate(up_troponin_v_controls = dplyr::case_when(
        signifficant == "Upregulated (π < 0.05)" ~ 1,
        signifficant == "Upregulated (FDR < 0.05)" ~ 1,
        TRUE ~ 0
    )) |>
    dplyr::select(up_troponin_v_controls, Genes)

down_troponin_v_controls <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls.csv")) |>
    dplyr::mutate(down_troponin_v_controls = dplyr::case_when(
        signifficant == "Downregulated (π < 0.05)" ~ 1,
        signifficant == "Downregulated (FDR < 0.05)" ~ 1,
        TRUE ~ 0
    )) |>
    dplyr::select(down_troponin_v_controls, Genes)


# Third comparison actin vs troponin --------------------------------------

down_troponin_v_actin <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin.csv")) |>
    dplyr::mutate(down_troponin_v_actin = dplyr::case_when(
        signifficant == "Upregulated (π < 0.05)" ~ 1,
        signifficant == "Upregulated (FDR < 0.05)" ~ 1,
        TRUE ~ 0
    )) |>
    dplyr::select(down_troponin_v_actin, Genes)

up_troponin_v_actin <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin.csv")) |>
    dplyr::mutate(up_troponin_v_actin = dplyr::case_when(
        signifficant == "Downregulated (π < 0.05)" ~ 1,
        signifficant == "Downregulated (FDR < 0.05)" ~ 1,
        TRUE ~ 0
    )) |>
    dplyr::select(up_troponin_v_actin, Genes)

# Upset plot for all comparisons

data_upset_all <- up_actin_v_controls |>
    dplyr::inner_join(down_actin_v_controls) |>
    dplyr::inner_join(up_troponin_v_controls) |>
    dplyr::inner_join(down_troponin_v_controls) |>
    dplyr::inner_join(up_troponin_v_actin) |>
    dplyr::inner_join(down_troponin_v_actin)

comparison_upset_all <- c("up_actin_v_controls",
                          "down_actin_v_controls",
                          "up_troponin_v_controls",
                          "down_troponin_v_controls",
                          "up_troponin_v_actin",
                          "down_troponin_v_actin")

upset_myopathies_vs_control_all <- ComplexUpset::upset(data_upset_all,
                                                       comparison_upset_all,
                                                       name = "Myopathy combinations",
                                                       width_ratio=0.1,
                                                       wrap = TRUE,
                                                       min_size = 1, # Min size of group to be included
                                                       min_degree = 1, # Min intersection size
                                                       keep_empty_groups = TRUE,
                                                       themes =  ComplexUpset::upset_default_themes(text=element_text(size = 6)),
                                                       labeller=ggplot2::as_labeller(c(
                                                           'up_actin_v_controls' = 'Up ACTA1-NM vs Control',
                                                           'down_actin_v_controls' = 'Down ACTA1-NM vs Control',
                                                           'up_troponin_v_controls' = 'Up TNNT1-NM vs Control',
                                                           'down_troponin_v_controls' = 'Down TNNT1-NM vs Control',
                                                           'up_troponin_v_actin' = 'Up TNNT1-NM vs ACTA1-NM',
                                                           'down_troponin_v_actin' = 'Up ACTA1-NM vs TNNT1-NM'
                                                       )),
                                                       base_annotations = list( # Anything related to bars
                                                           'Intersection size' = ComplexUpset::intersection_size(
                                                               bar_number_threshold=1,  # show all numbers on top of bars
                                                               width=0.5, # width of bars,
                                                               mapping=aes(fill='bars_color'),
                                                               text = list(
                                                                   size = 2
                                                               )
                                                           ) +
                                                               scale_fill_manual(values = c("bars_color" = "grey20")) +
                                                               ylab('Intersection size of DE proteins') +
                                                               ggtitle("DE proteins across myopathies") +
                                                               ylim(c(0,210)) +
                                                               theme(
                                                                   panel.grid = ggplot2::element_blank(), # hide grid lines
                                                                   axis.line=ggplot2::element_line(colour='black'),  # show axis lines
                                                                   text = ggplot2::element_text(face="bold", size=6, colour="black"),
                                                                   legend.position = "none",
                                                                   plot.title = ggplot2::element_text(hjust = 0.5, size = 8)
                                                               )

                                                       ),
                                                       set_sizes= FALSE, # Anything related to set sizes
                                                       # to prevent connectors from getting the colorured, use `fill` instead of `color`, together with `shape='circle filled'`
                                                       matrix=(
                                                           ComplexUpset::intersection_matrix(geom=geom_point(shape='circle filled', size=1.5))
                                                           + scale_color_manual(
                                                               values=c('up_actin_v_controls'='#E48C2AFF', 'down_actin_v_controls'='darkgrey', 'up_troponin_v_controls'='#d662c4',
                                                                        'down_troponin_v_controls'='lightgrey', 'up_troponin_v_actin'='red', 'down_troponin_v_actin'='blue'),
                                                               guide=guide_legend(override.aes=list(shape='circle'))
                                                           )
                                                       ),
                                                       queries=list(
                                                           ComplexUpset::upset_query(set='up_actin_v_controls', fill='#E48C2AFF'),
                                                           ComplexUpset::upset_query(set='down_actin_v_controls', fill='darkgrey'),
                                                           ComplexUpset::upset_query(set='up_troponin_v_controls', fill='#d662c4'),
                                                           ComplexUpset::upset_query(set='down_troponin_v_controls', fill='lightgrey'),
                                                           ComplexUpset::upset_query(set='up_troponin_v_actin', fill='red'),
                                                           ComplexUpset::upset_query(set='down_troponin_v_actin', fill='blue')
                                                       )
)

ggplot2::ggsave(upset_myopathies_vs_control_all,
                file = "doc/figures/figure_6_S2/figure_6_S2C.png",
                units = "mm",
                height = 60,
                width = 128)

################################################################################################################################################
##################################       Panels D, E, F, G, H and I      #####################################################
################################################################################################################################################

# fiber type-specific pseudobulk and limma analysis

pseudobulk_maker_fiber_type <- function(.data, metadata, subject_id, grouping, grouping_2, colname) {
    selection_vector <- metadata |>
        dplyr::filter(condition == grouping) |>
        dplyr::filter(fiber_type == grouping_2) |>
        dplyr::filter(subject == subject_id) |>
        dplyr::pull(fiber_ID)

    pseudobulk_median_calculation <- .data |>
        as.data.frame() |>
        dplyr::select(selection_vector) |>
        t() |>
        as.data.frame() |>
        dplyr::mutate(dplyr::across(.cols = everything(), median, na.rm = T)) |>
        unique() |>
        t()

    colnames(pseudobulk_median_calculation) <- colname
    return(pseudobulk_median_calculation)
}

data_pseudobulk_ft <- data.frame(
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 1",
        subject_id = "C1",
        colname = "C1_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 1",
        subject_id = "C2",
        colname = "C2_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 1",
        subject_id = "C3",
        colname = "C3_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 2A",
        subject_id = "C1",
        colname = "C1_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 2A",
        subject_id = "C2",
        colname = "C2_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 2A",
        subject_id = "C3",
        colname = "C3_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "A1",
        colname = "A1_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "A2",
        colname = "A2_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "A3",
        colname = "A3_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "A1",
        colname = "A1_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "A2",
        colname = "A2_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "A3",
        colname = "A3_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "T1",
        colname = "T1_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "T2",
        colname = "T2_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "T3",
        colname = "T3_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "T1",
        colname = "T1_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "T2",
        colname = "T2_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "T3",
        colname = "T3_T2A"
    )
)

data_limma <- data_pseudobulk_ft

data_grouping_ft <- data.frame(
    "fiber_ID" = colnames(data_limma),
    "condition" = c(
        "control",
        "control",
        "control",
        "control",
        "control",
        "control",
        "actin",
        "actin",
        "actin",
        "actin",
        "actin",
        "actin",
        "troponin",
        "troponin",
        "troponin",
        "troponin",
        "troponin",
        "troponin"
    ),
    "fiber_type" = c(
        "type_1",
        "type_1",
        "type_1",
        "type_2A",
        "type_2A",
        "type_2A",
        "type_1",
        "type_1",
        "type_1",
        "type_2A",
        "type_2A",
        "type_2A",
        "type_1",
        "type_1",
        "type_1",
        "type_2A",
        "type_2A",
        "type_2A"
    )
)

# Type 1 limma -----------------------------------------------

data_grouping_type_1 <- data_grouping_ft |>
    dplyr::filter(fiber_type == "type_1")

data_limma_type_1 <- data_limma |>
    dplyr::select(data_grouping_type_1$fiber_ID) |>
    limma::normalizeBetweenArrays()

vector_factor <- factor(data_grouping_type_1$condition,
                        levels = c(
                            "actin",
                            "control",
                            "troponin"
                        )
)

design_matrix <-
    model.matrix(
        ~ 0 + vector_factor,
        data_grouping_type_1
    )

colnames(design_matrix) <- c(
    "actin",
    "control",
    "troponin"
)

fit <- limma::lmFit(
    data_limma_type_1,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "actin_vs_control" = actin - control,
    "troponin_vs_control" = troponin - control,
    "actin_vs_troponin" = troponin - actin,
    levels = design_matrix
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp_slow <- limma::eBayes(tmp)

DE_actin_controls_type_1 <- limma::topTable(tmp_slow,
                                            coef = "actin_vs_control",
                                            sort.by = "P",
                                            n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_troponin_controls_type_1 <- limma::topTable(tmp_slow,
                                               coef = "troponin_vs_control",
                                               sort.by = "P",
                                               n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_actin_troponin_type_1 <- limma::topTable(tmp_slow,
                                            coef = "actin_vs_troponin",
                                            sort.by = "P",
                                            n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


# ACTA1-NM - control ------------------------------------------------------

DE_actin_controls_type_1 <- DE_actin_controls_type_1 |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "TPPP3",
            "CFL",
            "ANXA1",
            "S100A11",
            "RRAD",
            "ANXA5",
            "MRRF",
            "RPS11",
            "GFPT1",
            "NDUFA11",
            "RPS17",
            "RPL35",
            "NDUFS3"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_actin_controls_type_1 |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            color = signifficant,
            names = Genes
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(DE_actin_controls_type_1$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_actin_controls_type_1$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_actin_controls_type_1$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#969594",
        "lightgrey",
        "#E48C2AFF"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("ACTA1-NM vs control type 1") +
    ggplot2::xlab("log2FC (ACTA1-NM - control in type 1)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 5,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_actin_controls_type_1 |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 10
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#d5d4d4",
        "#f0bd87"
    ))

ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2D.pdf"),
                units = "mm",
                height = 50,
                width = 45)

# TNNT1-NM - control ------------------------------------------------------

DE_troponin_controls_type_1 <- DE_troponin_controls_type_1 |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "RPL23",
            "TPPP3",
            "ANXA5",
            "TNNC2",
            "FABP4",
            "ANXA2",
            "FASN",
            "COPZ1",
            "VPS37C",
            "NDUFS3",
            "CAV1",
            "UGDH",
            "RPL35",
            "RPS17"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_troponin_controls_type_1 |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            color = signifficant,
            names = Genes
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(DE_troponin_controls_type_1$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_troponin_controls_type_1$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_troponin_controls_type_1$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#969594",
        "lightgrey",
        "#d662c4"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs control type 1") +
    ggplot2::xlab("log2FC (TNNT1-NM - control in type 1)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 5,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_troponin_controls_type_1 |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 25
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#d5d4d4",
        "#eab0e1"
    ))

ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2F.pdf"),
                units = "mm",
                height = 50,
                width = 45)

# TNNT1-NM - ACTA1-NM ------------------------------------------------------

DE_actin_troponin_type_1 <- DE_actin_troponin_type_1 |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "HBG2",
            "MBP",
            "ATP5MF",
            "RPL23",
            "NDUFA11",
            "SDHC",
            "COPZ1",
            "PSME3",
            "RCN3",
            "TNNT1",
            "UGDH",
            "TPPP3"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_actin_troponin_type_1 |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            color = signifficant,
            names = Genes
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(DE_actin_troponin_type_1$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_actin_troponin_type_1$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_actin_troponin_type_1$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#E48C2AFF",
        "lightgrey",
        "#d662c4"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs ACTA1-NM type 1") +
    ggplot2::xlab("log2FC (TNNT1-NM - ACTA1-NM in type 1)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 5,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_actin_troponin_type_1 |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 10
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#f0bd87",
        "#eab0e1"
    ))

ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2H.pdf"),
                units = "mm",
                height = 50,
                width = 45)

# Type 2 limma -----------------------------------------------

data_grouping_type_2 <- data_grouping_ft |>
    dplyr::filter(fiber_type == "type_2A")

data_limma_type_2 <- data_limma |>
    dplyr::select(data_grouping_type_2$fiber_ID) |>
    limma::normalizeBetweenArrays()

vector_factor <- factor(data_grouping_type_2$condition,
                        levels = c(
                            "actin",
                            "control",
                            "troponin"
                        )
)

design_matrix <-
    model.matrix(
        ~ 0 + vector_factor,
        data_grouping_type_2
    )

colnames(design_matrix) <- c(
    "actin",
    "control",
    "troponin"
)

fit <- limma::lmFit(
    data_limma_type_2,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "actin_vs_control" = actin - control,
    "troponin_vs_control" = troponin - control,
    "actin_vs_troponin" = troponin - actin,
    levels = design_matrix
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp_fast <- limma::eBayes(tmp)

DE_actin_controls_type_2 <- limma::topTable(tmp_fast,
                                            coef = "actin_vs_control",
                                            sort.by = "P",
                                            n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_troponin_controls_type_2 <- limma::topTable(tmp_fast,
                                               coef = "troponin_vs_control",
                                               sort.by = "P",
                                               n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_actin_troponin_type_2 <- limma::topTable(tmp_fast,
                                            coef = "actin_vs_troponin",
                                            sort.by = "P",
                                            n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


# ACTA1-NM - control ------------------------------------------------------

DE_actin_controls_type_2 <- DE_actin_controls_type_2 |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "TPPP3",
            "EI3C",
            "BGN",
            "PI16",
            "ANXA1",
            "SOD3",
            "S100A11",
            "RPL37A",
            "ATP5PB",
            "AIFM1",
            "TXNRD2",
            "CKMT2",
            "MAIP1",
            "LACTB"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_actin_controls_type_2 |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            color = signifficant,
            names = Genes
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(DE_actin_controls_type_2$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_actin_controls_type_2$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_actin_controls_type_2$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#969594",
        "lightgrey",
        "#E48C2AFF",
        "#E48C2AFF"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("ACTA1-NM vs control type 2A") +
    ggplot2::xlab("log2FC (ACTA1-NM  - control in type 2A)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 5,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_actin_controls_type_2 |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 25
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#d5d4d4",
        "#f0bd87",
        "#f0bd87"
    ))

ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2E.pdf"),
                units = "mm",
                height = 50,
                width = 45)

# TNNT1-NM - control ------------------------------------------------------

DE_troponin_controls_type_2 <- DE_troponin_controls_type_2 |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "PPIB",
            "TPPP3",
            "DPYSL3",
            "OGN",
            "EPPK1",
            "ANXA5",
            "PPIA",
            "S100A11",
            "ATP5B",
            "SLC25A12",
            "NDUFA12",
            "MYL4",
            "MT-CO2",
            "VDAC1",
            "CMT2"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_troponin_controls_type_2 |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            color = signifficant,
            names = Genes
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(DE_troponin_controls_type_2$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_troponin_controls_type_2$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_troponin_controls_type_2$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#969594",
        "#969594",
        "lightgrey",
        "#d662c4",
        "#d662c4"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs control type 2A") +
    ggplot2::xlab("log2FC (TNNT1-NM - control in type 2A)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 5,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_troponin_controls_type_2 |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 50
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#d5d4d4",
        "#eab0e1"
    ))

ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2G.pdf"),
                units = "mm",
                height = 50,
                width = 45)

# TNNT1-NM - ACTA1-NM ------------------------------------------------------

DE_actin_troponin_type_2 <- DE_actin_troponin_type_2 |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "ATP2A3",
            "PPIB",
            "APOC3",
            "TMED9",
            "CALML3",
            "GMPPB",
            "NDUFA11",
            "ACTC1",
            "FHL3",
            "CRYM",
            "RPL8",
            "IGKV2D-24"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_actin_troponin_type_2 |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            color = signifficant,
            names = Genes
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(DE_actin_troponin_type_1$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_actin_troponin_type_1$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_actin_troponin_type_1$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#E48C2AFF",
        "lightgrey",
        "#d662c4"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs ACTA1-NM type 2A") +
    ggplot2::xlab("log2FC (TNNT1-NM - ACTA1-NM in type 2A)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 5,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_actin_troponin_type_2 |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 200
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#f0bd87",
        "#eab0e1"
    ))

ggplot2::ggsave(here::here("doc/figures/figure_6_S2/figure_6_S2I.pdf"),
                units = "mm",
                height = 50,
                width = 45)

readr::write_csv(DE_actin_controls_type_1,
                 file = here::here("data/DE_analysis_MD/actin_v_controls_type_1.csv"))

readr::write_csv(DE_troponin_controls_type_1,
                 file = here::here("data/DE_analysis_MD/troponin_v_controls_type_1.csv"))

readr::write_csv(DE_actin_troponin_type_1,
                 file = here::here("data/DE_analysis_MD/actin_v_troponin_type_1.csv"))

readr::write_csv(DE_actin_controls_type_2,
                 file = here::here("data/DE_analysis_MD/actin_v_controls_type_2.csv"))

readr::write_csv(DE_troponin_controls_type_2,
                 file = here::here("data/DE_analysis_MD/troponin_v_controls_type_2.csv"))

readr::write_csv(DE_actin_troponin_type_2,
                 file = here::here("data/DE_analysis_MD/actin_v_troponin_type_2.csv"))
