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

################################################################################################################################################
#################################################     Panel C  ##############################################################
################################################################################################################################################

muscle_disease_data <- vroom::vroom(
    here::here("data/data_muscle_disease.csv")
) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene_name") |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        log2
    )) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::arrange(desc(fiber_ID)) |>
    tibble::column_to_rownames("fiber_ID") |>
    t() |>
    as.data.frame()

number_of_quantified_proteins <- colSums(!is.na(muscle_disease_data)) |>
    as.data.frame() |>
    dplyr::rename(
        "number_quantified_proteins" = 1
    ) |>
    dplyr::filter(number_quantified_proteins > 1000) |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::pull("fiber_ID")

data_col_filtered <- muscle_disease_data |>
    dplyr::select(number_of_quantified_proteins)

metadata_MD_filtered <- metadata_MD |>
    dplyr::filter(fiber_ID %in% colnames(data_col_filtered)) |>
    dplyr::mutate(
condition = dplyr::case_when(
    condition == "TNNT1_nemaline_myopaty" ~ "TNNT1-NM",
    condition == "control" ~ "Control",
    TRUE ~ "ACTA1-NM"
)
    )

data_row_filtered <- data_col_filtered |>
    PhosR::selectGrps(
        grps = metadata_MD_filtered$condition,
        percent = 0.7
    )

data_wrangled <- data_row_filtered |>
    limma::normalizeQuantiles() |>
    PhosR::tImpute()

pca_object <- prcomp(t(data_wrangled), scale = TRUE)

factoextra::fviz_eig(pca_object)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2, PC3) |>
    tibble::add_column(
        fiber_type = metadata_MD_filtered$fiber_type,
        condition = metadata_MD_filtered$condition,
        fiberID = metadata_MD_filtered$fiber_ID,
        subject = as.factor(metadata_MD_filtered$subject)
    )

pca_object_MD <- pca_object

data_pca$fiber_type <- factor(data_pca$fiber_type,
                              levels = c("Hybrid 1/2A",
                                         "Hybrid 2A/2X",
                                         "Type 2A",
                                         "Type 1"))

# PCA disease:

data_pca |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = condition,
            names = fiberID
        )
    ) +
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
        legend.position = "bottom",
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
    ggplot2::xlab("PC1 (23.2%)") +
    ggplot2::ylab("PC2 (13.2%)") +
    ggplot2::scale_color_manual(
        "condition:",
        values = c("#E48C2AFF",
                   "#969594",
                   "#d662c4")
    ) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6C.png"),
                units = "mm",
                height = 55,
                width = 55)

################################################################################################################################################
#################################################     Panel D  ##############################################################
################################################################################################################################################

# Loading heterofiber data ------------------------------------------------

original_data <-vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    log2() |>
    as.data.frame()

metadata <- vroom::vroom(
    here::here("data/metadata_proteomics_fiber_type.csv")
) |>
    dplyr::select(!1)

metadata_MD <- metadata_MD |>
    dplyr::arrange(desc(fiber_ID))

disease_data <- data_row_filtered

# Arrange to same order and merge:

common_genes <- intersect(rownames(original_data), rownames(disease_data))

original_data <- original_data |>
    dplyr::filter(row.names(original_data) %in% common_genes) |>
    tibble::rownames_to_column("Gene.name")

disease_data <- disease_data |>
    dplyr::filter(row.names(disease_data) %in% common_genes) |>
    tibble::rownames_to_column("Gene.name")

both_datasets <- original_data |>
    dplyr::inner_join(disease_data) |>
    tibble::column_to_rownames("Gene.name")

# Quantile normalize both datasets together:

norm_both_datasets <- both_datasets |>
    limma::normalizeBetweenArrays(method = "quantile") |>
    as.data.frame() |>
    tibble::rownames_to_column("Gene.name")

# Running PCA only on 1000 fibers:

norm_original_data <- norm_both_datasets |>
    dplyr::select(colnames(original_data)) |>
    tibble::column_to_rownames("Gene.name") |>
    limma::removeBatchEffect(
        batch = metadata$digestion_batch,
        batch2 = metadata$MS_batch
    ) |>
    as.data.frame() |>
    PhosR::tImpute()

original_pca_object <- prcomp(t(norm_original_data),
                              scale. = TRUE)

pca_plot_data <- original_pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    dplyr::mutate(PC2 = PC2 * -1) |>
    tibble::add_column(
        fiber_type = metadata$fiber_type,
        fiberID = rownames(metadata),
        digestion_batch = as.factor(metadata$digestion_batch),
        subject = as.factor(metadata$subject),
        date = metadata$date_isolation,
        length = metadata$fiber_length
    )

# Calculating disease data projection:

norm_disease_data <- norm_both_datasets |>
    dplyr::select(colnames(disease_data)) |>
    tibble::column_to_rownames("Gene.name") |>
    PhosR::tImpute()

disease_projection <- scale(t(norm_disease_data),
                            center = original_pca_object$center,
                            scale = original_pca_object$scale) %*% original_pca_object$rotation |>
    as.data.frame() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample") |>
    dplyr::select(sample,
                  PC1,
                  PC2) |>
    dplyr::mutate(PC2 = PC2 * -1) |>
    dplyr::inner_join(metadata_MD_filtered |>
                          dplyr::rename(sample = fiber_ID) |>
                          dplyr::select(sample, condition, fiber_type)) |>
    dplyr::rename(fiberID = sample)

disease_projection$fiber_type <- factor(disease_projection$fiber_type,
                                        levels = c("Type 1",
                                                   "Hybrid 1/2A",
                                                   "Type 2A",
                                                   "Hybrid 2A/2X"))

# Plot PCA:

pca_plot_data$fiber_type <- factor(pca_plot_data$fiber_type,
                                   levels = c("Type 1",
                                              "Hybrid 1/2A",
                                              "Type 2A",
                                              "Hybrid 2A/2X"))

pca_plot_data |>

    ggplot2::ggplot(ggplot2::aes(x = PC1,
                                 y = PC2,
                                 color = fiber_type)) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Projection over 1000 fiber PCA plot") +
    ggplot2::stat_ellipse(
        geom = "polygon",
        ggplot2::aes(fill = fiber_type,
                     color = fiber_type),
        alpha = 0.15
    ) +
    ggplot2::scale_fill_manual(
        "Fiber type area\n1000 fiber dataset",
        values = c("#440154FF",
                   "#8CB3E8",
                   "#5DC863FF",
                   "#fdc325")
    ) +
    ggplot2::scale_color_manual(
        "outline",
        values = c(
            ggplot2::alpha("#8CB3E8", 0),
            ggplot2::alpha("#fdc325", 0),
            ggplot2::alpha("#5DC863FF", 0),
            ggplot2::alpha("#440154FF", 0)
        ),
        guide = "none"
    ) +
    ggnewscale::new_scale_colour() +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = PC1,
            y = PC2,
            color = fiber_type,
            fill = condition
        ),
        size = 1,
        stroke = 0.7,
        shape = 21,
        data = disease_projection
    ) +
    ggplot2::scale_color_manual(
        "Myopathy dataset \nsample fiber type",
        values = c("#440154FF",
                   "#8CB3E8",
                   "#5DC863FF",
                   "#fdc325")
    ) +
    ggplot2::scale_fill_manual("Myopathy dataset \nsample condition",
                               values = c("#E48C2AFF",
                                          "#969594",
                                          "#d662c4")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8,
                                                      face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      vjust = 1.5)) +
    ggplot2::theme(
        legend.position = "right",
        text = ggplot2::element_text(
            face = "bold",
            colour = "black",
            size = 8
        ),
        axis.text = ggplot2::element_text(size = 8),
    ) +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)")

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6D.png"),
                height = 60,
                width = 90,
                units = "mm")

################################################################################################################################################
#################################################     Panel E, F and G  ##############################################################
################################################################################################################################################

# Pseudobulk in limma -----------------------------------------------------

pseudobulk_maker <- function(.data, metadata, subject_id, grouping, colname) {
    selection_vector <- metadata |>
        dplyr::filter(condition == grouping) |>
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

data_pseudobulk <- data.frame(
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "A1",
        grouping = "ACTA1_nemaline_myopaty",
        colname = "actin_A1"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "A2",
        grouping = "ACTA1_nemaline_myopaty",
        colname = "actin_A2"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "A3",
        grouping = "ACTA1_nemaline_myopaty",
        colname = "actin_A3"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "C1",
        grouping = "control",
        colname = "control_C1"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "C2",
        grouping = "control",
        colname = "control_C2"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "C3",
        grouping = "control",
        colname = "control_C3"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "T1",
        grouping = "TNNT1_nemaline_myopaty",
        colname = "troponin_T1"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "T2",
        grouping = "TNNT1_nemaline_myopaty",
        colname = "troponin_T2"
    ),
    pseudobulk_maker(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        subject_id = "T3",
        grouping = "TNNT1_nemaline_myopaty",
        colname = "troponin_T3"
    )
)

export <- data_pseudobulk |>
    tibble::rownames_to_column("Genes")

# readr::write_csv(export,
#                  here::here("data/data_MD_pseudobulk.csv"))

# Running limma on pseudobulk data ----------------------------------------

data_limma <- data_pseudobulk |>
    limma::normalizeBetweenArrays(method = "quantile")

data_grouping <- data.frame(
    "fiber_ID" = colnames(data_pseudobulk),
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
)

vector_factor <- factor(data_grouping$condition,
                        levels = c(
                            "actin",
                            "control",
                            "troponin"
                        )
)

design_matrix <-
    model.matrix(
        ~ 0 + vector_factor,
        data_grouping
    )

colnames(design_matrix) <- c(
    "actin",
    "control",
    "troponin"
)

fit <- limma::lmFit(
    data_limma,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "actin_vs_control" = actin - control,
    "troponin_vs_control" = troponin - control,
    "troponin_vs_actin" = troponin - actin,
    levels = design_matrix
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp <- limma::eBayes(tmp)

DE_results <- limma::topTable(tmp,
                              sort.by = "F",
                              n = Inf
)

# ACTA1 - control ---------------------------------------------------------

DE_actin_controls <- limma::topTable(tmp,
                                     coef = "actin_vs_control",
                                     sort.by = "P",
                                     n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_actin_controls <- DE_actin_controls |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "S100A11",
            "ANXA2",
            "TPPP3",
            "SOD3",
            "ANXA1",
            "TNNT1",
            "ANXA5",
            "THBS1",
            "DLAT",
            "AIFM1",
            "NDUFS3",
            "RPS17",
            "CNP"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_actin_controls |>
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
        shape = ifelse(DE_actin_controls$signifficant == "not signifficant", 16, 1),
        alpha = ifelse(DE_actin_controls$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#969594",
        "lightgrey",
        "#E48C2AFF"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("ACTA1-NM vs control ") +
    ggplot2::xlab("log2FC (ACTA1-NM - control)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            size = 5,
            colour = "black",
            face = "bold"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_actin_controls |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.5,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 1
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#d5d4d4",
        "#f0bd87"
    )) +
    ggplot2::xlim(c(-4, 5.5)) +
    ggplot2::ylim(c(0, 5))

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6E.png"),
                units = "mm",
                height = 45,
                width = 45)


# TNNT1 - control ---------------------------------------------------------

DE_troponin_controls <- DE_troponin_controls |>
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
            "SOD3",
            "ANXA5",
            "TPPP3",
            "TNNT1",
            "S100A11",
            "ANXA2",
            "ANXA1",
            "PFKL",
            "CKMT2",
            "SKIC8",
            "TIMM8A",
            "NDUFAF5",
            "COQ5"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_troponin_controls |>
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
        shape = ifelse(DE_troponin_controls$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_troponin_controls$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_troponin_controls$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#969594",
        "#969594",
        "lightgrey",
        "#d662c4",
        "#d662c4"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs control") +
    ggplot2::xlab("log2FC (TNNT1-NM - control)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            size = 5,
            colour = "black",
            face = "bold"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_troponin_controls |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.5,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 40
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#d5d4d4",
        "#d5d4d4",
        "#eab0e1",
        "#eab0e1"
    )) +
    ggplot2::xlim(c(-4, 5.5)) +
    ggplot2::ylim(c(0, 5))

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6F.png"),
                units = "mm",
                height = 45,
                width = 45)

DE_troponin_controls <- limma::topTable(tmp,
                                        coef = "troponin_vs_control",
                                        sort.by = "P",
                                        n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_troponin_actin <- limma::topTable(tmp,
                                     coef = "troponin_vs_actin",
                                     sort.by = "P",
                                     n = Inf
) |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_troponin_actin <- DE_troponin_actin |>
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
            "TPM3",
            "EIF2B2",
            "TMED9",
            "TNNT1",
            "OGN",
            "UGDH",
            "UBE2G1",
            "HMGCS2",
            "SKIC8",
            "COL6A6",
            "MFAP4",
            "THBS1"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_troponin_actin |>
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
        shape = ifelse(DE_troponin_actin$signifficant == "Upregulated (π < 0.05)", 1,
                       ifelse(DE_troponin_actin$signifficant == "Downregulated (π < 0.05)", 1, 16)),
        alpha = ifelse(DE_troponin_actin$signifficant == "not signifficant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c(
        "#E48C2AFF",
        "lightgrey",
        "#d662c4"
    )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs ACTA1-NM") +
    ggplot2::xlab("log2FC (TNNT1-NM - ACTA1-NM") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            size = 5,
            colour = "black",
            face = "bold"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_troponin_actin |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 1.5,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 10
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#f0bd87",
        "#eab0e1"
    ))  +
    ggplot2::xlim(c(-4, 5.5)) +
    ggplot2::ylim(c(0, 5))

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6G.png"),
                units = "mm",
                height = 45,
                width = 45)
