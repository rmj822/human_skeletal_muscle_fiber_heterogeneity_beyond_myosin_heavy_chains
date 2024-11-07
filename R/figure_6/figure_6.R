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

# Save DE results:

write.csv(DE_actin_controls,
          here::here("data/DE_analysis_MD/actin_v_controls.csv"))

write.csv(DE_troponin_controls,
          here::here("data/DE_analysis_MD/troponin_v_controls.csv"))

write.csv(DE_troponin_actin,
          here::here("data/DE_analysis_MD/actin_v_troponin.csv"))


################################################################################################################################################
#################################################     Panel H  ##############################################################
################################################################################################################################################

library(org.Hs.eg.db)

# actin_vs_controls --------------------------------------------------------

actin_v_controls <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_actin_control <- clusterProfiler::gseGO(
    actin_v_controls,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_actin_control <- clusterProfiler::simplify(GSEA_actin_control,
                                                cutoff = 0.5,
                                                by = "p.adjust",
                                                select_fun = min)

dotplot_title <- c(
    `activated` = "Actin",
    `suppressed` = "Control"
)

clusterProfiler::dotplot(GSEA_actin_control,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("GSEA Actin vs Control") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

# troponin vs control --------------------------------------------------------

troponin_v_controls <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_troponin_control <- clusterProfiler::gseGO(
    troponin_v_controls,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_troponin_control <- clusterProfiler::simplify(GSEA_troponin_control,
                                                   cutoff = 0.5,
                                                   by = "p.adjust",
                                                   select_fun = min)

dotplot_title <- c(
    `activated` = "Troponin",
    `suppressed` = "Control"
)

clusterProfiler::dotplot(GSEA_troponin_control,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("Gene Set Enrichment Analysis") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


# troponin_vs_actin -------------------------------------------------------

troponin_v_actin <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_troponin_v_actin <- clusterProfiler::gseGO(
    troponin_v_actin,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

# GSEA_troponin_v_actin <- clusterProfiler::simplify(GSEA_troponin_v_actin,
#                                                    cutoff = 0.5,
#                                                    by = "p.adjust",
#                                                    select_fun = min)

dotplot_title <- c(
    `activated` = "Actin",
    `suppressed` = "Troponin"
)

clusterProfiler::dotplot(GSEA_troponin_v_actin,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("Gene Set Enrichment Analysis") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


# actin_vs_controls in type 1 --------------------------------------------------------

actin_v_controls_type_1 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls_type_1.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_actin_control_type_1 <- clusterProfiler::gseGO(
    actin_v_controls_type_1,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

# GSEA_actin_control_type_1 <- clusterProfiler::simplify(GSEA_actin_control_type_1,
#                                                 cutoff = 0.5,
#                                                 by = "p.adjust",
#                                                 select_fun = min)

dotplot_title <- c(
    `activated` = "Actin",
    `suppressed` = "Control"
)

clusterProfiler::dotplot(GSEA_actin_control_type_1,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("GSEA Actin vs Control Type 1 fibers") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

# troponin_vs_controls in type 1 --------------------------------------------------------

troponin_v_controls_type_1 <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls_type_1.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_troponin_control_type_1 <- clusterProfiler::gseGO(
    troponin_v_controls_type_1,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_troponin_control_type_1 <- clusterProfiler::simplify(GSEA_troponin_control_type_1,
                                                          cutoff = 0.5,
                                                          by = "p.adjust",
                                                          select_fun = min)

dotplot_title <- c(
    `activated` = "Troponin",
    `suppressed` = "Control"
)

clusterProfiler::dotplot(GSEA_troponin_control_type_1,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("GSEA Troponin vs Control Type 1 fibers") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

# troponin_vs_actin_type_1 -------------------------------------------------------

# It doesn't find any signifficantly enriched term
troponin_v_actin_type_1 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin_type_1.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_troponin_v_actin_type_1 <- clusterProfiler::gseGO(
    troponin_v_actin_type_1,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_troponin_v_actin_type_1 <- clusterProfiler::simplify(GSEA_troponin_v_actin_type_1,
                                                          cutoff = 0.5,
                                                          by = "p.adjust",
                                                          select_fun = min)
#
# dotplot_title <- c(
#     `activated` = "Actin",
#     `suppressed` = "Troponin"
# )
#
# clusterProfiler::dotplot(GSEA_troponin_v_actin_type_1,
#                          showCategory = 10,
#                          split = ".sign",
#                          font.size = 8) +
#     ggplot2::facet_grid(. ~ .sign,
#                         labeller = ggplot2::as_labeller(dotplot_title)) +
#     ggplot2::ggtitle("GSEA Troponin vs actin Type 1 fibers") +
#     ggplot2::scale_color_viridis_c(option = "plasma") +
#     ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
#     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

# actin_vs_controls in type 2 --------------------------------------------------------

actin_v_controls_type_2 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls_type_2.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_actin_control_type_2 <- clusterProfiler::gseGO(
    actin_v_controls_type_2,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

# GSEA_actin_control_type_2 <- clusterProfiler::simplify(GSEA_actin_control_type_2,
#                                                        cutoff = 0.5,
#                                                        by = "p.adjust",
#                                                        select_fun = min)

dotplot_title <- c(
    `activated` = "Actin",
    `suppressed` = "Control"
)

clusterProfiler::dotplot(GSEA_actin_control_type_2,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("GSEA Actin vs Control Type 2A fibers") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

# troponin_vs_controls in type 2 --------------------------------------------------------

troponin_v_controls_type_2 <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls_type_2.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_troponin_control_type_2 <- clusterProfiler::gseGO(
    troponin_v_controls_type_2,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_troponin_control_type_2 <- clusterProfiler::simplify(GSEA_troponin_control_type_2,
                                                          cutoff = 0.5,
                                                          by = "p.adjust",
                                                          select_fun = min)

dotplot_title <- c(
    `activated` = "Troponin",
    `suppressed` = "Control"
)

clusterProfiler::dotplot(GSEA_troponin_control_type_2,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("GSEA Troponin vs Control Type 2A fibers") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


# troponin_vs_actin_type_2 -------------------------------------------------------

troponin_v_actin_type_2 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin_type_2.csv")) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::pull(logFC, name = Genes)

GSEA_troponin_v_actin_type_2 <- clusterProfiler::gseGO(
    troponin_v_actin_type_2,
    keyType = "SYMBOL",
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_troponin_v_actin_type_2 <- clusterProfiler::simplify(GSEA_troponin_v_actin_type_2,
                                                          cutoff = 0.5,
                                                          by = "p.adjust",
                                                          select_fun = min)

dotplot_title <- c(
    `activated` = "Actin",
    `suppressed` = "Troponin"
)

clusterProfiler::dotplot(GSEA_troponin_v_actin_type_2,
                         showCategory = 10,
                         split = ".sign",
                         font.size = 8) +
    ggplot2::facet_grid(. ~ .sign,
                        labeller = ggplot2::as_labeller(dotplot_title)) +
    ggplot2::ggtitle("GSEA Troponin vs actin Type 2 fibers") +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

# Creating Combined GSEA plot ---------------------------------------------

result_GSEA_actin_control <- GSEA_actin_control@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_actin_control,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_control.csv"))

result_GSEA_actin_control_type_1 <- GSEA_actin_control_type_1@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_actin_control_type_1,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_control_type_1.csv"))

result_GSEA_actin_control_type_2 <- GSEA_actin_control_type_2@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_actin_control_type_2,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_control_type_2.csv"))

result_GSEA_troponin_control <- GSEA_troponin_control@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_troponin_control,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_troponin_control.csv"))

result_GSEA_troponin_control_type_1 <- GSEA_troponin_control_type_1@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_troponin_control_type_1,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_troponin_control_type_1.csv"))

result_GSEA_troponin_control_type_2 <- GSEA_troponin_control_type_2@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_troponin_control_type_2,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_troponin_control_type_2.csv"))

result_GSEA_actin_troponin <- GSEA_troponin_v_actin@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_actin_troponin,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_troponin.csv"))

result_GSEA_actin_troponin_type_1 <- GSEA_troponin_v_actin_type_1@result |>
    tibble::remove_rownames()
#
write.csv(result_GSEA_actin_troponin_type_1,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_troponin_type_1.csv"))

result_GSEA_actin_troponin_type_2 <- GSEA_troponin_v_actin_type_2@result |>
    tibble::remove_rownames()

write.csv(result_GSEA_actin_troponin_type_2,
          here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_troponin_type_2.csv"))

# Figure GSEA -------------------------------------------------------------

selected_terms <- c(
    "oxidative phosphorylation",
    "immune system process",
    "cell adhesion"
)


data_plot_actin_control <- result_GSEA_actin_control |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_actin",
        TRUE ~ "downregulated_actin"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

data_plot_actin_control_type_1 <- result_GSEA_actin_control_type_1 |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_actin_type_1",
        TRUE ~ "downregulated_actin_type_1"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

data_plot_actin_control_type_2 <- result_GSEA_actin_control_type_2 |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_actin_type_2",
        TRUE ~ "downregulated_actin_type_2"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

data_plot_troponin_control <- result_GSEA_troponin_control |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_troponin",
        TRUE ~ "downregulated_troponin"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

data_plot_troponin_control_type_1 <- result_GSEA_troponin_control_type_1 |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_troponin_type_1",
        TRUE ~ "downregulated_troponin_type_1"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

#creating empty column for upregulated_troponin_type_1 in dot plot
data_plot_troponin_control_type_1[2,1] <- "oxidative phosphorylation"
data_plot_troponin_control_type_1[2,2] <- "upregulated_troponin_type_1"
data_plot_troponin_control_type_1[2,3] <- NA
data_plot_troponin_control_type_1[2,4] <- NA

data_plot_troponin_control_type_2 <- result_GSEA_troponin_control_type_2 |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_troponin_type_2",
        TRUE ~ "downregulated_troponin_type_2"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

data_plot_actin_troponin <- result_GSEA_actin_troponin |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_troponin_vs_actin",
        TRUE ~ "upregulated_actin_vs_troponin"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

data_plot_actin_troponin_type_1 <- result_GSEA_actin_troponin_type_1 |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_troponin_vs_actin_type_1",
        TRUE ~ "upregulated_actin_vs_troponin_type_1"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

#creating empty columns for upregulated_troponin_vs_actin_type_1 and upregulated_actin_vs_troponin_type_1 in dot plot
data_plot_actin_troponin_type_1[1,1] <- "oxidative phosphorylation"
data_plot_actin_troponin_type_1[1,2] <- "upregulated_troponin_vs_actin_type_1"
data_plot_actin_troponin_type_1[1,3] <- NA
data_plot_actin_troponin_type_1[1,4] <- NA

data_plot_actin_troponin_type_1[2,1] <- "oxidative phosphorylation"
data_plot_actin_troponin_type_1[2,2] <- "upregulated_actin_vs_troponin_type_1"
data_plot_actin_troponin_type_1[2,3] <- NA
data_plot_actin_troponin_type_1[2,4] <- NA


data_plot_actin_troponin_type_2 <- result_GSEA_actin_troponin_type_2 |>
    dplyr::filter(Description %in% selected_terms) |>
    dplyr::mutate(condition = dplyr::case_when(
        enrichmentScore > 0 ~ "upregulated_troponin_vs_actin_type_2",
        TRUE ~ "upregulated_actin_vs_troponin_type_2"
    )) |>
    dplyr::select(Description, condition, p.adjust, enrichmentScore)

#creating empty columns for upregulated_troponin_vs_actin_type_2 in dot plot
data_plot_actin_troponin_type_2[2,1] <- "oxidative phosphorylation"
data_plot_actin_troponin_type_2[2,2] <- "upregulated_troponin_vs_actin_type_2"
data_plot_actin_troponin_type_2[2,3] <- NA
data_plot_actin_troponin_type_2[2,4] <- NA



data_plot_combined <- rbind(data_plot_actin_control,
                            data_plot_actin_control_type_1,
                            data_plot_actin_control_type_2,
                            data_plot_troponin_control,
                            data_plot_troponin_control_type_1,
                            data_plot_troponin_control_type_2,
                            data_plot_actin_troponin,
                            data_plot_actin_troponin_type_1,
                            data_plot_actin_troponin_type_2) |>
    dplyr::mutate(enrichmentScore = abs(enrichmentScore))


desired_order <- c(
    rep("downregulated_actin", 3),
    rep("downregulated_actin_type_1", 3),
    rep("downregulated_actin_type_2", 2),
    rep("upregulated_actin", 1),
    rep("upregulated_actin_type_1", 3),
    rep("upregulated_actin_type_2", 1),
    rep("downregulated_troponin", 3),
    rep("downregulated_troponin_type_1", 3),
    rep("downregulated_troponin_type_2", 3),
    rep("upregulated_troponin", 3),
    rep("upregulated_troponin_type_1", 1),
    rep("upregulated_troponin_type_2", 3),
    rep("upregulated_troponin_vs_actin", 2),
    rep("upregulated_troponin_vs_actin_type_1", 1),
    rep("upregulated_troponin_vs_actin_type_2", 1),
    rep("upregulated_actin_vs_troponin", 3),
    rep("upregulated_actin_vs_troponin_type_1", 1),
    rep("upregulated_actin_vs_troponin_type_2", 1)
)

data_plot_combined <- data_plot_combined |>
    dplyr::mutate(order = factor(condition,
                                 levels = unique(desired_order))) |>
    dplyr::mutate(order_terms = factor(Description,
                                       levels = c(
                                           "oxidative phosphorylation",
                                           "immune system process",
                                           "cell adhesion"
                                       ))) |>
    dplyr::arrange(order_terms) |>
    dplyr::arrange(order)


# forcats::fct_reorder(data_plot_combined$condition, desired_order)

data_plot_combined |>
    ggplot2::ggplot(
        ggplot2::aes(x = order,
                     y = order_terms,
                     size = enrichmentScore,
                     color = -log10(p.adjust))
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c("-log10\n(P.value)", option = "plasma") +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        text = ggplot2::element_text(size = 8),
        axis.title = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(4, "mm"),
        legend.position = "bottom",
        plot.margin=grid::unit(c(1,0,0,0), "mm")
    )

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6H.pdf"),
                height = 60,
                width = 128,
                units = "mm"
)

################################################################################################################################################
#################################################     Panel I and J  ##############################################################
################################################################################################################################################

# Coloring PCA by GO terms ------------------------------------------------

library(org.Hs.eg.db)
meta_go2gene <- as.list(org.Hs.egGO2EG)
meta_go2goName <- as.list(GO.db::GOTERM)

# naming the goID with a human readable title. AnnotationDBI handles these 'higher' level classes.
meta_go2goName <- sapply(
    meta_go2goName,
    AnnotationDbi::Term
)

# using a named vector to not mixup our names
names(meta_go2gene) <- meta_go2goName[names(meta_go2gene)]

# translating entrez ids to hugo symbols.
meta_entrez2symbol <- AnnotationDbi::select(
    org.Hs.eg.db,
    unique(unlist(meta_go2gene)),
    "SYMBOL",
    "ENTREZID"
)

# using the power of data.table to quickly match vectors in a list
meta_entrez2symbol <- data.table::as.data.table(meta_entrez2symbol)
data.table::setkey(meta_entrez2symbol, ENTREZID)

# data.table syntax is unique, but it is also blazing fast.
meta_go2gene <- lapply(meta_go2gene, \(i){
    meta_entrez2symbol[i, unique(SYMBOL)]
})
meta_go2gene <- reshape2::melt(meta_go2gene)
colnames(meta_go2gene) <- c("SYMBOL", "GO")

# Mitochondrial term ------------------------------------------------------

vector_mito <- meta_go2gene[grep("mitochondria",
                                 meta_go2gene$GO),
                            "SYMBOL"] |>
    unique()

mito_filtered <- data_wrangled |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_mito) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(data_wrangled, na.rm = T)

rel_abundance_mito <- (colSums(mito_filtered,
                               na.rm = T) / total_intensity) * 100

data_pca_mito <- as.data.frame(rel_abundance_mito) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca)


data_mito_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2) |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ .x * 50)) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% vector_mito)

data_mito_proteins <- data_mito_proteins |>
    dplyr::mutate("x.centroid" = 0) |>
    dplyr::mutate("y.centroid" = 0) |>
    dplyr::mutate(label = Genes)

data_pca_mito |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = rel_abundance_mito
    )) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_viridis_c("Rel. \nAbundance", option = "plasma") +
    # ggplot2::geom_segment(
    #     data = data_mito_proteins,
    #     ggplot2::aes(x = x.centroid,
    #                  y = y.centroid,
    #                  xend = Dim.1,
    #                  yend = Dim.2),
    #     color = "black",
    #     alpha = 0.5,
    #     arrow = grid::arrow(type = "closed",
    #                         length = ggplot2::unit(1, "mm"))
    # ) +
    # ggrepel::geom_label_repel(data = data_mito_proteins,
    #                           ggplot2::aes(
    #                               x = Dim.1,
    #                               y = Dim.2,
    #                               label = label),
    #                           size = 2,
    #                           color = "black",
    #                           label.padding = 0.1,
    #                           min.segment.length = 0.1,
    #                           segment.size = 0.2,
    #                           force = 20,
    #                           max.overlaps = Inf) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA by mitochondria GO term") +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=7),
        axis.text = ggplot2::element_text(size=7),
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           vjust = 3),
        legend.position = "bottom",
        legend.key.height = ggplot2::unit(2, "mm")
    )

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6I.png"),
                height = 60,
                width = 60,
                units = "mm")

# Coloring by ECM ---------------------------------------------------------

vector_ecm <- meta_go2gene[grep("extracellular matrix",
                                meta_go2gene$GO),
                           "SYMBOL"] |>
    unique()

ecm_filtered <- data_wrangled |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_ecm) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(data_wrangled, na.rm = T)

rel_abundance_ecm <- (colSums(ecm_filtered,
                              na.rm = T) / total_intensity) * 100

data_pca_ecm <- as.data.frame(rel_abundance_ecm) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca)


data_ecm_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2) |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ .x * 50)) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% vector_ecm)

data_ecm_proteins <- data_ecm_proteins |>
    dplyr::mutate("x.centroid" = 0) |>
    dplyr::mutate("y.centroid" = 0) |>
    dplyr::mutate(label = Genes)

data_pca_ecm |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = rel_abundance_ecm
    )) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_viridis_c("Rel. \nAbundance", option = "plasma") +
    # ggplot2::geom_segment(
    #     data = data_ecm_proteins,
    #     ggplot2::aes(x = x.centroid,
    #                  y = y.centroid,
    #                  xend = Dim.1,
    #                  yend = Dim.2),
    #     alpha = 0.5,
    #     color = "black",
    #     arrow = grid::arrow(type = "closed",
    #                         length = ggplot2::unit(1, "mm"))
    # ) +
    # ggrepel::geom_label_repel(data = data_ecm_proteins,
    #                           ggplot2::aes(
    #                               x = Dim.1,
    #                               y = Dim.2,
    #                               label = label),
    #                           size = 2,
    #                           color = "black",
    #                           label.padding = 0.1,
    #                           min.segment.length = 0.1,
    #                           segment.size = 0.2,
    #                           force = 20,
    #                           max.overlaps = Inf) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA by ECM GO term") +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=7),
        axis.text = ggplot2::element_text(size=7),
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           vjust = 3),
        legend.position = "bottom",
        legend.key.height = ggplot2::unit(2, "mm")
    )

ggplot2::ggsave(here::here("doc/figures/figure_6/figure_6J.png"),
                height = 60,
                width = 60,
                units = "mm")
