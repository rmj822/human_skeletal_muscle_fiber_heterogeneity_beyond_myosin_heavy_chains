# Loading proteomics data -------------------------------------------------

proteomics_data <- vroom::vroom(
    here::here("data/data_proteomics_filtered.csv")
) |>
    dplyr::select(!1) |>
    as.data.frame() |>
    dplyr::filter(
        !Gene.name == "",
        !duplicated(Gene.name)
    ) |>
    tibble::column_to_rownames("Gene.name")

metadata <- vroom::vroom(
    here::here("data/metadata_proteomics.csv")
) |>
    dplyr::rename("fiberID" = 1)

heterofiber_number_quantified_proteins <- colSums(!is.na(proteomics_data)) |>
    as.data.frame()

colnames(heterofiber_number_quantified_proteins) <- "number_of_proteins"

heterofiber_number_quantified_proteins <- heterofiber_number_quantified_proteins |>
    dplyr::mutate(study = rep("heterofiber",
                              ncol(proteomics_data))) |>
    tibble::rownames_to_column("fiberID")

heterofiber_number_quantified_proteins <- metadata |>
    dplyr::filter(fiberID %in% heterofiber_number_quantified_proteins$fiberID) |>
    dplyr::filter(!duplicated(fiberID)) |>
    dplyr::select(fiberID,
                  subject) |>
    dplyr::inner_join(heterofiber_number_quantified_proteins) |>
    tibble::column_to_rownames("fiberID")

data_plot <- heterofiber_number_quantified_proteins |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::filter(!number_of_proteins < 1000) |>
    dplyr::mutate(project = dplyr::case_when(
        study == "heterofiber" ~ "proteomics",
        TRUE ~ "muscle_disease"
    ))

# Violins 1000 fiber proteome ---------------------------------------------
data_plot_heterofiber <- data_plot |>
    dplyr::filter(project == "proteomics")

data_plot_heterofiber$subject <- factor(data_plot_heterofiber$subject,
                                        levels = c("P1",
                                                   "P2",
                                                   "P3",
                                                   "P4",
                                                   "P5"))

violins_heterofiber <- data_plot_heterofiber |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = subject,
            y = number_of_proteins,
            fill = subject
        )
    ) +
    ggplot2::geom_violin(alpha = 0.75) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual(values = c(
        "#08519c",
        "#3182bd",
        "#6baed6",
        "#bdd7e7",
        "#eff3ff"
    ),
    guide = "none") +
    ggplot2::ylim(1000, 2500) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1000 fiber proteome",
                  y = "Number of quantified proteins") +
    ggplot2::theme(text = ggplot2::element_text(size = 9.5)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_discrete(labels = c("P1",
                                         "P2",
                                         "P3",
                                         "P4",
                                         "P5"))

# ggplot2::ggsave(here::here("doc/figures/figure_1/violins_heterofiber_participants.png"),
#                 height = 60,
#                 width = 60,
#                 units = "mm")

# CV ----------------------------------------------------------------------

proteomics_raw <- utils::read.table(
    here::here(
        "data-raw/loading_controls_proteomics.txt"
    ),
    header = TRUE,
    sep = "\t"
)

contaminants <- proteomics_raw |>
    dplyr::distinct(Gene.name, .keep_all = TRUE) |>
    dplyr::filter(!Gene.name == "") |>
    tibble::column_to_rownames("Gene.name") |>
    t() |>
    as.data.frame() |>
    dplyr::select(dplyr::contains("KRT"))

contaminants <- colnames(contaminants)

proteomics_raw <- proteomics_raw |>
    dplyr::filter(!Gene.name %in% contaminants)

Genes <- proteomics_raw |>
    dplyr::pull(Gene.name)


# Renaming columns and filtering for technical controls -------------------

old_column_names <- colnames(proteomics_raw)

new_column_names <- gsub(
    pattern = ".*DIA_",
    replacement = "",
    old_column_names
)

new_column_names <- gsub(
    pattern = "_.*",
    replacement = "",
    new_column_names
)

colnames(proteomics_raw) <- new_column_names

proteomics_controls <- proteomics_raw |>
    dplyr::select(dplyr::starts_with("C")) |>
    dplyr::bind_cols("Gene.name" = Genes) |>
    dplyr::select(!C11) |>
    dplyr::filter(!Gene.name == "",
                  !duplicated(Gene.name)) |>
    tibble::column_to_rownames("Gene.name")


# Calculating mean of intensities and CVs ---------------------------------


proteomics_controls <- proteomics_controls |>
    dplyr::mutate(coefficient_variation = apply(proteomics_controls,
                                                1,
                                                FUN = function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)) |>
    dplyr::mutate(mean_LFQ = rowMeans(proteomics_controls, na.rm = TRUE))


# Loading filtered dataset and annotating CV dataframe --------------------

proteomics_filtered <- vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1)

data_plot <- proteomics_controls |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::mutate(filtered = dplyr::case_when(
        Gene.name %in% proteomics_filtered$Gene.name ~ "kept",
        TRUE ~ "dropped"
    )) |>
    dplyr::mutate(perc_CVs = dplyr::case_when(
        coefficient_variation <= 30 & coefficient_variation >= 20 ~ "under_30",
        coefficient_variation <= 20 & coefficient_variation >= 10~ "under_20",
        coefficient_variation <= 10 ~ "under_10",
        TRUE ~ "high_perc"
    ))

table(data_plot$perc_CVs)


data_plot |>
    ggplot2::ggplot(
        ggplot2::aes(x = log2(mean_LFQ),
                     y = coefficient_variation,
                     color = perc_CVs)
    ) +
    ggplot2::geom_point(size = 0.5,
                        alpha = 0.5,
                        na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(
        values = c("#bdbdbd",
                   "#08519c",
                   "#3182bd",
                   "#6baed6")
    ) +
    ggplot2::geom_hline(
        yintercept = 10,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 24.5,
        y = 5,
        fontface = 2,
        label = "182",
        size = 2.25
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 23.5,
        y = 15,
        fontface = 2,
        label = "1153",
        size = 2.25
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 22.5,
        y = 25,
        fontface = 2,
        label = "2048",
        size = 2.25
    ) +
    ggplot2::geom_hline(
        yintercept = 20,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = 30,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::ylab("Mean coefficient of variation (%)") +
    ggplot2::xlab("LFQ intensities (Log2)") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size = 8),
        legend.position = "none",
    )
# ggplot2::coord_cartesian(xlim = c(7.55, 23), clip = "off")

ggplot2::ggsave(filename = here::here("doc/figures/figure_1/CVs_proteomics.png"),
                height = 60,
                width = 60,
                units = "mm")


# Dynamic range plot ------------------------------------------------------

sum_of_intensities <- colSums(proteomics_filtered,
                              na.rm = TRUE)

rel_abundance <- proteomics_filtered |>
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

mean_rel_abundance <- mean_rel_abundance |>
    dplyr::arrange(desc(log10_rel_abundance)) |>
    dplyr::mutate(
        order = 1:nrow(proteomics_filtered)
    )

Gene_names <- rownames(mean_rel_abundance)

mean_expression_all <- mean_rel_abundance |>
    tibble::rownames_to_column("gene")

ggplot2::ggplot() +

    # Add all genes
    ggplot2::geom_point(data = mean_expression_all,
                        aes(x=order,
                            y=log10_rel_abundance),
                        colour = "#045a8d",
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
                     mapping = aes(order, log10_rel_abundance, label = gene),
                     size = 1.8, label.padding=0.1, max.overlaps = Inf, min.segment.length=0.1, segment.size=0.2, force = 10) +

# Change design
    ylab("% total intensities, 10log") +
    xlab("Protein rank") +
    theme_classic() +
    ggtitle("Proteomics") +
    theme(
        text = element_text(face="bold", colour="black", size = 6),
        plot.title = element_text(face = "bold", color = "black", size = 8, hjust = 0.5),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        legend.position = "none",
    )

# ggsave(here::here("doc/figures/figure_1/dynamic_range_proteome.png"),
#        units = "mm",
#        height = 60,
#        width = 60)
