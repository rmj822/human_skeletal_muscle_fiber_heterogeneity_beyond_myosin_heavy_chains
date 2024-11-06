
# Loading and filtering raw data ------------------------------------------

proteomics_raw <- utils::read.table(
    here::here(
        {{"C:/Users/jns822/Desktop/Scripts/Heterofiber/data-raw/annotated_lib_based_heterofiber_all_samples.txt"}}
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

proteomics_filtered <- vroom::vroom(here::here("data/proteomics_working_data.csv")) |>
    dplyr::rename("gene" = "...1")

data_plot <- proteomics_controls |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::mutate(filtered = dplyr::case_when(
        Gene.name %in% proteomics_filtered$gene ~ "kept",
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
