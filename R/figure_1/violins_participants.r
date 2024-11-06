
# Loading proteomics data -------------------------------------------------

proteomics_data <- vroom::vroom(
    here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_proteomics_filtered.csv"),
    col_select = !c(1)
) |>
    as.data.frame() |>
    dplyr::filter(
        !Gene.name == "",
        !duplicated(Gene.name)
    ) |>
    tibble::column_to_rownames("Gene.name")

metadata <- vroom::vroom(
    "C:/Users/jns822/Desktop/Scripts/Heterofiber/data-raw/metadata.txt"
)


# Loading muscle disease data ---------------------------------------------

muscle_disease_data <- vroom::vroom(
    here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_muscle_disease.csv")
) |>
    tibble::column_to_rownames("...1")

metadata_muscle_disease <- vroom::vroom(here::here("data/metadata_MD_w_fiber_type_w_anonim.csv")) |>
    dplyr::mutate(anonimized_subject = dplyr::case_when(
        anonimized_subject == "T1" ~ "TM1",
        anonimized_subject == "T2" ~ "TM2",
        anonimized_subject == "T3" ~ "TM3",
        anonimized_subject == "A1" ~ "AM1",
        anonimized_subject == "A2" ~ "AM2",
        anonimized_subject == "A3" ~ "AM3",
        anonimized_subject == "C1" ~ "C1",
        anonimized_subject == "C2" ~ "C2",
        anonimized_subject == "C3" ~ "C3",
        TRUE ~ "error"
    ))

# Calculating number quantified proteins/study ----------------------------

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

muscle_disease_number_quantified_proteins <- colSums(!is.na(muscle_disease_data)) |>
    as.data.frame() |>
    tibble::rownames_to_column("SampleID")

colnames(muscle_disease_number_quantified_proteins) <- c(
    "fiber_ID",
    "number_of_proteins"
)

muscle_disease_number_quantified_proteins <- metadata_muscle_disease |>
    dplyr::select(fiber_ID, condition, anonimized_subject) |>
    dplyr::inner_join(muscle_disease_number_quantified_proteins) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::rename("study" = "condition")

# Merging and plotting ----------------------------------------------------

data_plot <- heterofiber_number_quantified_proteins |>
    dplyr::bind_rows(
        muscle_disease_number_quantified_proteins
    ) |>
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
                                        levels = c("FOR2",
                                                   "FOR4",
                                                   "FOR9",
                                                   "FOR10",
                                                   "FOR11"))

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

violins_muscle_disease <- data_plot |>
    dplyr::filter(project == "muscle_disease") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = anonimized_subject,
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

# ggplot2::ggsave(here::here("doc/figures/figure_6/violins_disease_participants.png"),
#                 height = 60,
#                 width = 128,
#                 units = "mm")

ggpubr::ggarrange(violins_heterofiber,
                  violins_muscle_disease,
                  legend = "right")

# ggplot2::ggsave(here::here("doc/figures/figure_1/violins_proteomics_participants.png"),
#                 height = 60,
#                 width = 195,
#                 units = "mm")
