
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

metadata_muscle_disease <- vroom::vroom(
  here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data-raw/muscle_disease/metadata_muscle_disease.txt")
)


# Loading Transcriptomics data --------------------------------------------

load(here::here("data-raw/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata"))

transcriptomics_data <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data$nGene |>
  as.data.frame()

colnames(transcriptomics_data) <- "number_of_proteins"

# Calculating number quantified proteins/study ----------------------------

heterofiber_number_quantified_proteins <- colSums(!is.na(proteomics_data)) |>
  as.data.frame()

colnames(heterofiber_number_quantified_proteins) <- "number_of_proteins"

heterofiber_number_quantified_proteins <- heterofiber_number_quantified_proteins |>
  dplyr::mutate(study = rep("heterofiber", ncol(proteomics_data)))

muscle_disease_number_quantified_proteins <- colSums(!is.na(muscle_disease_data)) |>
  as.data.frame() |>
  tibble::rownames_to_column("SampleID")

colnames(muscle_disease_number_quantified_proteins) <- c(
  "SampleID",
  "number_of_proteins"
)

muscle_disease_number_quantified_proteins <- metadata_muscle_disease |>
  dplyr::select(SampleID, condition) |>
  dplyr::inner_join(muscle_disease_number_quantified_proteins) |>
  tibble::column_to_rownames("SampleID") |>
  dplyr::rename("study" = "condition")

# Transcriptomics

transcriptomics_data <- transcriptomics_data |>
  dplyr::mutate(study = rep("transcriptomics", nrow(transcriptomics_data)))

# Merging and plotting ----------------------------------------------------

heterofiber_number_quantified_proteins |>
  dplyr::bind_rows(
    muscle_disease_number_quantified_proteins,
    transcriptomics_data
  ) |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::filter(!number_of_proteins < 1000) |>
  dplyr::mutate(project = dplyr::case_when(
    study == "heterofiber" ~ "proteomics",
    study == "transcriptomics" ~ "transcriptomics",
    TRUE ~ "muscle_disease"
  )) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = forcats::fct_infreq(study),
      y = number_of_proteins,
      fill = study
    )
  ) +
  ggplot2::geom_violin(alpha = 0.5) +
  ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3)) +
  ggplot2::scale_fill_manual("Study", values = c(
    "#E48C2AFF",
    "#969594",
    "#68b9ec",
    "#66db5b",
    "#d662c4"
  )) +
  ggplot2::theme_minimal()

# Making three plots and merging ------------------------------------------

heterofiber_violin <- heterofiber_number_quantified_proteins |>
    dplyr::bind_rows(
        muscle_disease_number_quantified_proteins,
        transcriptomics_data
    ) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::filter(!number_of_proteins < 1000) |>
    dplyr::mutate(project = dplyr::case_when(
        study == "heterofiber" ~ "proteomics",
        study == "transcriptomics" ~ "transcriptomics",
        TRUE ~ "muscle_disease"
    )) |>
    dplyr::filter(study ==  "heterofiber") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = study,
            y = number_of_proteins,
            fill = study
        )
    ) +
    ggplot2::geom_violin(alpha = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual("Study", values = c(
        "#68b9ec")) +
    ggplot2::ylim(1000, 2500) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1000 fiber proteome",
                  y = "Number of quantified proteins") +
    ggplot2::scale_x_discrete(breaks = c("heterofiber"),
                     labels = c("")) +
    ggplot2::theme(text = ggplot2::element_text(size = 8))

transcriptomics_violin <- heterofiber_number_quantified_proteins |>
    dplyr::bind_rows(
        muscle_disease_number_quantified_proteins,
        transcriptomics_data
    ) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::filter(!number_of_proteins < 1000) |>
    dplyr::mutate(project = dplyr::case_when(
        study == "heterofiber" ~ "proteomics",
        study == "transcriptomics" ~ "transcriptomics",
        TRUE ~ "muscle_disease"
    )) |>
    dplyr::filter(study ==  "transcriptomics") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = study,
            y = number_of_proteins,
            fill = study
        )
    ) +
    ggplot2::geom_violin(alpha = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual("Study", values = c(
        "#66db5b")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1000 fiber transcriptome",
                  y = "Number of quantified genes") +
    ggpubr::rremove("x.text")  +
    ggplot2::theme(text = ggplot2::element_text(size = 8))

muscle_disease_violin <- heterofiber_number_quantified_proteins |>
    dplyr::bind_rows(muscle_disease_number_quantified_proteins,
                     transcriptomics_data) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::filter(!number_of_proteins < 1000) |>
    dplyr::mutate(
        project = dplyr::case_when(
            study == "heterofiber" ~ "proteomics",
            study == "transcriptomics" ~ "transcriptomics",
            TRUE ~ "muscle_disease"
        )
    ) |>
    dplyr::filter(project ==  "muscle_disease") |>
    ggplot2::ggplot(ggplot2::aes(x = study,
                                 y = number_of_proteins,
                                 fill = study)) +
    ggplot2::geom_violin(alpha = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual("Study", values = c("#E48C2AFF",
                                                   "#969594",
                                                   "#d662c4")) +
    ggplot2::theme_minimal() +
    ggplot2::ylim(1000, 2500) +
    ggpubr::rremove("y.axis") +
    ggpubr::rremove("ylab") +
    ggpubr::rremove("y.text") +
    ggplot2::labs(x = "muscle disease") +
    ggplot2::theme(text = ggplot2::element_text(size = 8))

ggpubr::ggarrange(transcriptomics_violin,
                  heterofiber_violin,
                  muscle_disease_violin,
                  ncol = 3,
                  common.legend = FALSE,
                  legend = "none")

ggplot2::ggsave(here::here("doc/figures/figure_1/violin_n_quantified_genes_proteins.png"),
                height = 60,
                width = 120,
                units = "mm")
