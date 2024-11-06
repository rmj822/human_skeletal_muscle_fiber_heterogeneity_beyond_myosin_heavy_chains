proteomics_data <- vroom::vroom(
    here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_proteomics_filtered.csv"),
    col_select = !c(1)
) |>
    as.data.frame() |>
    tibble::column_to_rownames("Gene.name")

data_log <- proteomics_data |>
    log10()

metadata <- vroom::vroom(
  "C:/Users/jns822/Desktop/Scripts/Heterofiber/data-raw/metadata.txt"
)


# Making the function -----------------------------------------------------


vln_percentage_protein <- function(.data,
                                   .metadata,
                                   protein_of_interest) {
  protein <- as.character(protein_of_interest)

  sum_of_intensities <- colSums(proteomics_data,
    na.rm = TRUE
  )

  perc_expression_protein <- proteomics_data |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name == protein) |>
    imputeTS::na_replace(0) |>
    tibble::column_to_rownames("Gene.name")

  perc_expression_protein <- perc_expression_protein / sum_of_intensities * 100

  metadata_subject <- metadata |>
    dplyr::select(subject, fiberID)

  perc_expression_protein |>
    t() |>
    as.data.frame() |>
    dplyr::rename(
      "percentage" = protein
    ) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(metadata_subject) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = subject,
        y = percentage,
        fill = subject
      )
    ) +
    ggplot2::geom_violin(alpha = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3)) +
    ggplot2::scale_fill_viridis_d("Subject", option = "turbo") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste("Percentage of Expression of", protein)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 13, face = "bold")) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(vjust = -0.35),
      axis.title.y = ggplot2::element_text(vjust = 0.35)
    ) +
    ggplot2::labs(x = "Subject", y = "Percentage")
}


# Percentage of expression satellite cells markers: -----------------------

# PAX7 was not found in the dataset

# NCAM1

perc_NCAM1 <- vln_percentage_protein(
  .data = proteomics_data,
  .metadata = metadata,
  protein_of_interest = "NCAM1"
)

ggplot2::ggsave(perc_NCAM1,
  device = "png",
  filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_NCAM1.png")
)


# Percentage of expression macrophage markers ------------------------------------

# MRC 1 and C1QA were not present in the dataset

# Percentage of expression FAP cell markers -------------------------------

# DCN

perc_DCN <- vln_percentage_protein(
    .data = proteomics_data,
    .metadata = metadata,
    protein_of_interest = "DCN"
)

ggplot2::ggsave(perc_DCN,
                device = "png",
                filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_DCN.png")
)

# CD34

perc_CD34 <- vln_percentage_protein(
    .data = proteomics_data,
    .metadata = metadata,
    protein_of_interest = "CD34"
)

ggplot2::ggsave(perc_CD34,
                device = "png",
                filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_CD34.png")
)

# Percentage of expression endothelial cells markers: ---------------------

# CDH5

perc_CDH5 <- vln_percentage_protein(
    .data = proteomics_data,
    .metadata = metadata,
    protein_of_interest = "CDH5"
)

ggplot2::ggsave(perc_CDH5,
                device = "png",
                filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_CDH5.png")
)

# PECAM1

perc_PECAM1 <- vln_percentage_protein(
    .data = proteomics_data,
    .metadata = metadata,
    protein_of_interest = "PECAM1"
)

ggplot2::ggsave(perc_PECAM1,
                device = "png",
                filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_PECAM1.png")
)

# Percentage of expression smooth muscle cells markers --------------------

# ACTA2

perc_ACTA2 <- vln_percentage_protein(
    .data = proteomics_data,
    .metadata = metadata,
    protein_of_interest = "ACTA2"
)

ggplot2::ggsave(perc_ACTA2,
                device = "png",
                filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_ACTA2.png")
)

# MYH11

perc_MYH11 <- vln_percentage_protein(
    .data = proteomics_data,
    .metadata = metadata,
    protein_of_interest = "MYH11"
)

ggplot2::ggsave(perc_MYH11,
                device = "png",
                filename = here::here("doc/figures/figure_1/violins_non_muscle_markers/violin_MYH11.png")
)
