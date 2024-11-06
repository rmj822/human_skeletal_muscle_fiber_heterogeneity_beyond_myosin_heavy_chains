
################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Matrix)
library(Seurat)
library(viridis)
library(rstatix)
library(ggnewscale)
library(ComplexUpset)
library(viridis)
library(ggpubr)
library(ggdist)
library(janitor)
library(DropletUtils)
library(ggrepel)
library(Hmisc)
library(grid)


#' Filtering missing values from rows
#'
#' @param .data dataset to filter
#' @param percentage_accepted_missing % of accepted missing values
#'
#' @return a dataframe
#' @export
#'
#' @examples
filtering_rows_Na <- function(.data, percentage_accepted_missing) {
    row_keep_vector <- .data |>
        is.na() |>
        rowSums()

    row_keep_vector <- row_keep_vector / ncol(.data)

    row_keep_vector <- row_keep_vector <= percentage_accepted_missing

    data_filtered <- .data |>
        tibble::add_column(row_keep_vector) |>
        dplyr::filter(row_keep_vector == T) |>
        dplyr::select(!starts_with("row"))

    return(data_filtered)
}


################################################################################################################################################
################################################       LOAD DATA      ########################################################################
################################################################################################################################################

# Load filtered Seurat object ---------------------------------------------
load("/Users/thibauxvds/Library/CloudStorage/OneDrive-UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq/8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

# Extract counts and metadata
counts <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")
counts <- as.data.frame(counts)

# Filter genes with % expression lower than 30 (not done, check if needed)
counts <- counts %>%
    tibble::rownames_to_column("Genes")

# Futher data prep
list_transcriptomics <- counts |>
    dplyr::group_split(
        Genes
    )

names(list_transcriptomics) <- counts |>
    dplyr::arrange(Genes) |>
    dplyr::pull(Genes)

gene_names <- counts |>
    dplyr::arrange(Genes) |>
    dplyr::pull(Genes)


################################################################################################################################################
########################################       Create correlation matrix      ##################################################################
################################################################################################################################################

correlation_function <- function(.list, .gene, .myh_to_correlate) {

    x <- .list[[.myh_to_correlate]] |>
        tibble::column_to_rownames("Genes") |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Genes") |>
        dplyr::pull(2, name = Genes)

    y <- .list[[.gene]] |>
        tibble::column_to_rownames("Genes") |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Genes") |>
        dplyr::pull(2, name = Genes)

    test <- cor.test(x,
                     y,
                     na.rm = TRUE,
                     method = "pearson")

    result <- data.frame(gene = .list[[.gene]]$Genes,
                         statistic = test$statistic,
                         correlation = test$estimate,
                         p_val = test$p.value) |>
        tibble::remove_rownames()

    return(result)
}

num_cores <- parallel::detectCores() - 1  # Adjust the number of cores as needed
doParallel::registerDoParallel(cores = num_cores)

correlation_matrix_MYH7 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_transcriptomics,
        .myh_to_correlate = "MYH7",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

correlation_matrix_MYH2 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_transcriptomics,
        .myh_to_correlate = "MYH2",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

correlation_matrix_MYH1 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_transcriptomics,
        .myh_to_correlate = "MYH1",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

# Stop parallel processing
doParallel::stopImplicitCluster()

################################################################################################################################################
#####################################################       MYH7     ###########################################################################
################################################################################################################################################

correlation_matrix_MYH7$log_pvalue <- -log10(correlation_matrix_MYH7$p_val)

correlation_matrix_MYH7 <- correlation_matrix_MYH7 %>%
    dplyr::mutate(
        log_pvalue = dplyr::case_when(
            log_pvalue == Inf ~ 349,
            TRUE ~ log_pvalue
        )
    )

labels_MYH7 <- c(correlation_matrix_MYH7 |>
                     dplyr::arrange(desc(correlation)) |>
                     dplyr::slice_head(n = 10) |>
                     dplyr::pull(gene),
                 correlation_matrix_MYH7 |>
                     dplyr::arrange(correlation) |>
                     dplyr::slice_head(n = 10) |>
                     dplyr::pull(gene)
)

correlation_matrix_MYH7 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.5 ~ "positive",
            fdr < 0.05 & correlation < -0.5 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = log_pvalue,
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH7 |>
            dplyr::mutate(
                label = dplyr::case_when(
                    gene %in% labels_MYH7 ~ gene,
                )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.5 ~ "positive",
                    fdr < 0.05 & correlation < -0.5 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = log_pvalue,
            fill = color,
            label = label
        ),
        max.overlaps = Inf,
        size = 2,
        color = "black",
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::ggtitle("Transcripts correlating with MYH7") +
    ggplot2::ylim(c(0, 350)) +
    ggplot2::xlim(c(-1.01, 1.01)) +
    ggplot2::scale_color_manual(
        values = c("#990000",
                   "lightgray",
                   "#0570b0")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#ebcccc",
        "lightgrey",
        "#cde2ef"
    )) +
    ggplot2::xlab("Pearson's r") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6),
        plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5)
    )

ggplot2::ggsave(here::here("doc/figures/figure_1/transcriptomics_correlation_MYH7.pdf"),
                units = "mm",
                height = 60,
                width = 60)


################################################################################################################################################
#####################################################       MYH2     ###########################################################################
################################################################################################################################################

correlation_matrix_MYH2$log_pvalue <- -log10(correlation_matrix_MYH2$p_val)

correlation_matrix_MYH2 <- correlation_matrix_MYH2 %>%
    dplyr::mutate(
        log_pvalue = dplyr::case_when(
            log_pvalue == Inf ~ 349,
            TRUE ~ log_pvalue
        )
    )

labels_MYH2 <- c(correlation_matrix_MYH2 |>
                     dplyr::arrange(desc(correlation)) |>
                     dplyr::slice_head(n = 10) |>
                     dplyr::pull(gene),
                 correlation_matrix_MYH2 |>
                     dplyr::arrange(correlation) |>
                     dplyr::slice_head(n = 10) |>
                     dplyr::pull(gene))

correlation_matrix_MYH2 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.5 ~ "positive",
            fdr < 0.05 & correlation < -0.5 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = log_pvalue,
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH2 |>
            dplyr::mutate(
                label = dplyr::case_when(
                    gene %in% labels_MYH2 ~ gene,
                )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.5 ~ "positive",
                    fdr < 0.05 & correlation < -0.5 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = log_pvalue,
            fill = color,
            label = label
        ),
        max.overlaps = Inf,
        size = 2,
        color = "black",
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::ggtitle("Transcripts correlating with MYH2") +
    ggplot2::ylim(c(0, 350)) +
    ggplot2::xlim(c(-1.01, 1.01)) +
    ggplot2::scale_color_manual(
        values = c("#990000",
                   "lightgray",
                   "#0570b0")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#ebcccc",
        "lightgrey",
        "#cde2ef"
    )) +
    ggplot2::xlab("Pearson's r") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6),
        plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5)
    )

ggplot2::ggsave(here::here("doc/figures/figure_1/transcriptomics_correlation_MYH2.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
#####################################################       MYH1     ###########################################################################
################################################################################################################################################

correlation_matrix_MYH1$log_pvalue <- -log10(correlation_matrix_MYH1$p_val)

correlation_matrix_MYH1 <- correlation_matrix_MYH1 %>%
    dplyr::mutate(
        log_pvalue = dplyr::case_when(
            log_pvalue == Inf ~ 349,
            TRUE ~ log_pvalue
        )
    )

labels_MYH1 <- c(correlation_matrix_MYH1 |>
                     dplyr::arrange(desc(correlation)) |>
                     dplyr::slice_head(n = 3) |>
                     dplyr::pull(gene),
                 correlation_matrix_MYH1 |>
                     dplyr::arrange(correlation) |>
                     dplyr::slice_head(n = 1) |>
                     dplyr::pull(gene))


correlation_matrix_MYH1 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.5 ~ "positive",
            fdr < 0.05 & correlation < -0.5 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = log_pvalue,
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH1 |>
            dplyr::mutate(
                label = dplyr::case_when(
                    gene %in% labels_MYH1 ~ gene,
                )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.5 ~ "positive",
                    fdr < 0.05 & correlation < -0.5 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = log_pvalue,
            fill = color,
            label = label
        ),
        max.overlaps = Inf,
        size = 2,
        color = "black",
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::ggtitle("Transcripts correlating with MYH1") +
    ggplot2::ylim(c(0, 350)) +
    ggplot2::xlim(c(-1.01, 1.01)) +
    ggplot2::scale_color_manual(
        values = c("#990000",
                   "lightgray",
                   "#0570b0")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#ebcccc",
        "lightgrey",
        "#cde2ef"
    )) +
    ggplot2::xlab("Pearson's r") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6),
        plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5)
    )

ggplot2::ggsave(here::here("doc/figures/figure_1/transcriptomics_correlation_MYH1.pdf"),
                units = "mm",
                height = 60,
                width = 60)

