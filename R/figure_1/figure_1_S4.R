
################################################################################################################################################
########################################################       FIGURE 1 S4A     ###################################################################
################################################################################################################################################

# Load filtered Seurat object ---------------------------------------------
load(here::here("data/figure_4/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata"))

# Assign identity of resolution 0.7 to clusters ---------------------------------------------
Idents(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"

filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data$final_cluster

# Assign slow vs fast clustering
metadata <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data |>
    dplyr::mutate(fiber_type_seurat = dplyr::case_when(
        final_cluster == "Slow1" ~ "slow",
        final_cluster == "Slow2" ~ "slow",
        final_cluster == "Fast1" ~ "fast",
        final_cluster == "Fast2" ~ "fast",
        final_cluster == "Fast3" ~ "fast",
        final_cluster == "Intermediate1" ~ "hybrid",
        TRUE ~ "NA"
    ))


# ACTN3 -------------------------------------------------------------------

feature_plot_ACTN3 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                          features = c("ACTN3"),
                                          reduction = "umap",
                                          pt.size = 0.3,
                                          label.size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c("Norm.\ncounts", option = "plasma") +
    ggplot2::ggtitle("UMAP transcriptomics\nby ACTN3 expression") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.1, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

# MYLK2 -------------------------------------------------------------------

feature_plot_MYLK2 <- Seurat::FeaturePlot(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                          features = c("MYLK2"),
                                          reduction = "umap",
                                          pt.size = 0.3,
                                          label.size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c("Norm.\ncounts", option = "plasma") +
    ggplot2::ggtitle("UMAP transcriptomics\nby MYLK2 expression") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.1, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        plot.margin = grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

patchwork::wrap_plots(feature_plot_ACTN3, feature_plot_MYLK2)

ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4A.png"),
       width = 99,
       height = 45,
       units="mm")

################################################################################################################################################
########################################################       FIGURE 1 S4B     ###################################################################
################################################################################################################################################

data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins
data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")
metadata <- vroom::vroom(
    here::here("data/metadata_proteomics_fiber_type.csv")
) |>
    dplyr::select(!1)
seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                              meta.data = metadata)
# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                                selection.method = "vst")
# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)
#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,  features = Seurat::VariableFeatures(object = seurat_proteome))
# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)
# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))
# Dimensionality reduction
# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6)

# Feature plot ACTN3 ------------------------------------------------------

feature_plot_ACTN3 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("ACTN3"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_ACTN3 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("ACTN3"),
                                          pt.size = 0.3,
                                          label.size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c("LFQ\nintensity\n(log2)", option = "plasma") +
    ggplot2::ggtitle("UMAP proteomics\nby ACTN3 abundance") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.1, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

# MYLK2 -------------------------------------------------------------------

feature_plot_MYLK2 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("MYLK2"),
                                          pt.size = 0.3,
                                          label.size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c("LFQ\nintensity\n(log2)", option = "plasma") +
    ggplot2::ggtitle("UMAP proteomics\nby MYLK2 abundance") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.1, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        plot.margin = grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

patchwork::wrap_plots(feature_plot_ACTN3, feature_plot_MYLK2)

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4B.pdf"),
                units = "mm",
                height = 45,
                width = 99)

################################################################################################################################################
########################################################       FIGURE 1 S4C, D and E    ###################################################################
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

counts <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")
counts <- as.data.frame(counts)

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
#####################################################       Panel C     ###########################################################################
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

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4C.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
#####################################################       panel D     ###########################################################################
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

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4D.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
#####################################################       Panel E     ###########################################################################
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

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4E.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
########################################################       FIGURE 1 S4F, G and H    ###################################################################
################################################################################################################################################

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


# Load data ---------------------------------------------------------------

data_proteomics <-vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    filtering_rows_Na(percentage_accepted_missing = 0.3) |>
    log2() |>
    as.data.frame() |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate(
        Genes = gsub(
            pattern = "-",
            replacement = ".",
            Genes
        )
    )

metadata <- vroom::vroom(here::here("data/metadata_proteomics_seurat_clusters.csv"))|>
    # dplyr::rename("fiberID" = 1) |>
    tibble::column_to_rownames("fiberID")

list_proteomics <- data_proteomics |>
    dplyr::group_split(
        Genes
    )

names(list_proteomics) <- data_proteomics |>
    dplyr::arrange(Genes) |>
    dplyr::pull(Genes)

gene_names <- data_proteomics |>
    dplyr::arrange(Genes) |>
    dplyr::pull(Genes)

# Create correlation matrix -----------------------------------------------

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


doParallel::registerDoParallel(cores = num_cores)

correlation_matrix_MYH7 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_proteomics,
        .myh_to_correlate = "MYH7",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

correlation_matrix_MYH2 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_proteomics,
        .myh_to_correlate = "MYH2",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

correlation_matrix_MYH1 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_proteomics,
        .myh_to_correlate = "MYH1",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

# Stop parallel processing
doParallel::stopImplicitCluster()

################################################################################################################################################
########################################################       FIGURE 1 S4F   ###################################################################
################################################################################################################################################

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
            fdr < 0.05 & correlation > 0.7 ~ "positive",
            fdr < 0.05 & correlation < -0.7 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
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
                    fdr < 0.05 & correlation > 0.7 ~ "positive",
                    fdr < 0.05 & correlation < -0.7 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
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
    ggplot2::ggtitle("Proteins correlating with MYH7") +
    ggplot2::ylim(c(0, 350)) +
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

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4F.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
########################################################       FIGURE 1 S4G   ###################################################################
################################################################################################################################################

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
            fdr < 0.05 & correlation > 0.7 ~ "positive",
            fdr < 0.05 & correlation < -0.7 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
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
                    fdr < 0.05 & correlation > 0.7 ~ "positive",
                    fdr < 0.05 & correlation < -0.7 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
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
    ggplot2::ggtitle("Proteins correlating with MYH2") +
    ggplot2::ylim(c(0, 350)) +
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

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4G.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
########################################################       FIGURE 1 S4H   ###################################################################
################################################################################################################################################

correlation_matrix_MYH1 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.7 ~ "positive",
            fdr < 0.05 & correlation < -0.7 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH1 |>
            dplyr::mutate(
                label = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.7 ~ gene,
                    fdr < 0.05 & correlation < -0.7 ~ gene,
                    TRUE ~ ""
                )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.7 ~ "positive",
                    fdr < 0.05 & correlation < -0.7 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
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
    ggplot2::ggtitle("Proteins correlating with MYH1") +
    ggplot2::ylim(c(0, 350)) +
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

ggplot2::ggsave(here::here("doc/figures/figure_1_S4/figure_1_S4H.pdf"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
########################################################       FIGURE 1 S4I ###################################################################
################################################################################################################################################

# Get counts data
data_transcriptomics <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                     assay = "SCT",
                                     slot = "counts")

data_transcriptomics <- as.data.frame(data_transcriptomics)

# Extract MYH counts
MYH1 <- c("MYH1")
counts_MYH1 <- data_transcriptomics[MYH1, ]
counts_MYH1 <- as.data.frame(as.matrix(counts_MYH1))
counts_MYH1 <- t(counts_MYH1)

MYH2 <- c("MYH2")
counts_MYH2 <- data_transcriptomics[MYH2, ]
counts_MYH2 <- as.data.frame(as.matrix(counts_MYH2))
counts_MYH2 <- t(counts_MYH2)

MYH7 <- c("MYH7")
counts_MYH7 <- data_transcriptomics[MYH7, ]
counts_MYH7 <- as.data.frame(as.matrix(counts_MYH7))
counts_MYH7 <- as.data.frame(t(counts_MYH7))
counts_MYH7$Row.names <- rownames(counts_MYH7)

# Merge MYH dataframes
counts_MYH <- merge(counts_MYH1, counts_MYH2, by="row.names")
counts_MYH <- merge(counts_MYH, counts_MYH7, by="Row.names")

# Delete Row.names column
rownames(counts_MYH) <- counts_MYH$Row.names
counts_MYH <- counts_MYH %>% dplyr::select(MYH1, MYH2, MYH7)

# Calculate percentages for each MYH
counts_MYH$sum_MYH <- counts_MYH$MYH7 + counts_MYH$MYH2 + counts_MYH$MYH1

counts_MYH$MYH7_fraction <- counts_MYH$MYH7 / counts_MYH$sum_MYH * 100
counts_MYH$MYH2_fraction <- counts_MYH$MYH2 / counts_MYH$sum_MYH * 100
counts_MYH$MYH1_fraction <- counts_MYH$MYH1 / counts_MYH$sum_MYH * 100

# Calculate difference between MYH2 and MYH1
counts_MYH$MYH2_MYH1_diff <- counts_MYH$MYH2_fraction - counts_MYH$MYH1_fraction
counts_MYH$MYHs_fast <- counts_MYH$MYH2 + counts_MYH$MYH1

# Extract UMAP coordinates and add MYH diff colouring
df_UMAP_blended <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest[["umap"]]@cell.embeddings %>%
    as.data.frame() %>%
    rownames_to_column("row.names")

df_UMAP_blended <- counts_MYH %>%
    dplyr::select(MYH2_MYH1_diff, MYHs_fast) %>%
    rownames_to_column("row.names") %>%
    left_join(df_UMAP_blended, by = "row.names")

# UMAP coloured by MYH
UMAP_blended <- df_UMAP_blended %>%
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color =MYH2_MYH1_diff)
    ) +
    ggplot2::geom_point(aes(alpha = MYHs_fast),
                        size = 1) +
    scale_alpha(range = c(0.2, 1),
                guide = "none") +
    ggplot2::scale_color_gradient2(
        "MYH2 % \n     - \nMYH1 %",
        low = "#D2631C",
        mid = "lightgrey",
        high = "#23CE6B",
        midpoint = 0,
        limits = c(-100, 100)
    ) +
    ggplot2::ggtitle("UMAP transcriptomics - by fast MYH blended expression") +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6, face = "bold"),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, "mm"),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 7, face = "bold")
        # plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm")
    )

perc_MYHs <- counts_MYH |>
    rownames_to_column("fiber_ID") %>%
    dplyr::select(c(fiber_ID, MYH7_fraction, MYH2_fraction, MYH1_fraction)) |>
    dplyr::mutate(
        dplyr::across(
            .cols = !fiber_ID,
            ~ .x/100
        )
    )

perc_list <- perc_MYHs |>
    dplyr::group_split(
        fiber_ID
    )

names(perc_list) <- perc_MYHs |>
    dplyr::arrange(fiber_ID) |>
    dplyr::pull(fiber_ID)

rgb_maker <- function(.data, .fiber_ID){
    tmp2 <- .data[[.fiber_ID]]

    tmp2 <- tmp2 |>
        dplyr::mutate(
            color = rgb(
                red = tmp2$MYH1_fraction,
                green = tmp2$MYH2_fraction,
                blue = tmp2$MYH7_fraction,
                alpha = 0.7
            ))

    return(tmp2)
}

color_list <- purrr::map(
    names(perc_list),
    ~ rgb_maker(
        .data = perc_list,
        .fiber_ID = .x
    )
) |>
    purrr::list_rbind() |>
    tibble::column_to_rownames("fiber_ID")

# Extract UMAP coordinates and add MYH diff colouring
df_UMAP_blended_rgb <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest[["umap"]]@cell.embeddings %>%
    as.data.frame() %>%
    rownames_to_column("row.names")

df_UMAP_blended_rgb <- color_list %>%
    dplyr::select(color) %>%
    rownames_to_column("row.names") %>%
    left_join(df_UMAP_blended_rgb, by = "row.names")

df_UMAP_blended_rgb$color


# UMAP coloured by MYH RGB
df_UMAP_blended_rgb |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = color)
    ) +
    ggplot2::geom_point(
        # ggplot2::aes(alpha = coloring$MYHs_fast),
        size = 1,
        alpha = 0.65) +
    ggplot2::scale_color_identity() +
    # ggplot2::scale_alpha(range = c(0.2, 1),
    #                      guide = "none") +
    # ggplot2::scale_color_gradient2(
    #     "MYH2 % \n     - \nMYH1 %",
    #     low = "#D2631C",
    #     mid = "lightgrey",
    #     high = "#23CE6B",
    #     midpoint = 0,
    #     limits = c(-100, 100),
    #     breaks = c(-100, -50, 0, 50, 100)
    # ) +
    ggplot2::ggtitle("UMAP transcriptomics - by MYH blended expression") +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6, face = "bold"),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, "mm"),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 7, face = "bold")
        # plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm")
    )

ggplot2::ggsave(
    here::here("doc/figures/figure_1_S4/figure_1_S4I.pdf"),
    units = "mm",
    height = 60,
    width = 80
)

################################################################################################################################################
########################################################       FIGURE 1 S4J ###################################################################
################################################################################################################################################
perc_MYHs <- vroom::vroom(here::here("data/perc_MYH_proteomics/perc_MYH_proteomics.csv")) |>
    dplyr::select(!1)

# Replace patterns in fiberID
perc_MYHs$fiber_ID <- gsub("FOR2", "P1", perc_MYHs$fiber_ID)
perc_MYHs$fiber_ID <- gsub("FOR4", "P2", perc_MYHs$fiber_ID)
perc_MYHs$fiber_ID <- gsub("FOR9", "P3", perc_MYHs$fiber_ID)
perc_MYHs$fiber_ID <- gsub("FOR10", "P4", perc_MYHs$fiber_ID)
perc_MYHs$fiber_ID <- gsub("FOR11", "P5", perc_MYHs$fiber_ID)

perc_MYHs <- perc_MYHs |>
    dplyr::select(c(fiber_ID, MYHs, values)) |>
    dplyr::filter(fiber_ID %in% colnames(data_proteomics)) |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    dplyr::mutate(
        dplyr::across(
            .cols = !fiber_ID,
            ~ .x/100
        )
    )

perc_list <- perc_MYHs |>
    dplyr::group_split(
        fiber_ID
    )

names(perc_list) <- perc_MYHs |>
    dplyr::arrange(fiber_ID) |>
    dplyr::pull(fiber_ID)

rgb_maker <- function(.data, .fiber_ID){
    tmp2 <- .data[[.fiber_ID]]

    tmp2 <- tmp2 |>
        dplyr::mutate(
            color = rgb(
                red = tmp2$MYH1,
                green = tmp2$MYH2,
                blue = tmp2$MYH7,
                alpha = 0.7
            ))

    return(tmp2)
}

color_list <- purrr::map(
    names(perc_list),
    ~ rgb_maker(
        .data = perc_list,
        .fiber_ID = .x
    )
) |>
    purrr::list_rbind() |>
    tibble::column_to_rownames("fiber_ID")

seurat_proteome[["umap"]]@cell.embeddings |>
    as.data.frame() |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::arrange(fiber_ID) |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = color_list$color)
    ) +
    ggplot2::geom_point(
        # ggplot2::aes(alpha = coloring$MYHs_fast),
        size = 1,
        alpha = 0.65) +
    ggplot2::scale_color_identity() +
    # ggplot2::scale_alpha(range = c(0.2, 1),
    #                      guide = "none") +
    # ggplot2::scale_color_gradient2(
    #     "MYH2 % \n     - \nMYH1 %",
    #     low = "#D2631C",
    #     mid = "lightgrey",
    #     high = "#23CE6B",
    #     midpoint = 0,
    #     limits = c(-100, 100),
    #     breaks = c(-100, -50, 0, 50, 100)
    # ) +
    ggplot2::ggtitle("UMAP proteomics - by MYH blended expression") +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6, face = "bold"),
        legend.position = "none",
        legend.key.size = ggplot2::unit(2, "mm"),
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 7, face = "bold")
        # plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm")
    )

ggplot2::ggsave(
    here::here("doc/figures/figure_1_S4/figure_1_S4J.pdf"),
    units = "mm",
    height = 60,
    width = 80
)
