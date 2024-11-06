
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

metadata <- vroom::vroom(here::here("data/metadata_proteomics_seurat_clusters.csv")) |>
    dplyr::rename("sample_id" = 1) |>
    dplyr::mutate(subject = dplyr::case_when(
        subject == "FOR2" ~ "P1",
        subject == "FOR4" ~ "P2",
        subject == "FOR9" ~ "P3",
        subject == "FOR10" ~ "P4",
        subject == "FOR11" ~ "P5",
        TRUE ~ "0"
    )) |>
    dplyr::mutate(fiber_type =
                      dplyr::case_when(
                          seurat_clusters == "3" ~ "slow",
                          seurat_clusters == "1" ~ "slow",
                          TRUE ~ "fast"
                      ))

proteomics_data <-vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    filtering_rows_Na(percentage_accepted_missing = 0) |>
    log2() |>
    as.data.frame() |>
    limma::normalizeBetweenArrays(method = "quantile") |>
    limma::removeBatchEffect(
        metadata$MS_batch
    ) |>
    as.data.frame()

data_heatmap <- proteomics_data |>
    t() |>
    scale() |>
    t() |>
    as.matrix()

color_annotations_fiber_type <- c(
    "slow" = "#440154FF",
    "fast" = "#5DC863FF"
)

color_annotations_subject <- c(
    "P1" = "#30123BFF",
    "P2" = "#28BBECFF",
    "P3" = "#A2FC3CFF",
    "P4" = "#FB8022FF",
    "P5" = "#7A0403FF"
)

color_annotations_batch <- c(
    "1" = "#0D0887FF",
    "2" = "#CC4678FF",
    "3" = "#F0F921FF"
)

annotations <- as.data.frame(colnames(data_heatmap)) |>
    dplyr::rename("sample_id" = 1) |>
    dplyr::inner_join(metadata |>
                          dplyr::select(c(fiber_type,
                                          subject,
                                          # MS_batch,
                                          # cluster_name,
                                          sample_id))) |>
    # dplyr::mutate(MS_batch = as.factor(MS_batch)) |>
    tibble::column_to_rownames("sample_id")

num_cores <- detectCores() - 1  # Adjust the number of cores as needed
registerDoParallel(cores = num_cores)

heatmap <- pheatmap::pheatmap(
    color = colorRampPalette(c("darkblue", "white", "red"))(100),
    # color = viridisLite::viridis(n = 100),
    breaks = seq(from = -5, to = 5, by = 0.1),
    data_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # legend_breaks = c(1, 0.9825, 0.965),
    fontsize = 6,
    show_colnames = F,
    show_rownames = F,
    annotation_col = annotations,
    annotation_colors = list(fiber_type = color_annotations_fiber_type,
                             subject = color_annotations_subject
                             # MS_batch = color_annotations_batch
                             ),
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    treeheight_row = 20,
    treeheight_col = 20
)
ggplotify::as.ggplot(heatmap)

ggplot2::ggsave(here::here("doc/figures/figure_2/full_heatmap_proteomics.png"),
                units = "mm",
                height = 60,
                width = 120)

# Stop parallel processing
stopImplicitCluster()

protein_clusters <- cutree(heatmap$tree_row, 5) |>
    as.data.frame()


# Correlation heatmap -----------------------------------------------------

cor_matrix <- cor(proteomics_data, use = "complete.obs")

color_annotations_fiber_type <- c(
    "slow" = "#440154FF",
    "fast" = "#5DC863FF"
)

color_annotations_subject <- c(
    "P1" = "#30123BFF",
    "P2" = "#28BBECFF",
    "P3" = "#A2FC3CFF",
    "P4" = "#FB8022FF",
    "P5" = "#7A0403FF"
)

annotations <- as.data.frame(colnames(cor_matrix)) |>
    dplyr::rename("sample_id" = 1) |>
    dplyr::inner_join(metadata |>
                          dplyr::select(c(fiber_type,
                                          subject,
                                          # cluster_name,
                                          sample_id))) |>
    tibble::column_to_rownames("sample_id")

heatmap <- pheatmap::pheatmap(
    color = colorRampPalette(c("#fffdfc", "#4292c6", "#03004d"))(200),
    cor_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # legend_breaks = c(1, 0.9825, 0.965),
    fontsize = 6,
    show_colnames = F,
    show_rownames = F,
    annotation_col = annotations,
    annotation_colors = list(fiber_type = color_annotations_fiber_type,
                             subject = color_annotations_subject),
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    treeheight_row = 20,
    treeheight_col = 20
)

ggplotify::as.ggplot(heatmap)

ggplot2::ggsave(here::here("doc/figures/figure_2/correlation_heatmap_proteomics.png"),
                units = "mm",
                height = 60,
                width = 90)


# Extract clusters --------------------------------------------------------

correlation_clusters <- cutree(heatmap$tree_row, 2) |>
    as.data.frame() |>
    dplyr::rename("correlation_clusters" = 1) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::inner_join(metadata)


# Slow cluster ------------------------------------------------------------

slow_metadata <- correlation_clusters |>
    dplyr::filter(correlation_clusters == 2)

data_heatmap_slow <-
    vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    filtering_rows_Na(percentage_accepted_missing = 0.3) |>
    log2() |>
    as.data.frame() |>
    limma::normalizeBetweenArrays(method = "quantile") |>
    limma::removeBatchEffect(metadata$MS_batch) |>
    as.data.frame() |>
    # dplyr::select(slow_metadata$sample_id) |>
    t() |>
    scale() |>
    t() |>
    as.data.frame()
    # dplyr::select(slow_metadata$sample_id)

heatmap_slow <- pheatmap::pheatmap(
    mat = data_heatmap_slow,
    color = colorRampPalette(c("darkblue", "white", "red"))(100),,
    breaks = seq(from = -5, to = 5, by = 0.1),
    show_rownames = T,
    show_colnames = F,
    cluster_cols = T,
    fontsize = 4,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2"
    # cutree_rows = 5
)
