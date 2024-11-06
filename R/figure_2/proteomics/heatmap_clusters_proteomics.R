
# Loading data ------------------------------------------------------------

proteomics_data <- read.csv(
  here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_proteomics_filtered.csv")
) |>
  as.data.frame() |>
  dplyr::select(!X) |>
  tibble::column_to_rownames("Gene.name")

metadata <- vroom::vroom(
  "data/metadata_proteomics.csv"
) |>
  dplyr::rename(fiberID = 1)

seurat_clusters <- vroom::vroom(
  "data/proteomics clustering/seurat_clusters_6PC_res04.csv"
)

metadata_proteomics <- metadata |>
  dplyr::inner_join(seurat_clusters)


# Data preparation and pseudobulk -----------------------------------------

data_proteomics <- proteomics_data |>
  log2()

pseudobulk_maker <- function(.data, metadata, subject_id, grouping, colname) {
  selection_vector <- metadata |>
    dplyr::filter(seurat_clusters == grouping) |>
    dplyr::filter(subject == subject_id) |>
    dplyr::pull("fiberID")

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

pseudobulking <- function() {
  return(
    data_pseudobulk <- data.frame(
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR2",
        grouping = 0,
        colname = "FOR2_cluster_0"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR2",
        grouping = 1,
        colname = "FOR2_cluster_1"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR2",
        grouping = 2,
        colname = "FOR2_cluster_2"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR2",
        grouping = 3,
        colname = "FOR2_cluster_3"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR2",
        grouping = 4,
        colname = "FOR2_cluster_4"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR2",
        grouping = 5,
        colname = "FOR2_cluster_5"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR4",
        grouping = 0,
        colname = "FOR4_cluster_0"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR4",
        grouping = 1,
        colname = "FOR4_cluster_1"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR4",
        grouping = 2,
        colname = "FOR4_cluster_2"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR4",
        grouping = 3,
        colname = "FOR4_cluster_3"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR4",
        grouping = 4,
        colname = "FOR4_cluster_4"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR4",
        grouping = 5,
        colname = "FOR4_cluster_5"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR9",
        grouping = 0,
        colname = "FOR9_cluster_0"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR9",
        grouping = 1,
        colname = "FOR9_cluster_1"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR9",
        grouping = 2,
        colname = "FOR9_cluster_2"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR9",
        grouping = 3,
        colname = "FOR9_cluster_3"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR9",
        grouping = 4,
        colname = "FOR9_cluster_4"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR9",
        grouping = 5,
        colname = "FOR9_cluster_5"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR10",
        grouping = 0,
        colname = "FOR10_cluster_0"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR10",
        grouping = 1,
        colname = "FOR10_cluster_1"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR10",
        grouping = 2,
        colname = "FOR10_cluster_2"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR10",
        grouping = 3,
        colname = "FOR10_cluster_3"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR10",
        grouping = 4,
        colname = "FOR10_cluster_4"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR10",
        grouping = 5,
        colname = "FOR10_cluster_5"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR11",
        grouping = 0,
        colname = "FOR11_cluster_0"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR11",
        grouping = 1,
        colname = "FOR11_cluster_1"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR11",
        grouping = 2,
        colname = "FOR11_cluster_2"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR11",
        grouping = 3,
        colname = "FOR11_cluster_3"
      ),
      pseudobulk_maker(
        .data = as.data.frame(data_proteomics),
        metadata = metadata_proteomics,
        subject_id = "FOR11",
        grouping = 4,
        colname = "FOR11_cluster_4"
      )
      # pseudobulk_maker(
      #     .data = as.data.frame(data_proteomics),
      #     metadata = metadata_proteomics,
      #     subject_id = "FOR11",
      #     grouping = 5,
      #     colname = "FOR11_cluster_5"
      # )
    )
  )
}

data_pseudobulk <- pseudobulking()

data_grouping <- data.frame(
  "sample_id" = colnames(data_pseudobulk),
  "cluster" = c(
    rep(
      c(
        "cluster_0",
        "cluster_1",
        "cluster_2",
        "cluster_3",
        "cluster_4",
        "cluster_5"
      ),
      times = 4
    ),
    "cluster_0",
    "cluster_1",
    "cluster_2",
    "cluster_3",
    "cluster_4"
  ),
  "subject" = c(
    rep("FOR2", times = 6),
    rep("FOR4", times = 6),
    rep("FOR9", times = 6),
    rep("FOR10", times = 6),
    rep("FOR11", times = 5)
  )
)


# Limma repeated measures ANOVA -------------------------------------------

data_limma <- data_pseudobulk |>
    limma::normalizeBetweenArrays()

vector_factor_clusters <- factor(data_grouping$cluster,
                                 levels = c(
                                     "cluster_0",
                                     "cluster_1",
                                     "cluster_2",
                                     "cluster_3",
                                     "cluster_4",
                                     "cluster_5"
                                 )
)

matrix_design <- model.matrix(
    ~ 0 + vector_factor_clusters + subject,
    data_grouping
)

colnames(matrix_design) <- c(
    "cluster_0",
    "cluster_1",
    "cluster_2",
    "cluster_3",
    "cluster_4",
    "cluster_5",
    "FOR11",
    "FOR2",
    "FOR4",
    "FOR9"
)

contrast_matrix <- limma::makeContrasts(
    cluster_0 - cluster_5,
    cluster_0 - cluster_2,
    cluster_4 - cluster_0,
    cluster_1 - cluster_4,
    cluster_3 - cluster_1,
    levels = matrix_design
)

fit <- limma::lmFit(
    data_limma,
    matrix_design
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp <- limma::eBayes(tmp)

F_table <- limma::topTable(tmp,
                           coef = NULL,
                           p.value = 0.05,
                           n = Inf,
                           lfc = 2)

colnames(F_table) <- paste("F", colnames(F_table), sep = "_")

F_table <- F_table |>
    tibble::rownames_to_column("Gene.name")

# I produced a heatmap of the > 2 logFC proteins and then I manually selected 5 proteins from each cluster:


F_table <- F_table |>
    dplyr::filter(Gene.name %in% c(
        "MYH7",
        "TNNT1",
        "GOLGA4",
        "CASQ2",
        "PDLIM1",
        "PHIP",
        "UGDH",
        "AK4",
        "KLHL36",
        "RPL32",
        "MYH2",
        "USP28",
        "TNNT3",
        "KIF1C",
        "DNAH7",
        "HIST1H1B",
        "MYL4",
        "SART3",
        "FAM213A",
        "COPG1",
        "HIST1H2AB",
        "RPL35",
        "KDELC2",
        "PIP",
        "TPPP"
    ))

# Making the heatmap ------------------------------------------------------

data_heatmap <- data_pseudobulk |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% F_table$Gene.name) |>
    dplyr::inner_join(F_table) |>
    tibble::column_to_rownames("Gene.name")

data_grouped_subject <- list(
    data_scale_FOR2 <- data_heatmap |>
        dplyr::select(dplyr::contains("FOR2")) |>
        t() |>
        scale() |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Gene.name"),
    data_scale_FOR4 <- data_heatmap |>
        dplyr::select(dplyr::contains("FOR4")) |>
        t() |>
        scale() |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Gene.name"),
    data_scale_FOR9 <- data_heatmap |>
        dplyr::select(dplyr::contains("FOR9")) |>
        t() |>
        scale() |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Gene.name"),
    data_scale_FOR10 <- data_heatmap |>
        dplyr::select(dplyr::contains("FOR10")) |>
        t() |>
        scale() |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Gene.name"),
    data_scale_FOR11 <- data_heatmap |>
        dplyr::select(dplyr::contains("FOR11")) |>
        t() |>
        scale() |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("Gene.name")
)

names(data_grouped_subject) <- c(
    "data_scale_FOR2",
    "data_scale_FOR4",
    "data_scale_FOR9",
    "data_scale_FOR10",
    "data_scale_FOR11"
)

data_scaled <- cbind(
    data_grouped_subject$data_scale_FOR2,
    data_grouped_subject$data_scale_FOR4,
    data_grouped_subject$data_scale_FOR9,
    data_grouped_subject$data_scale_FOR10,
    data_grouped_subject$data_scale_FOR11,
    make.row.names = "Gene.name"
) |>
    tibble::column_to_rownames("Gene.name") |>
    dplyr::select(dplyr::contains("FOR")) |>
    t() |>
    scale() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_name") |>
    dplyr::mutate(
        "cluster" = dplyr::case_when(
            sample_name == "FOR2_cluster_0" ~ 0,
            sample_name == "FOR4_cluster_0" ~ 0,
            sample_name == "FOR9_cluster_0" ~ 0,
            sample_name == "FOR10_cluster_0" ~ 0,
            sample_name == "FOR11_cluster_0" ~ 0,
            sample_name == "FOR2_cluster_1" ~ 1,
            sample_name == "FOR4_cluster_1" ~ 1,
            sample_name == "FOR9_cluster_1" ~ 1,
            sample_name == "FOR10_cluster_1" ~ 1,
            sample_name == "FOR11_cluster_1" ~ 1,
            sample_name == "FOR2_cluster_2" ~ 2,
            sample_name == "FOR4_cluster_2" ~ 2,
            sample_name == "FOR9_cluster_2" ~ 2,
            sample_name == "FOR10_cluster_2" ~ 2,
            sample_name == "FOR11_cluster_2" ~ 2,
            sample_name == "FOR2_cluster_3" ~ 3,
            sample_name == "FOR4_cluster_3" ~ 3,
            sample_name == "FOR9_cluster_3" ~ 3,
            sample_name == "FOR10_cluster_3" ~ 3,
            sample_name == "FOR11_cluster_3" ~ 3,
            sample_name == "FOR2_cluster_4" ~ 4,
            sample_name == "FOR4_cluster_4" ~ 4,
            sample_name == "FOR9_cluster_4" ~ 4,
            sample_name == "FOR10_cluster_4" ~ 4,
            sample_name == "FOR11_cluster_4" ~ 4,
            sample_name == "FOR2_cluster_5" ~ 5,
            sample_name == "FOR4_cluster_5" ~ 5,
            sample_name == "FOR9_cluster_5" ~ 5,
            sample_name == "FOR10_cluster_5" ~ 5,
            sample_name == "FOR11_cluster_5" ~ 5
        )
    )

data_grouped_cluster <- list(
    data_scaled_cluster0 <- data_scaled |>
        dplyr::filter(cluster == 0) |>
        dplyr::select(!cluster) |>
        tibble::column_to_rownames("sample_name") |>
        t() |>
        as.data.frame(),
    data_scaled_cluster1 <- data_scaled |>
        dplyr::filter(cluster == 1) |>
        dplyr::select(!cluster) |>
        tibble::column_to_rownames("sample_name") |>
        t() |>
        as.data.frame(),
    data_scaled_cluster2 <- data_scaled |>
        dplyr::filter(cluster == 2) |>
        dplyr::select(!cluster) |>
        tibble::column_to_rownames("sample_name") |>
        t() |>
        as.data.frame(),
    data_scaled_cluster3 <- data_scaled |>
        dplyr::filter(cluster == 3) |>
        dplyr::select(!cluster) |>
        tibble::column_to_rownames("sample_name") |>
        t() |>
        as.data.frame(),
    data_scaled_cluster4 <- data_scaled |>
        dplyr::filter(cluster == 4) |>
        dplyr::select(!cluster) |>
        tibble::column_to_rownames("sample_name") |>
        t() |>
        as.data.frame(),
    data_scaled_cluster5 <- data_scaled |>
        dplyr::filter(cluster == 5) |>
        dplyr::select(!cluster) |>
        tibble::column_to_rownames("sample_name") |>
        t() |>
        as.data.frame()
)

names(data_grouped_cluster) <- c(
    "data_scaled_cluster0",
    "data_scaled_cluster1",
    "data_scaled_cluster2",
    "data_scaled_cluster3",
    "data_scaled_cluster4",
    "data_scaled_cluster5"
)

data_heatmap <- cbind(
    data_grouped_cluster$data_scaled_cluster5,
    data_grouped_cluster$data_scaled_cluster2,
    data_grouped_cluster$data_scaled_cluster0,
    data_grouped_cluster$data_scaled_cluster4,
    data_grouped_cluster$data_scaled_cluster1,
    data_grouped_cluster$data_scaled_cluster3
)

heatmap <- pheatmap::pheatmap(
    mat = data_heatmap,
    color = colorRampPalette(c("#3B049AFF", "white", "#FDC328FF"))(100),
    breaks = seq(from = -5, to = 5, by = 0.1),
    show_rownames = T,
    show_colnames = T,
    cluster_cols = F,
    fontsize = 5,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    cutree_rows = 5
)

protein_clusters <- cutree(heatmap$tree_row, 5) |>
    as.data.frame()

colnames(protein_clusters) <- "heatmap_clusters_continuous"

protein_clusters <- protein_clusters |>
    dplyr::mutate(
        "heatmap_clusters" = dplyr::case_when(
            heatmap_clusters_continuous == 5 ~ "cluster_A",
            heatmap_clusters_continuous == 1 ~ "cluster_B",
            heatmap_clusters_continuous == 4 ~ "cluster_C",
            heatmap_clusters_continuous == 2 ~ "cluster_D",
            heatmap_clusters_continuous == 3 ~ "cluster_E"
        )
    ) |>
    dplyr::select(heatmap_clusters)

color_annotations_heatmap_clusters <- c(
    "cluster_A" = "#30123BFF",
    "cluster_B" = "#28BBECFF",
    "cluster_C" = "#A2FC3CFF",
    "cluster_D" = "#FB8022FF",
    "cluster_E" = "#7A0403FF"
)

seurat_clusters <- as.data.frame(colnames(data_heatmap)) |>
    dplyr::rename("sample_id" = "colnames(data_heatmap)") |>
    dplyr::inner_join(data_grouping) |>
    dplyr::select(!subject) |>
    tibble::column_to_rownames("sample_id") |>
    dplyr::rename("seurat_clusters" = "cluster")

color_annotations_seurat_clusters <- c(
    "cluster_0" = "#657060FF",
    "cluster_1" = "#60CEACFF",
    "cluster_2" = "#e5c494",
    "cluster_3" = "#A11A5BFF",
    "cluster_4" = "#F05B12FF",
    "cluster_5" = "#ffff33"
)

heatmap <- pheatmap::pheatmap(
    mat = data_heatmap,
    color = colorRampPalette(c("#3B049AFF", "white", "#FDC328FF"))(100),
    breaks = seq(from = -5, to = 5, by = 0.1),
    show_rownames = T,
    show_colnames = F,
    cellheight = 5,
    legend = F,
    annotation_legend = F,
    cellwidth = 5,
    cluster_cols = F,
    treeheight_row = 10,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    fontsize = 4,
    annotation_row = protein_clusters,
    annotation_col = seurat_clusters,
    annotation_colors = list(heatmap_clusters = color_annotations_heatmap_clusters,
                             seurat_clusters = color_annotations_seurat_clusters),
    annotation_names_row = F,
    cutree_rows = 5
)

no_legend <- ggplotify::as.ggplot(heatmap)

# ggplot2::ggsave(here::here("doc/figures/figure_4/heatmap_proteomics_no_legend.png"),
#                 height = 60,
#                 width = 90,
#                 units = "mm")

legend <- heatmap <- pheatmap::pheatmap(
    mat = data_heatmap,
    color = colorRampPalette(c("#3B049AFF", "white", "#FDC328FF"))(100),
    breaks = seq(from = -5, to = 5, by = 0.1),
    show_rownames = T,
    show_colnames = F,
    cellheight = 5,
    legend = T,
    annotation_legend = T,
    cellwidth = 5,
    cluster_cols = F,
    treeheight_row = 10,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    fontsize = 4,
    annotation_row = protein_clusters,
    annotation_col = seurat_clusters,
    annotation_colors = list(heatmap_clusters = color_annotations_heatmap_clusters,
                             seurat_clusters = color_annotations_seurat_clusters),
    annotation_names_row = F,
    cutree_rows = 5
)

legend <- ggplotify::as.ggplot(legend)

# ggplot2::ggsave(here::here("doc/figures/figure_4/heatmap_proteomics_legend.png"),
#                 height = 60,
#                 width = 90,
#                 units = "mm")
