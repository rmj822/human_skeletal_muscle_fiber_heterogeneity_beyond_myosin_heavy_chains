# Load proteomics data ----------------------------------------------------

data_proteomics <- read.csv(here::here("data/data_pca_proteomics.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

seurat_clusters <-
    vroom::vroom(here::here("data/proteomics clustering/seurat_clusters_6PC_res04.csv"))

metadata <- vroom::vroom(here::here("data/metadata_proteomics.csv"))|>
    dplyr::rename("fiberID" = 1) |>
    dplyr::inner_join(seurat_clusters) |>
    dplyr::mutate(fiber_type_seurat = dplyr::case_when(
        seurat_clusters == "3" ~ "slow",
        seurat_clusters == "1" ~ "slow",
        TRUE ~ "fast"
    )) |>
    dplyr::mutate("dataset" = dplyr::case_when(
        digestion_batch == 1 ~ "before_Xmas",
        digestion_batch == 2 ~ "before_Xmas",
        digestion_batch == 3 ~ "after_Xmas",
        digestion_batch == 4 ~ "after_Xmas",
        digestion_batch == 5 ~ "after_Xmas",
        digestion_batch == 6 ~ "after_Xmas"
    )) |>
    tibble::column_to_rownames("fiberID")


# Seurat workflow ---------------------------------------------------------

seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                              meta.data = metadata)

seurat_proteome[["RNA"]]$data <- seurat_proteome[["RNA"]]$counts

################################################################################################################################################
########################################################      PCA   ############################################################################
################################################################################################################################################

# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                                selection.method = "vst")

# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)

#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,
                                  features = Seurat::VariableFeatures(object = seurat_proteome))

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))

################################################################################################################################################
################################################      DIMENSIONALITY REDUCTION   ##############################################################
################################################################################################################################################

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6)

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 1,
    pt.size = 0.05,
    group.by = "fiber_type_seurat"
)

seurat_proteome[["umap"]]@cell.embeddings |>
    as.data.frame() |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = metadata$fiber_type_seurat)
    ) +
    ggplot2::geom_point(size = 0.025,
                        alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(
        "#5DC863FF",
        "#440154FF"
    )) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        # text = ggplot2::element_blank(),
        legend.position = "none",
        # axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

# Make MYH coloring -------------------------------------------------------

perc_MYHs <- vroom::vroom(here::here("data/perc_MYH_proteomics.csv")) |>
    dplyr::select(!1)

coloring <- perc_MYHs |>
    dplyr::select(fiber_ID, MYHs, values) |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    dplyr::filter(fiber_ID %in% colnames(data_proteomics)) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::mutate(
        merged_fast_MYHs = MYH2 - MYH1,
    MYHs_fast = MYH2 + MYH1) |>
    as.data.frame()


test <- data_proteomics |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% c("MYH7", "MYH2", "MYH1")) |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    dplyr::mutate(
        MYH2_subs = MYH7 - MYH2,
        MYH1_subs = MYH7 - MYH1,
        merged_fast_MYH = MYH2_subs - MYH1_subs,
        merged_fast_MYH_test = (MYH2/(MYH7 + MYH2 + MYH1) * 100) - (MYH1/(MYH7 + MYH2 + MYH1)*100)
    ) |>
    # scale() |>
    as.data.frame()

seurat_proteome[["umap"]]@cell.embeddings |>
    as.data.frame() |>
    ggplot2::ggplot(
        ggplot2::aes(x = UMAP_1,
                     y = UMAP_2,
                     color = coloring$merged_fast_MYHs)
    ) +
    ggplot2::geom_point(ggplot2::aes(alpha = coloring$MYHs_fast),
                        size = 1) +
    ggplot2::scale_alpha(range = c(0.2, 1),
                guide = "none") +
    ggplot2::scale_color_gradient2(
        "MYH2 % \n     - \nMYH1 %",
        low = "#D2631C",
        mid = "lightgrey",
        high = "#23CE6B",
        midpoint = 0,
        limits = c(-100, 100),
        breaks = c(-100, -50, 0, 50, 100)
    ) +
    ggplot2::ggtitle("UMAP proteomics - by fast MYH blended expression") +
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

# ggplot2::ggsave(
#     here::here("doc/figures/figure_1/fast_MYHs_blended_proteomics_scaled.png"),
#     units = "mm",
#     height = 60,
#     width = 90
# )


# Coloring using RGB ------------------------------------------------------

perc_MYHs <- perc_MYHs |>
    dplyr::select(c(fiber_ID, MYHs, values)) |>
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
    here::here("doc/figures/figure_1/all_MYHs_blended_proteomics_scaled.pdf"),
    units = "mm",
    height = 60,
    width = 80
)
