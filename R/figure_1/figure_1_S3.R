################################################################################################################################################
########################################################      Panel B   ############################################################################
################################################################################################################################################

# Load proteomics data and create Seurat object ---------------------------

data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

pca_object <- prcomp(t(data_proteomics), scale = TRUE)


# Extract %s and create scree plot ----------------------------------------

variance <- factoextra::get_eig(pca_object) |>
    as.data.frame() |>
    dplyr::select(variance.percent) |>
    dplyr::mutate(rank = 1:length(variance.percent)) |>
    dplyr::slice_head(n = 50)

# Nice Elbow plot for paper

ElbowPlot <- ggplot2::ggplot(variance,
                             ggplot2::aes(x = rank,
                                          y = variance.percent,
                                          color = rank > 6)) +
    ggplot2::annotate("rect",
                      xmin=-Inf,
                      xmax=6.5,
                      ymin=-Inf,
                      ymax=Inf,
                      alpha=0.2,
                      fill= "#045a8d") +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c("#045a8d", "#9ecae1")) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Principal Component (1-50)") +
    ggplot2::ylab("Variance per PC (%)") +
    ggplot2::ggtitle("Scree plot proteomics") +
    ggplot2::theme(
        text = ggplot2::element_text(face = "bold",size = 8, colour = "black"),
        axis.title = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = "bold")
    ) +
    ggplot2::scale_x_continuous(
        breaks = c(0,
                   10,
                   20,
                   30,
                   40,
                   50),
        labels = c(0,
                   10,
                   20,
                   30,
                   40,
                   50)
    )

ElbowPlot

ggplot2::ggsave(here::here("doc/figures/figure_1_S3/scree_plot_proteomics.png"),
                units = "mm",
                height = 60,
                width = 90)

################################################################################################################################################
########################################################      Panel D   ############################################################################
################################################################################################################################################



# Filtered and processed data
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

# Graph-based clustering using K-nearest neighbor graph  ---------------------------------

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6, seed.use = 42)

my_cols <- viridisLite::turbo(n = 5)

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 5,
    pt.size = 0.5,
    cols = ggplot2::alpha(my_cols, 0.65),
    group.by = "subject") +
    ggplot2::ggtitle("Resolution 0.4") +
    # scale_color_viridis_d(option = "magma")
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle("UMAP Proteomics - by subject") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        axis.text = ggplot2::element_text(size=6),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.text= ggplot2::element_text(size=4),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position = "right"
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_1_S3/UMAP_proteomics_subject.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )
