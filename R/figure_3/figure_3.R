################################################################################################################################################
########################################################      Panel A   ############################################################################
################################################################################################################################################

data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

metadata <- vroom::vroom(
    here::here("data/metadata_proteomics_fiber_type.csv")
) |>
    dplyr::select(!1)

pca_object <- prcomp(t(data_proteomics), scale = T)

all(colnames(data_proteomics) %in% rownames(metadata))
all(colnames(data_proteomics) == rownames(metadata))

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    tibble::add_column(
        fiber_type = metadata$fiber_type,
        digestion_batch = as.factor(metadata$digestion_batch),
        subject = as.factor(metadata$subject),
        date = metadata$date_isolation,
        length = metadata$fiber_length
    )

# Coloring PCA by GO terms ------------------------------------------------

library(org.Hs.eg.db)
meta_go2gene <- as.list(org.Hs.egGO2EG)
meta_go2goName <- as.list(GO.db::GOTERM)

# naming the goID with a human readable title. AnnotationDBI handles these 'higher' level classes.
meta_go2goName <- sapply(
    meta_go2goName,
    AnnotationDbi::Term
)

# using a named vector to not mixup our names
names(meta_go2gene) <- meta_go2goName[names(meta_go2gene)]

# translating entrez ids to hugo symbols.
meta_entrez2symbol <- AnnotationDbi::select(
    org.Hs.eg.db,
    unique(unlist(meta_go2gene)),
    "SYMBOL",
    "ENTREZID"
)

# using the power of data.table to quickly match vectors in a list
meta_entrez2symbol <- data.table::as.data.table(meta_entrez2symbol)
data.table::setkey(meta_entrez2symbol, ENTREZID)

# data.table syntax is unique, but it is also blazing fast.
meta_go2gene <- lapply(meta_go2gene, \(i){
    meta_entrez2symbol[i, unique(SYMBOL)]
})
meta_go2gene <- reshape2::melt(meta_go2gene)
colnames(meta_go2gene) <- c("SYMBOL", "GO")

vector_ribosome <- meta_go2gene[grep("cytosolic ribosome",
                                     meta_go2gene$GO),
                                "SYMBOL"] |>
    unique()

ribosome_filtered <- data_proteomics |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_ribosome) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(data_proteomics, na.rm = T)

rel_abundance_ribosomes <- (colSums(ribosome_filtered,
                                    na.rm = T) / total_intensity) * 100

data_pca_ribosomes <- as.data.frame(rel_abundance_ribosomes) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca |>
                          tibble::rownames_to_column("fiberID"))

factoextra::fviz_pca_biplot(pca_object,
                            label = "var",
                            geom = "point",
                            pointshape = 16,
                            alpha.ind = 0.65,
                            pointsize = 1,
                            col.ind = data_pca_ribosomes$rel_abundance_ribosomes,
                            col.var = "black",
                            alpha.var = 1,
                            mean.point = F,
                            addEllipses = F,
                            axes.linetype = "dashed",
                            select.var = list(name = vector_ribosome),
                            repel = TRUE
) +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Cytosolic ribosome") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Rel. \nabundance")

data_ribosomal_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2) |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ .x * 50)) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% vector_ribosome)

data_ribosomal_proteins <- data_ribosomal_proteins |>
    dplyr::mutate("x.centroid" = 0) |>
    dplyr::mutate("y.centroid" = 0) |>
    dplyr::mutate(label = dplyr::case_when(
        Genes %in% c(
            "RPL18",
            "RPS13",
            "RPS18",
            "RPL38",
            "RPL35",
            "RPL31"
        ) ~ Genes,
        TRUE ~ ""))

data_pca_ribosomes |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = rel_abundance_ribosomes
    )) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_viridis_c("Rel. \nAbundance", option = "plasma") +
    ggplot2::geom_segment(
        data = data_ribosomal_proteins,
        ggplot2::aes(x = x.centroid,
                     y = y.centroid,
                     xend = Dim.1,
                     yend = Dim.2),
        color = "black",
        arrow = grid::arrow(type = "closed",
                            length = unit(1, "mm"))
    ) +
    ggrepel::geom_label_repel(data = data_ribosomal_proteins,
                              ggplot2::aes(
                                  x = Dim.1,
                                  y = Dim.2,
                                  label = label),
                              size = 2,
                              color = "black",
                              label.padding = 0.1,
                              min.segment.length = 0.1,
                              segment.size = 0.2,
                              force = 1,
                              max.overlaps = Inf) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA by Cytosolic Ribosome GO term\n Proteomics") +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    ggplot2::theme(
        text = element_text(face="bold", colour="black", size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 3),
        legend.position = "right",
        legend.key.width = ggplot2::unit(2, "mm")
    )

ggplot2::ggsave(here::here("doc/figures/figure_3/figure_3A.pdf"),
                units = "mm",
                height = 60,
                width = 90)

################################################################################################################################################
########################################################      Panel B and C   ############################################################################
################################################################################################################################################

seurat_object <- data_proteomics |>
    Seurat::CreateSeuratObject() |>
    Seurat::FindVariableFeatures() |>
    Seurat::ScaleData() |>
    Seurat::RunPCA()

plot_RPL38 <- seurat_object |>
    Seurat::FeaturePlot("RPL38",
                        reduction = "pca",
                        cols = ggplot2::alpha(c("#cccccc", "#045a8d"), 0.5),
                        pt.size = 0.25) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 5),
        plot.title = element_text(hjust = 0.5,
                                  size = 7),
        legend.position = "none",
        # legend.key.width = ggplot2::unit(1, units = "mm"),
        legend.key.size = ggplot2::unit(2, units = "mm"),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_gradient2("LFQ intensity",
                                   low = ggplot2::alpha("darkblue", 0.7),
                                   mid = ggplot2::alpha("grey", 0.7),
                                   high = ggplot2::alpha("red", 0.7),
                                   limits = c(8, 15),
                                   midpoint = 11.5,
                                   n.breaks = 4) +
    scale_y_reverse()

# ggplot2::ggsave(here::here("doc/figures/figure_5/RPL38_proteomics.png"),
#                 units = "mm",
#                 height = 30,
#                 width = 40)

# RPL18:


plot_RPS13 <- seurat_object |>
    Seurat::FeaturePlot("RPS13",
                        reduction = "pca",
                        cols = ggplot2::alpha(c("#cccccc", "#045a8d"), 0.5),
                        pt.size = 0.25) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 5),
        plot.title = element_text(hjust = 0.5,
                                  size = 7),
        legend.position = "bottom",
        legend.title.align = 0.4,
        legend.title = ggplot2::element_text(size = 5),
        # legend.key.width = ggplot2::unit(1, units = "mm"),
        legend.key.height = ggplot2::unit(2, units = "mm"),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_gradient2("LFQ intensity",
                                   low = ggplot2::alpha("darkblue", 0.7),
                                   mid = ggplot2::alpha("grey", 0.7),
                                   high = ggplot2::alpha("red", 0.7),
                                   limits = c(8, 15),
                                   midpoint = 11.5,
                                   n.breaks = 4) +
    scale_y_reverse()

legend_plots <- ggpubr::get_legend(plot_RPS13)

plot_legend <- ggpubr::as_ggplot(legend_plots)

plot_RPS13_no_legend <- plot_RPS13 + ggplot2::theme(legend.position = "none")

plot_RPS13_no_legend + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggplot2::ggsave(here::here("doc/figures/figure_3/figure_3B.png"),
                units = "mm",
                height = 30,
                width = 40)

plot_RPL38 + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggplot2::ggsave(here::here("doc/figures/figure_3/figure_3C.png"),
                units = "mm",
                height = 30,
                width = 40)

plot_legend + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggplot2::ggsave(here::here("doc/figures/figure_3/feature_plot_legend.png"),
#                 units = "mm",
#                 height = 5,
#                 width = 40)

################################################################################################################################################
########################################################      Panel D   ############################################################################
################################################################################################################################################

metadata <- vroom::vroom(here::here("data/metadata_proteomics_seurat_clusters.csv")) |>
    dplyr::rename("sample_id" = 1)

data_proteomics <- vroom::vroom(here::here("data/data_proteomics_filtered.csv"))

# Extracting ribosomal proteins from GO term ------------------------------

library(org.Hs.eg.db)
meta_go2gene <- as.list(org.Hs.egGO2EG)
meta_go2goName <- as.list(GO.db::GOTERM)

# naming the goID with a human readable title. AnnotationDBI handles these 'higher' level classes.
meta_go2goName <- sapply(
    meta_go2goName,
    AnnotationDbi::Term
)

# using a named vector to not mixup our names
names(meta_go2gene) <- meta_go2goName[names(meta_go2gene)]

# translating entrez ids to hugo symbols.
meta_entrez2symbol <- AnnotationDbi::select(
    org.Hs.eg.db,
    unique(unlist(meta_go2gene)),
    "SYMBOL",
    "ENTREZID"
)

# using the power of data.table to quickly match vectors in a list
meta_entrez2symbol <- data.table::as.data.table(meta_entrez2symbol)
data.table::setkey(meta_entrez2symbol, ENTREZID)

# data.table syntax is unique, but it is also blazing fast.
meta_go2gene <- lapply(meta_go2gene, \(i){
    meta_entrez2symbol[i, unique(SYMBOL)]
})
meta_go2gene <- reshape2::melt(meta_go2gene)
colnames(meta_go2gene) <- c("SYMBOL", "GO")

vector_ribosome <- meta_go2gene[grep("cytosolic ribosome",
                                     meta_go2gene$GO),
                                "SYMBOL"] |>
    unique()


# Filtering dataset for ribosomal proteins --------------------------------

data_ribosomes <- data_proteomics |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    log2() |>
    # limma::normalizeQuantiles() |>
    as.data.frame() |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_ribosome) |>
    tibble::column_to_rownames("Gene.name") |>
    as.data.frame()

data_ribosomes <- vroom::vroom(here::here("data/data_ribosomes_proteomics_updated.csv")) |>
    tibble::column_to_rownames("...1")

data_heatmap <- data_ribosomes |>
    t() |>
    scale() |>
    t() |>
    as.matrix()

# metadata <- metadata |>
#     dplyr::mutate(cluster_name = dplyr::case_when(
#         seurat_clusters == "0" ~ "cluster_0",
#         seurat_clusters == "1" ~ "cluster_1",
#         seurat_clusters == "2" ~ "cluster_2",
#         seurat_clusters == "3" ~ "cluster_3",
#         seurat_clusters == "4" ~ "cluster_4",
#         seurat_clusters == "5" ~ "cluster_5",
#         TRUE ~ "0"
#     ))

# metadata <- metadata |>
#     dplyr::mutate(subject = dplyr::case_when(
#         subject == "FOR2" ~ "P1",
#         subject == "FOR4" ~ "P2",
#         subject == "FOR9" ~ "P3",
#         subject == "FOR10" ~ "P4",
#         subject == "FOR11" ~ "P5",
#         TRUE ~ "0"
#     ))

color_annotations_subject <- c(
    "P1" = "#30123BFF",
    "P2" = "#28BBECFF",
    "P3" = "#A2FC3CFF",
    "P4" = "#FB8022FF",
    "P5" = "#7A0403FF"
)

color_annotations_seurat_clusters <- c(
    "cluster_0" = "#657060FF",
    "cluster_1" = "#60CEACFF",
    "cluster_2" = "#e5c494",
    "cluster_3" = "#A11A5BFF",
    "cluster_4" = "#F05B12FF",
    "cluster_5" = "#ffff33"
)

color_annotations_fiber_type <- c(
    "Type 1" = "#440154FF",
    "Hybrid 1/2A" = "#8CB3E8",
    "Type 2A" = "#5DC863FF",
    "Hybrid 2A/2X" = "#fdc325"
)

annotations <- as.data.frame(colnames(data_heatmap)) |>
    dplyr::rename("sample_id" = 1) |>
    dplyr::inner_join(metadata |>
                          dplyr::select(c(fiber_type,
                                          subject,
                                          # cluster_name,
                                          sample_id))) |>
    tibble::column_to_rownames("sample_id")

heatmap <- pheatmap::pheatmap(
    mat = data_heatmap,
    color = colorRampPalette(c("darkblue", "white", "red"))(100),
    # color = viridisLite::viridis(n = 100),
    breaks = seq(from = -5, to = 5, by = 0.1),
    show_rownames = T,
    show_colnames = F,
    cluster_cols = T,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    annotation_col = annotations,
    annotation_colors = list(
        fiber_type = color_annotations_fiber_type,
        subject = color_annotations_subject
        # cluster_name = color_annotations_seurat_clusters
    ),
    cutree_rows = 3,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize = 5
)

ggplotify::as.ggplot(heatmap)

# Add color annotations to rows -------------------------------------------

ribosomal_clusters <- cutree(heatmap$tree_row, 3) |>
    as.data.frame() |>
    dplyr::rename("ribosomal_clusters" = 1) |>
    dplyr::mutate(ribosomal_clusters = dplyr::case_when(
        ribosomal_clusters == 3 ~ "ribosomal_cluster_1",
        ribosomal_clusters == 2 ~ "ribosomal_cluster_2",
        TRUE ~ "ribosomal_cluster_3"
    ))

ribosomal_clusters |>
    tibble::rownames_to_column("Gene.name") |>
    write.csv(here::here("data/ribosomal_clusters.csv"))

color_annotations_ribosomal_clusters <- c(
    "ribosomal_cluster_1" = "#008080",
    "ribosomal_cluster_2" = "#FB8022FF",
    "ribosomal_cluster_3" = "grey"
)

heatmap <- pheatmap::pheatmap(
    mat = data_heatmap,
    color = colorRampPalette(c("darkblue", "white", "red"))(100),
    # color = viridisLite::viridis(n = 100),
    breaks = seq(from = -5, to = 5, by = 0.1),
    show_rownames = T,
    show_colnames = F,
    cluster_cols = T,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    annotation_col = annotations,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    annotation_row = ribosomal_clusters,
    annotation_colors = list(
        fiber_type = color_annotations_fiber_type,
        subject = color_annotations_subject,
        ribosomal_clusters = color_annotations_ribosomal_clusters
        # cluster_name = color_annotations_seurat_clusters
    ),
    cutree_rows = 3,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize = 5
)

ggplotify::as.ggplot(heatmap) +
    ggplot2::theme(
        legend.position = "bottom"
    )

ggplot2::ggsave(here::here("doc/figures/figure_3/figure_3D.png"),
                height = 110,
                width = 180,
                units = "mm")
