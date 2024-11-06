source(here::here("R/figure_1/MYH_curves.r"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

data_proteomics <-read.csv(here::here("data/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "X") |>
    tibble::column_to_rownames("Gene.name")

metadata <-read.csv(here::here("data/metadata_proteomics.csv"))|>
    dplyr::rename("Sample" = "X") |>
    tibble::column_to_rownames("Sample")

pca_object <- prcomp(t(data_proteomics), scale = T)

all(colnames(data_proteomics) %in% rownames(metadata))
all(colnames(data_proteomics) == rownames(metadata))

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    tibble::add_column(
        fiber_type = metadata$fiber_type,
        fiberID = rownames(metadata),
        digestion_batch = as.factor(metadata$digestion_batch),
        subject = as.factor(metadata$subject),
        date = metadata$date_isolation,
        length = metadata$fiber_length
    )

data_pca$fiber_type <- factor(data_pca$fiber_type,
                             levels = c("Hybrid 1/2A",
                                        "Hybrid 2A/2X",
                                        "Type 2A",
                                        "Type 1"))

pca_proteomics <- data_pca |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = fiber_type,
            names = fiberID
        )
    ) +
    ggplot2::geom_point(
        size = 1,
        alpha = 0.65,
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA Proteomics") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        )
    ) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           vjust = 1.5)
    ) +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        axis.text = ggplot2::element_text(size=8),
        ) +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    ggplot2::scale_color_manual(
        "Fiber type",
        values = c("#8CB3E8",
                   "#fdc325",
                   "#5DC863FF",
                   "#440154FF")
    )

# ggplot2::ggsave(here::here("doc/figures/figure_3/PCA_proteomics_fiber_type.png"),
#                 units = "mm",
#                 height = 60,
#                 width = 90)


# PCA by subject proteomics -----------------------------------------------
data_pca <- data_pca %>%
    dplyr::mutate(subject = dplyr::case_when(
        subject == "FOR2" ~ "P1",
        subject == "FOR4" ~ "P2",
        subject == "FOR9" ~ "P3",
        subject == "FOR10" ~ "P4",
        subject == "FOR11" ~ "P5",
        TRUE ~ "unknown"
    ))

data_pca$subject <- factor(data_pca$subject, levels = c("P1", "P2", "P3", "P4", "P5"))

PCA_subject <- Seurat::DimPlot(seurat_by_subject,
                               label = FALSE,
                               reduction = "pca",
                               pt.size = 0.5,
                               group.by = "subject") +
    ggtitle("PCA Transcriptomics - by subject") +
    xlab("PC1 (11.1%)") +
    ylab("PC2 (3.5%)") +
    guides(color = guide_legend(override.aes = list(size=1), ncol=2)) +
     theme(
        text = element_text(face="bold", colour="black", size=8),
        axis.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = "right",
        legend.text=element_text(size=4)
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    )

data_pca |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = subject,
            names = fiberID
        )
    ) +
    ggplot2::geom_point(
        size = 0.5,
        alpha = 1
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("PCA Proteomics by subject") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        )
    ) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::theme(
        legend.position = "right",
        text = ggplot2::element_text(face="bold", colour="black", size=8),
        axis.text = ggplot2::element_text(size=8),
        legend.text=element_text(size=4)
    ) +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    guides(color = guide_legend(override.aes = list(size=1), ncol=1)) +
    scale_color_manual(values = c("#786A7F",
                                  "#8FCAE8",
                                  "#C8D685",
                                  "#EAAE7C",
                                  "#9F6461"),
                       name = "")


 ggplot2::ggsave(here::here("doc/figures/figure_2/PCA_proteomics_subject.png"),
                 units = "mm",
                 height = 60,
                 width = 90)


# Combined PCA plots ------------------------------------------------------

load(here::here("data/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata"))

factor(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids,
       levels = c("Hybrid 1/2A", "Type 1", "Type 2A", "Type 2X","Hybrid 2A/2X"))

filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids <-
    factor(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids,
           levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))

Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "fiber_type_MYH_hybrids"

my_cols <- c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325", "#D2631C")

pca_transcriptomics <- Seurat::DimPlot(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                       label = FALSE,
                                       reduction = "pca",
                                       cols = ggplot2::alpha(my_cols, 0.65),
                                       pt.size = 1,
                                       order = c("Hybrid 1/2A",
                                                 "Type 2X",
                                                 "Hybrid 2A/2X",
                                                 "Type 2A",
                                                 "Type 1")) +
    guides(color = guide_legend(override.aes = list(size=1))) +
    theme_minimal() +
    ggtitle("PCA Transcriptomics") +
    xlab("PC1 (11.1%)") +
    ylab("PC2 (3.5%)") +
    theme(
        text = element_text(face="bold", colour="black", size=8),
        axis.text = element_text(size=8),
        # strip.text = element_text(colour = "white"),
        # strip.background = element_rect(fill="black"),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 3),
        panel.background = element_rect(fill='transparent', color = "transparent"), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg

    ) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 8,
            face = "bold"
        )) +
    scale_color_manual(values = ggplot2::alpha(my_cols, 0.65),
                       limits = c("Type 1",
                                  "Hybrid 1/2A",
                                  "Type 2A",
                                  "Hybrid 2A/2X",
                                  "Type 2X")) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    )


PCA_combined <- pca_transcriptomics + pca_proteomics + plot_layout(guides = "collect")

# ggsave(PCA_combined,
#        filename = "~/single_fiber_heterogeneity/doc/figures/figure_3/PCA_combined.png",
#        width = 195,
#        height = 60,
#        units="mm")


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
    dplyr::inner_join(data_pca)

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

ggplot2::ggsave(here::here("doc/figures/figure_3/PCA_proteomics_ribosomes.pdf"),
                units = "mm",
                height = 60,
                width = 90)


# PCA plot by costamere -----------------------------------------------

vector_costamere <- meta_go2gene[grep("costamere",
                                     meta_go2gene$GO),
                                "SYMBOL"] |>
    unique()

costamere_filtered <- data_proteomics |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_costamere) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(data_proteomics, na.rm = T)

rel_abundance_costamere <- (colSums(costamere_filtered,
                                    na.rm = T) / total_intensity) * 100

data_pca_costamere <- as.data.frame(rel_abundance_costamere) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca)

factoextra::fviz_pca_biplot(pca_object,
                            label = "var",
                            geom = "point",
                            pointshape = 16,
                            alpha.ind = 0.65,
                            pointsize = 1,
                            col.ind = data_pca_costamere$rel_abundance_costamere,
                            col.var = "black",
                            alpha.var = 1,
                            mean.point = F,
                            addEllipses = F,
                            axes.linetype = "dashed",
                            select.var = list(name = vector_costamere),
                            repel = TRUE
) +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Cell costamere") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8, face = "bold")) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Rel. \nabundance")

data_costamere_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2) |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ .x * 50)) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% vector_costamere)

data_costamere_proteins <- data_costamere_proteins |>
    dplyr::mutate("x.centroid" = 0) |>
    dplyr::mutate("y.centroid" = 0) |>
    dplyr::mutate(label = dplyr::case_when(
        Genes %in% c(
            "SVIL",
            "DMD",
            "SYNM",
            "VCL",
            "PLEC",
            "ANK2"
        ) ~ Genes,
        TRUE ~ ""))

data_pca_ribosomes |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = rel_abundance_costamere
    )) +
    ggplot2::geom_point(alpha = 0.5,
                        size = 1) +
    ggplot2::scale_color_viridis_c("Rel. Abundance", option = "plasma") +
    ggplot2::geom_segment(
        data = data_costamere_proteins,
        ggplot2::aes(x = x.centroid,
                     y = y.centroid,
                     xend = Dim.1,
                     yend = Dim.2),
        color = "black",
        arrow = grid::arrow(type = "closed",
                            length = unit(1, "mm"))
    ) +
    ggrepel::geom_label_repel(data = data_costamere_proteins,
                              ggplot2::aes(
                                  x = Dim.1,
                                  y = Dim.2,
                                  label = label),
                              size = 2,
                              color = "black",
                              label.padding = 0.1,
                              min.segment.length = 0.1,
                              segment.size = 0.2,
                              force = 20,
                              max.overlaps = Inf) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA by Costamere \nGO term - Proteomics") +
    ggplot2::xlab("PC1 (12.4%)") +
    ggplot2::ylab("PC2 (10.6%)") +
    ggplot2::theme(
        text = element_text(face="bold", colour="black", size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 3),
        legend.position = "bottom",
        legend.key.height = ggplot2::unit(2, "mm")
    )

# ggplot2::ggsave(here::here("doc/figures/figure_4/PCA_proteomics_costamere.png"),
#                 units = "mm",
#                 height = 65,
#                 width = 60)

# Feature Plots -----------------------------------------------------------

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
# ggplot2::ggsave(here::here("doc/figures/figure_5/RPS13_proteomics.png"),
#                 units = "mm",
#                 height = 30,
#                 width = 40)

plot_RPL38 + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggplot2::ggsave(here::here("doc/figures/figure_5/RPL38_proteomics.png"),
#                 units = "mm",
#                 height = 30,
#                 width = 40)

plot_legend + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggplot2::ggsave(here::here("doc/figures/figure_5/feature_plot_legend.png"),
#                 units = "mm",
#                 height = 5,
#                 width = 40)

# Plotting them merged ----------------------------------------------------

# ggpubr::ggarrange(list(plot_RPL38,
#                   plot_RPL18_no_legend,
#                   plot_legend),
#                   nrow = 2,
#                   ncol = 1,
#                   legend = "none"
#                   )
#
# plot_RPL38 + plot_RPL18 + plot_layout(ncol = 1, nrow = 2, guides = "collect")
#
#
# cowplot::plot_grid(
#     list(plot_RPL38,
#          plot_RPL18_no_legend,
#          legend_plots),
#     nrow = 3,
#     ncol = 1,
#     rel_heights = c(1, 1, 0.2)
# )

# Changing log scale to raw intensities -----------------------------------

seurat_object@assays$RNA@counts <- seurat_object@assays$RNA@counts |>
    as.data.frame() |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ 2^.x)) |>
    as.matrix()

new_seurat_object <- Seurat::CreateSeuratObject(seurat_object@assays$RNA@counts)

new_seurat_object@assays$RNA@counts

seurat_object@assays$RNA@counts <- new_seurat_object@assays$RNA@counts

Seurat::FeaturePlot(seurat_object,
                    "RPL38",
                    reduction = "pca",
                    slot = "counts",
                    cols = ggplot2::alpha(c("#cccccc", "#045a8d"), 0.5),
                    pt.size = 0.25) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 5),
        plot.title = element_text(hjust = 0.5,
                                  size = 7),
        # legend.key.width = ggplot2::unit(1, units = "mm"),
        legend.key.size = ggplot2::unit(2, units = "mm"),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
    ) +
    scale_y_reverse()

# ggplot2::ggsave(here::here("doc/figures/figure_3/RPL38_proteomics_raw_color.png"),
#                 units = "mm",
#                 height = 30,
#                 width = 40)

Seurat::FeaturePlot(seurat_object,
                    "RPL18",
                    reduction = "pca",
                    slot = "counts",
                    cols = ggplot2::alpha(c("#cccccc", "#045a8d"), 0.5),
                    pt.size = 0.25) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 5),
        plot.title = element_text(hjust = 0.5,
                                  size = 7),
        # legend.key.width = ggplot2::unit(1, units = "mm"),
        legend.key.size = ggplot2::unit(2, units = "mm"),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
    ) +
    scale_y_reverse()

# ggplot2::ggsave(here::here("doc/figures/figure_3/RPL18_proteomics_raw_color.png"),
#                 units = "mm",
#                 height = 30,
#                 width = 40)


# PCA_drivers -------------------------------------------------------------

PC_genes <- pca_object$rotation |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    dplyr::rename("PC_1" = PC1,
                  "PC_2" = PC2) |>
    tibble::rownames_to_column("Gene")

# Data wrangling
Top_PC1_pos <- PC_genes |>
    dplyr::arrange(desc(PC_1)) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_1) |>
    dplyr::arrange(PC_1)

Top_PC1_neg <- PC_genes  |>
    dplyr::arrange(PC_1) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_1)

Top_PC2_pos <- PC_genes |>
    dplyr::arrange(desc(PC_2)) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_2)  |>
    dplyr::arrange(PC_2)

Top_PC2_neg <- PC_genes |>
    dplyr::arrange(PC_2) |>
    dplyr::slice(1:8) |>
    dplyr::select(Gene, PC_2)

PC_df <- data.frame(
    gene = c(Top_PC1_neg$Gene,
             Top_PC1_pos$Gene,
             Top_PC2_neg$Gene,
             Top_PC2_pos$Gene),
    PC_score = c(Top_PC1_neg$PC_1,
                 Top_PC1_pos$PC_1,
                 Top_PC2_neg$PC_2,
                 Top_PC2_pos$PC_2),
    PC = c(rep("PC1", 16),
           (rep("PC2", 16))),
    order = c(1:32)
)

# Create plot
PC_drivers_plot <- ggplot2::ggplot(PC_df,
                                   ggplot2::aes(x = PC_score,
                                                y = forcats::fct_reorder(gene,
                                                                         order,
                                                                         .desc = TRUE))) +
    ggplot2::geom_col(fill = "#045a8d") +
    ggplot2::facet_grid(~ PC, scales="free") +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = 0, colour="black") +
    ggplot2::ggtitle("PC drivers proteomics") +
    ggplot2::ylab("Genes") +
    ggplot2::ylab("PC score") +
    ggplot2::theme(
        text = ggplot2::element_text(face = "bold",size = 7, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face="bold"),
        legend.position = "none",
        strip.background = ggplot2::element_rect(fill= c("#9ecae1")),
        strip.text.x = ggplot2::element_text(size = 6, face = "bold", hjust=0.5),
        axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
        axis.title.y = ggplot2::element_blank()
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
    )

# ggsave(PC_drivers_plot,
#        filename = here::here("doc/figures/figure_4/PC_drivers_proteomics.png"),
#        width = 60,
#        height = 90,
#        units="mm")


