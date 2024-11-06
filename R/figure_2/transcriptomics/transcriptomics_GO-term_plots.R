
library(factoextra)
library(org.Hs.eg.db)

################################################################################################################################################
###########################################      GET NORMALIZED COUNTS      ####################################################################
################################################################################################################################################

# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object (in data-raw folder) ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

# Extract SCT scaled data
transcriptomics_scaledata <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "scale.data")

# Extract SCT normalized but not log-transformed counts
transcriptomics_counts <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")


################################################################################################################################################
###########################################      PERFORM PCA WITH PRCOMP      ##################################################################
################################################################################################################################################

pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2)

data_pca$fiberID <- rownames(data_pca)

################################################################################################################################################
################################################      GO TERM PREPARATION     ##################################################################
################################################################################################################################################

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


################################################################################################################################################
################################################      BY GO CELL JUNCTION TERM      ################################################################
################################################################################################################################################

vector_junction <- meta_go2gene[grep("cell-cell adhesion",
                                     meta_go2gene$GO),
                                "SYMBOL"] |>
    unique()

junction_filtered <- transcriptomics_counts |>
    as.data.frame() %>%
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_junction) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(transcriptomics_counts, na.rm = T)

rel_abundance_junction <- (colSums(junction_filtered,
                                    na.rm = T) / total_intensity) * 100

data_pca_junction <- as.data.frame(rel_abundance_junction) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca)

data_junction_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2) |>
    dplyr::mutate(dplyr::across(.cols = everything(),
                                ~ .x * 60)) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% vector_junction)

data_junction_proteins <- data_junction_proteins |>
    dplyr::mutate("x.centroid" = 0) |>
    dplyr::mutate("y.centroid" = 0) |>
    dplyr::mutate(label = dplyr::case_when(
        Genes %in% c(
            "LRRC7",
            "MAGI1",
            "ADGRL3",
            "TENM2",
            "CDH19"
        ) ~ Genes,
        TRUE ~ ""))

custom_junction_plot_with_legend <- data_pca_junction |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC2,
        color = rel_abundance_junction
    )) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::scale_color_viridis_c("Rel. Abundance", option = "plasma") +
    ggplot2::geom_segment(
        data = data_junction_proteins,
        ggplot2::aes(x = x.centroid,
                     y = y.centroid,
                     xend = Dim.1,
                     yend = Dim.2),
        color = "black",
        arrow = grid::arrow(type = "closed",
                            length = unit(1, "mm"))
    ) +
    ggrepel::geom_label_repel(data = data_junction_proteins,
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
    ggplot2::ggtitle("PCA colored by Cell-cell adhesion \n GO term - Transcriptomics") +
    ggplot2::xlab("PC1 (11.1%)") +
    ggplot2::ylab("PC2 (3.5%)") +
    ggplot2::theme(
        text = element_text(face="bold", colour="black", size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 3),
        legend.position = "bottom",
        legend.key.height = ggplot2::unit(2, "mm")
    ) +
    scale_y_continuous(trans = "reverse",
                       breaks=c(-40, -20, 0, 20),
                       labels=c("40", "20", "0", "-20")
    ) +
    scale_x_continuous(trans = "reverse",
                       breaks=c(-50, 0, 50, 100),
                       labels=c("50", "0", "-50", "-100")
    )


ggsave(custom_junction_plot_with_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/PCA_transcriptomics_junction.png",
       width = 60,
       height = 65,
       units="mm")

################################################################################################################################################
################################################      FEATURE PLOTS RIBO     ###################################################################
################################################################################################################################################

# RPL27 -------------------------------------------------------------------
featureplot_RPL27 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                 features = c("RPL27"),
                                 order = TRUE,
                                 min.cutoff = 'q10',
                                 label = F,
                                 repel = TRUE,
                                 reduction = "pca",
                                 cols = c("white", "#354D3C"),
                                 ncol = 1,
                                 pt.size = 0.25) +
    xlab("PC1") +
    ylab("PC2") +
    theme(
        text = element_text(face="bold", colour="black", size=5),
        plot.title = element_text(hjust = 0.5, size=7),
        #legend.key.size = unit(0.3, 'cm')
        legend.position = "none",
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
    ) +
    scale_y_reverse()

ggsave(featureplot_RPL27, filename = "~/single_fiber_heterogeneity/doc/figures/figure_3/Featureplot_RPL27_transcriptomics.png", width = 45, height = 30, units="mm")

# RPL31 -------------------------------------------------------------------
featureplot_RPL31 <- FeaturePlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                 features = c("RPL31"),
                                 order = TRUE,
                                 min.cutoff = 'q10',
                                 label = F,
                                 repel = TRUE,
                                 reduction = "pca",
                                 cols = c("white", "#354D3C"),
                                 ncol = 1,
                                 pt.size = 0.25) +
    xlab("PC1") +
    ylab("PC2") +
    theme(
        text = element_text(face="bold", colour="black", size=5),
        plot.title = element_text(hjust = 0.5, size=7),
        #legend.key.size = unit(0.3, 'cm')
        legend.position = "none",
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
    ) +
    scale_y_reverse()

ggsave(featureplot_RPL31, filename = "~/single_fiber_heterogeneity/doc/figures/figure_3/Featureplot_RPL31_transcriptomics.png", width = 45, height = 30, units="mm")
