
library(Seurat)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(ggh4x)

################################################################################################################################################
################################################       PROTEOMICS DATA WRANGLING       ####@####################################################
################################################################################################################################################

data_proteomics <-read.csv(here::here("data/figure_2/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "X") |>
    tibble::column_to_rownames("Gene.name")

metadata <-read.csv(here::here("data/figure_2/metadata_proteomics.csv"))|>
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


################################################################################################################################################
################################################       TRANSCRIPTOMICS DATA WRANGLING       ####@####################################################
################################################################################################################################################

load(here::here("data/figure_2/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata"))

factor(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids,
       levels = c("Hybrid 1/2A", "Type 1", "Type 2A", "Type 2X","Hybrid 2A/2X"))

filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids <-
    factor(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data$fiber_type_MYH_hybrids,
           levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))

Idents(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest) <- "fiber_type_MYH_hybrids"

my_cols <- c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325", "#D2631C")


################################################################################################################################################
################################################       FIGURE 2A       ####@####################################################
################################################################################################################################################

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


pca_transcriptomics + pca_proteomics + plot_layout(guides = "collect")

ggplot2::ggsave(here::here("doc/figures/figure_2/figure_2A.png"),
                units = "mm",
                height = 195,
                width = 60)

################################################################################################################################################
###############################################      FIGURE 2B - FIGURE CREATION      ##########################################################
################################################################################################################################################

Tx_PC1_pos <- read.csv(here::here("data/figure_2/PC1 pos_T.csv"))
Tx_PC1_neg <- read.csv(here::here("data/figure_2/PC1 neg_T.csv"))
Tx_PC2_pos <- read.csv(here::here("data/figure_2/PC2 pos_T.csv"))
Tx_PC2_neg <- read.csv(here::here("data/figure_2/PC2 neg_T.csv"))

Px_PC1_pos <- read.csv(here::here("data/figure_2/results_PC1_positive_P.csv"))
Px_PC1_neg <- read.csv(here::here("data/figure_2/results_PC1_negative_P.csv"))
Px_PC2_pos <- read.csv(here::here("data/figure_2/results_PC2_positive_P.csv"))
Px_PC2_neg <- read.csv(here::here("data/figure_2/results_PC2_negative_P.csv"))

# Calculate foldEnrich for proteomics
Px_PC1_pos <- mutate(Px_PC1_pos, foldEnrich =
                         (as.numeric(sub("/\\d+", "", Px_PC1_pos$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC1_pos$GeneRatio))) /
                         (as.numeric(sub("/\\d+", "", Px_PC1_pos$BgRatio)) / as.numeric(sub(".*/", "", Px_PC1_pos$BgRatio))))

Px_PC1_neg <- mutate(Px_PC1_neg, foldEnrich =
                         (as.numeric(sub("/\\d+", "", Px_PC1_neg$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC1_neg$GeneRatio))) /
                         (as.numeric(sub("/\\d+", "", Px_PC1_neg$BgRatio)) / as.numeric(sub(".*/", "", Px_PC1_neg$BgRatio))))

Px_PC2_pos <- mutate(Px_PC2_pos, foldEnrich =
                         (as.numeric(sub("/\\d+", "", Px_PC2_pos$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC2_pos$GeneRatio))) /
                         (as.numeric(sub("/\\d+", "", Px_PC2_pos$BgRatio)) / as.numeric(sub(".*/", "", Px_PC2_pos$BgRatio))))

Px_PC2_neg <- mutate(Px_PC2_neg, foldEnrich =
                         (as.numeric(sub("/\\d+", "", Px_PC2_neg$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC2_neg$GeneRatio))) /
                         (as.numeric(sub("/\\d+", "", Px_PC2_neg$BgRatio)) / as.numeric(sub(".*/", "", Px_PC2_neg$BgRatio))))


# Select categories for PC1 positive - Tx and Px
Px_PC1_pos_plot <- Px_PC1_pos %>%
    dplyr::filter(
        Description == "cytosolic ribosome" |
            Description == "cytoplasmic translation" |
            Description == "structural constituent of ribosome" |
            Description == "mitochondrial respirasome"
    )

Tx_PC1_pos_plot <- Tx_PC1_pos %>%
    dplyr::filter(
        Description == "cytosolic ribosome" |
            Description == "cytoplasmic translation" |
            Description == "structural constituent of ribosome" |
            Description == "mitochondrial respirasome"
    )

# Select categories for PC1 negative - Tx and Px
Px_PC1_neg_plot <- Px_PC1_neg %>%
    dplyr::filter(
        Description == "cell junction" |
            Description == "costamere" |
            Description == "enzyme binding"
    )

Tx_PC1_neg_plot <- Tx_PC1_neg %>%
    dplyr::filter(
        Description == "GTPase regulator activity" |
            Description == "cell-cell adhesion" |
            Description == "cell-cell junction"
    )

# Select categories for PC2 positive - Tx and Px
Px_PC2_pos_plot <- Px_PC2_pos %>%
    dplyr::filter(
        Description == "fatty-acyl-CoA binding" |
            Description == "oxidoreductase activity" |
            Description == "contractile fiber"
    )

Tx_PC2_pos_plot <- Tx_PC2_neg %>%
    dplyr::filter(
        Description == "fatty-acyl-CoA binding" |
            Description == "oxidoreductase activity" |
            Description == "contractile fiber"
    )

# Select categories for PC2 negative - Tx and Px
Px_PC2_neg_plot <- Px_PC2_neg %>%
    dplyr::filter(
        Description == "carbohydrate binding" |
            Description == "nucleotide phosphorylation" |
            Description == "contractile fiber"
    )

Tx_PC2_neg_plot <- Tx_PC2_pos %>%
    dplyr::filter(
        Description == "carbohydrate binding" |
            Description == "nucleotide phosphorylation" |
            Description == "contractile fiber"
    )

# Add required columns for plotting
Px_PC1_pos_plot$PCDir <- rep("PC1 - Pos", nrow(Px_PC1_pos_plot))
Px_PC1_pos_plot$graphdir <- rep(3, nrow(Px_PC1_pos_plot))
Px_PC1_pos_plot$data <- rep("Px", nrow(Px_PC1_pos_plot))

Tx_PC1_pos_plot$PCDir <- rep("PC1 - Pos", nrow(Tx_PC1_pos_plot))
Tx_PC1_pos_plot$graphdir <- rep(4, nrow(Tx_PC1_pos_plot))
Tx_PC1_pos_plot$data <- rep("Tx", nrow(Tx_PC1_pos_plot))

Px_PC1_neg_plot$PCDir <- rep("PC1 - Neg", nrow(Px_PC1_neg_plot))
Px_PC1_neg_plot$graphdir <- rep(1, nrow(Px_PC1_neg_plot))
Px_PC1_neg_plot$data <- rep("Px", nrow(Px_PC1_neg_plot))

Tx_PC1_neg_plot$PCDir <- rep("PC1 - Neg", nrow(Tx_PC1_neg_plot))
Tx_PC1_neg_plot$graphdir <- rep(2, nrow(Tx_PC1_neg_plot))
Tx_PC1_neg_plot$data <- rep("Tx", nrow(Tx_PC1_neg_plot))

Px_PC2_pos_plot$PCDir <- rep("PC2 - Pos", nrow(Px_PC2_pos_plot))
Px_PC2_pos_plot$graphdir <- rep(7, nrow(Px_PC2_pos_plot))
Px_PC2_pos_plot$data <- rep("Px", nrow(Px_PC2_pos_plot))

Tx_PC2_pos_plot$PCDir <- rep("PC2 - Pos", nrow(Tx_PC2_pos_plot))
Tx_PC2_pos_plot$graphdir <- rep(8, nrow(Tx_PC2_pos_plot))
Tx_PC2_pos_plot$data <- rep("Tx", nrow(Tx_PC2_pos_plot))

Px_PC2_neg_plot$PCDir <- rep("PC2 - Neg", nrow(Px_PC2_neg_plot))
Px_PC2_neg_plot$graphdir <- rep(5, nrow(Px_PC2_neg_plot))
Px_PC2_neg_plot$data <- rep("Px", nrow(Px_PC2_neg_plot))

Tx_PC2_neg_plot$PCDir <- rep("PC2 - Neg", nrow(Tx_PC2_neg_plot))
Tx_PC2_neg_plot$graphdir <- rep(6, nrow(Tx_PC2_neg_plot))
Tx_PC2_neg_plot$data <- rep("Tx", nrow(Tx_PC2_neg_plot))

# Combine all into one dataframe
PC_plot <- bind_rows(Tx_PC2_pos_plot, Px_PC2_pos_plot, Tx_PC2_neg_plot, Px_PC2_neg_plot, Tx_PC1_pos_plot, Px_PC1_pos_plot, Tx_PC1_neg_plot, Px_PC1_neg_plot)
PC_plot$graph <- rep(1, nrow(PC_plot))
PC_plot$PCDir = factor(PC_plot$PCDir, levels=c("PC2 - Pos", "PC2 - Neg", "PC1 - Pos", "PC1 - Neg"))
PC_plot$data = factor(PC_plot$data, levels=c("Tx", "Px"))

# Plot
ggplot(PC_plot, aes(graph, fct_reorder(Description, graphdir))) +
    geom_point(aes(color= data, alpha=foldEnrich, size = -log10(p.adjust))) +
    scale_color_manual(values = c("#354D3C", "#045a8d", "#354D3C", "#045a8d", "#354D3C","#045a8d", "#354D3C", "#045a8d")) +
    scale_alpha_continuous(
        name = "Fold enrichment",
        limits = c(2.0004, 12.574),
        breaks = c(2.0004, 12.574),
        range = c(0.2, 1)) +
    scale_size_continuous(
        name = "-log10(p.adjust)",
        limits = c(1.430, 33.73),
        breaks = c(1.430, 33.73),
        range=c(1, 3.5)) +
    ggh4x::facet_nested_wrap(~PCDir + data , nrow = 1, ncol=8) +
    theme_classic() +
    scale_x_continuous(limits = c(0.9,1.1), expand = c(0, 0)) +
    theme(
        text = ggplot2::element_text(face = "bold",size = 8, colour = "black"),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x= element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x =element_blank(),
        #axis.text.y = element_text(size=8),
        # plot.title = element_text(hjust = 0.5, face="bold", vjust=-7),
        plot.title = element_text(hjust = 1, face="bold"),
        legend.position = "none",
        strip.background = element_rect(colour=NA, fill="white"),
        strip.text = element_text(colour = "white")
    )

ggsave(here::here("doc/figures/figure_2/figure_2B.png"),
       width = 80,
       height = 100,
       units="mm")


################################################################################################################################################
###############################################      FIGURE 2C-D - DATA PREP      ##########################################################
################################################################################################################################################

# Coloring PCA by GO terms ------------------------------------------------
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
#######################################################      FIGURE 2C      ####################################################################
################################################################################################################################################

# Extract SCT scaled data
transcriptomics_scaledata <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "scale.data")

# Extract SCT normalized but not log-transformed counts
transcriptomics_counts <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")

pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2)

data_pca$fiberID <- rownames(data_pca)

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

data_pca_junction |>
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

ggsave(here::here("doc/figures/figure_2/figure_2C.png"),
       width = 60,
       height = 65,
       units="mm")


################################################################################################################################################
###############################################      FIGURE 2D   ##########################################################
################################################################################################################################################

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
    dplyr::mutate(Dim.1 = -Dim.1) %>%
    dplyr::mutate(Dim.2 = -Dim.2) %>%
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

data_pca_costamere <- data_pca_costamere %>%
    dplyr::mutate(PC1 = -PC1) %>%
    dplyr::mutate(PC2 = -PC2)

data_pca_costamere |>
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

ggsave(here::here("doc/figures/figure_2/figure_2D.png"),
       width = 60,
       height = 65,
       units="mm")

################################################################################################################################################
###############################################      FIGURE 2E   ##########################################################
################################################################################################################################################


feature_plot_SLIT3 <- Seurat::FeaturePlot(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                          reduction = "umap",
                                          features = c("SLIT3"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin = ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    scale_x_reverse()

feature_plot_CHCHD10 <- Seurat::FeaturePlot(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                            features = c("CHCHD10"),
                                            pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin= ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    scale_x_reverse()


feature_plot_CTDNEP1 <- Seurat::FeaturePlot(filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest,
                                            features = c("CTDNEP1"),
                                            pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin= ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    scale_x_reverse()


final_plot <- ggpubr::ggarrange(feature_plot_CHCHD10,
                                feature_plot_SLIT3,
                                feature_plot_CTDNEP1,
                                ncol = 3,
                                nrow = 1) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggpubr::annotate_figure(final_plot, top = ggpubr::text_grob("Transcriptomics",
                                                                          color = "black", face = "bold", size = 8))

ggsave(here::here("doc/figures/figure_2/figure_2E.png"),
       width = 120,
       height = 35,
       units="mm")


################################################################################################################################################
###############################################      FIGURE 2F   ##########################################################
################################################################################################################################################

data_proteomics <- read.csv(here::here("data/data_pca_proteomics.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

# Proteome metadata
metadata_proteomics <- read.csv(here::here("data/metadata_proteomics.csv"))

metadata_proteomics <- metadata_proteomics |>  dplyr::rename("fiberID" = X)
rownames(metadata_proteomics) <- metadata_proteomics$fiberID

################################################################################################################################################
################################################       CREATE SEURAT OBJECT      ###############################################################
################################################################################################################################################

seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                              meta.data = metadata_proteomics)


################################################################################################################################################
########################################################      PCA   ############################################################################
################################################################################################################################################

# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                                selection.method = "vst")

# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)

#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,  features = Seurat::VariableFeatures(object = seurat_proteome))

# Explore PCA plots
PCA_PC1_2 <- Seurat::DimPlot(seurat_proteome, reduction = "pca", dims = c(1,2), group.by = "fiber_type") # No separation of fiber type by PC1, good separation by PC2
PCA_PC1_3 <- Seurat::DimPlot(seurat_proteome, reduction = "pca", dims = c(1,3), group.by = "fiber_type") # No separation at all by PC3

# ggsave(PCA_PC1_2, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/PCA_proteome_PC1_PC2.png",  width = 128, height = 90, units="mm")
# ggsave(PCA_PC1_3, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/PCA_proteome_PC1_PC3.png",  width = 128, height = 90, units="mm")


################################################################################################################################################
################################################      Identify significant PCS    ##############################################################
################################################################################################################################################

# Explore heatmap of PCs (try to find PC where heatmap starts to look fuzzy, not so distinct between groups) --------

PCA_heatmaps <- Seurat::DimHeatmap(seurat_proteome,
                                   dims = 1:15,
                                   cells = 500,  # number of cells with most negative or positive PCA scores to use for plotting
                                   balanced = TRUE) # Still decent clustering up until PC 15

# ggsave(PCA_heatmaps, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/PCA_heatmaps_proteome.png",  width = 128, height = 90, units="mm")


# Elbow plot: visualizes SD of each PC, search for where SD begins to plateau --------
elbow_plot <- Seurat::ElbowPlot(object = seurat_proteome,
                                ndims = 40) # Plateau only after about 30 PCs

# ggsave(elbow_plot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/proteomics_all_proteins/elbowplot_proteome.png",  width = 128, height = 90, units="mm")


# Visual inspection: 40 PCs

# Quantitative determination elbow plateau (two metrics, choose lowest value)--------------------------------

# Metric 1: New PC only contributes 5%, and all other cumulatively contribute 90%
pct <- seurat_proteome[["pca"]]@stdev / sum(seurat_proteome[["pca"]]@stdev) * 100 # Calculate percent of variation for each PC
cumu <- cumsum(pct) # Calculate cumulative percents with each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 # PC 41

# Metric 2: PC where percent change to next PC is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
co2 # PC 23

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Conclusion: use first 9 PCs to cluster (number less important with new versions of Seurat, biggest difference is time to compute with more PCs)
# Other option: just use first 40 clusters with sctransform normalization, since the more clusters, the more variation is accounted for (only downside: computing time)



################################################################################################################################################
####################################################      CLUSTERING      ######################################################################
################################################################################################################################################


# Graph-based clustering using K-nearest neighbor graph  ---------------------------------

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))

################################################################################################################################################
################################################      DIMENSIONALITY REDUCTION   ##############################################################
################################################################################################################################################

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6)

# Run T-SNE ----------------------------------------------------------------
seurat_proteome <- Seurat::RunTSNE(seurat_proteome, dims = 1:6)



################################################################################################################################################
#################################################     Choosing Resolution 0.4  ##############################################################
################################################################################################################################################

# Resolution 0.4 -------------------------------

# UMAP_res0.4 <-
Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 5,
    pt.size = 0.5,
    group.by = "RNA_snn_res.0.4"
) +
    ggplot2::ggtitle("Resolution 0.4") +
    # scale_color_viridis_d(option = "magma")
    ggplot2::scale_color_manual(values = c(
        # "#e5c494", "#60CEACFF", "#657060FF", "#A11A5BFF", "#F05B12FF", "#ffff33"
        "#657060FF",
        "#60CEACFF",
        "#e5c494",
        "#A11A5BFF",
        "#F05B12FF",
        "#ffff33"
    )) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle("UMAP proteomics") +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            colour = "black",
            size = 8
        ),
        axis.text = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           vjust = 3)
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_4/UMAP_proteomics_res_04.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )

# "#30123BFF" "#3E9BFEFF" "#46F884FF" "#E1DD37FF" "#F05B12FF" "#7A0403FF"
#
# "#0B0405FF" "#382A54FF" "#395D9CFF" "#3497A9FF" "#60CEACFF" "#DEF5E5FF"
#
# "#03051AFF" "#4C1D4BFF" "#A11A5BFF" "#E83F3FFF" "#F69C73FF" "#FAEBDDFF"

my_cols <- c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325")

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 5,
    pt.size = 0.3,
    cols = ggplot2::alpha(my_cols, 1),
    group.by = "fiber_type",
    order = c("Hybrid 2A/2X",
              "Type 2A",
              "Hybrid 1/2A",
              "Type 1")) +
    ggplot2::ggtitle("Resolution 0.4") +
    # scale_color_viridis_d(option = "magma")
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle("UMAP Proteomics - by fiber type") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        axis.text = ggplot2::element_text(size=6),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.text= ggplot2::element_text(size=4),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position = "none"
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_1/UMAP_proteomics_fiber_type.png"),
#     height = 55,
#     width = 65,
#     units = "mm"
# )

# UMAP subject ------------------------------------------------------------

seurat_proteome@meta.data <- seurat_proteome@meta.data |>
    dplyr::mutate(subject = dplyr::case_when(
        subject == "FOR2" ~ "P1",
        subject == "FOR4" ~ "P2",
        subject == "FOR9" ~ "P3",
        subject == "FOR10" ~ "P4",
        subject == "FOR11" ~ "P5",
        TRUE ~ "NAN"
    ))

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
#     here::here("doc/figures/figure_2/UMAP_proteomics_subject.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )

# Monocle Pseudotime ------------------------------------------------------

# Converting a Seurat object into monocle:

data_monocle <-
    SeuratWrappers::as.cell_data_set(seurat_proteome)

# Visualization of the object metadata:

head(SummarizedExperiment::colData(data_monocle))

# Double checking that the gene names are there:

rownames(monocle3::fData(data_monocle))[1:10]

# Recreating object clusters:

recreate.partitions <-
    SummarizedExperiment::colData(data_monocle) |>
    as.data.frame() |>
    dplyr::pull(seurat_clusters)

names(recreate.partitions) <- data_monocle@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

data_monocle@clusters@listData[["UMAP"]][["partitions"]] <-
    recreate.partitions

list.cluster <- seurat_proteome@active.ident

list.cluster <- seurat_proteome@active.ident
data_monocle@clusters@listData[["UMAP"]][["clusters"]] <-
    list.cluster
data_monocle@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-
    seurat_proteome@reductions$umap@cell.embeddings

# Showing cluster of cells before trajectory:

monocle3::plot_cells(
    data_monocle,
    color_cells_by = "cluster",
    label_groups_by_cluster = F,
    group_label_size = 10,
    cell_size = 0.75
) +
    ggplot2::scale_color_viridis_d(option = "plasma") +
    ggplot2::theme(legend.position = "right") +
    ggplot2::theme_minimal()

# Calculate pseudotime:

data_monocle <-
    monocle3::learn_graph(data_monocle, use_partition = F)

monocle3::plot_cells(
    data_monocle,
    color_cells_by = "cluster",
    label_groups_by_cluster = F,
    group_label_size = 7,
    cell_size = 1,
    label_cell_groups = FALSE,
    trajectory_graph_color = "black",
    alpha = 0.5,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE
) +
    ggplot2::scale_color_manual(values = c(
        "#657060FF",
        "#60CEACFF",
        "#e5c494",
        "#A11A5BFF",
        "#F05B12FF",
        "#ffff33"
    )) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
                   text = ggplot2::element_text(size = 8))

# ggplot2::ggsave(
#     here::here("doc/figures/figure_4/UMAP_pseudotime_proteomics.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )

seurat_clusters <- seurat_proteome@meta.data |>
    as.data.frame() |>
    dplyr::select(fiberID, RNA_snn_res.0.4) |>
    dplyr::rename("seurat_clusters" = RNA_snn_res.0.4) |>
    tibble::remove_rownames()

# readr::write_csv(seurat_clusters,
#                  here::here("data/proteomics clustering/seurat_clusters_6PC_res04.csv"))


# Feature plots -----------------------------------------------------------

feature_plot_MYH2 <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("MYH2"),
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

feature_plot_MYH7 <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("MYH7"),
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

feature_plot_UGDH <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("UGDH"),
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

feature_plot_RPL35 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("RPL35"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")


ggpubr::ggarrange(feature_plot_MYH2,
                  feature_plot_MYH7,
                  feature_plot_UGDH,
                  feature_plot_RPL35,
                  ncol = 2,
                  nrow = 2)

# ggplot2::ggsave(here::here("doc/figures/figure_4/feature_plots_proteomics.png"),
#                 width = 90,
#                 height = 60,
#                 units = "mm")

feature_plot_MYH1 <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("MYH1"),
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
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent') #transparent legend panel
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot <- ggpubr::ggarrange(feature_plot_MYH7,
                                  feature_plot_MYH2,
                                  feature_plot_MYH1,
                                  ncol = 3,
                                  nrow = 1,
                                  legend = "none") +
    ggplot2::ggtitle("Proteomics") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 8,
                                           face = "bold"),
    )

ggplot2::ggsave(here::here("doc/figures/umaps_now_in_fig1/feature_plot_MYHs.png"),
                width = 100,
                height = 33,
                units = "mm")

legend_feature <- ggpubr::get_legend(feature_plot_MYH7)
legend_feature <- ggpubr::as_ggplot(legend_feature)

# ggplot2::ggsave(here::here("doc/figures/figure_2/legend_MYH7_proteomics.png"),
#                 units = "mm",
#                 height = 5,
#                 width = 10)

legend_feature_MYH2 <- ggpubr::get_legend(feature_plot_MYH2)
legend_feature_MYH2 <- ggpubr::as_ggplot(legend_feature_MYH2)
legend_feature_MYH2

# ggplot2::ggsave(here::here("doc/figures/figure_2/legend_MYH2_proteomics.png"),
#                 units = "mm",
#                 height = 5,
#                 width = 10)

legend_feature_MYH1 <- ggpubr::get_legend(feature_plot_MYH1)
legend_feature_MYH1 <- ggpubr::as_ggplot(legend_feature_MYH1)
legend_feature_MYH1

# ggplot2::ggsave(here::here("doc/figures/figure_2/legend_MYH1_proteomics.png"),
#                 units = "mm",
#                 height = 5,
#                 width = 10)

# Feature plots RPL38 and RPS13 -------------------------------------------

feature_plot_RPL38 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("RPL38"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

# ggplot2::ggsave(here::here("doc/figures/figure_5/feature_plot_RPL38.png"),
#                 width = 40,
#                 height = 40,
#                 units = "mm")

feature_plot_RPS13 <- Seurat::FeaturePlot(seurat_proteome,
                                          features = c("RPS13"),
                                          pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

# ggplot2::ggsave(here::here("doc/figures/figure_5/feature_plot_RPS13.png"),
#                 width = 40,
#                 height = 40,
#                 units = "mm")



# feature plot UGDH -------------------------------------------------------

feature_plot_UGDH <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("UGDH"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin = ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_PHIP <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("PHIP"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin= ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_hist <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("HIST1H2AB"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5, vjust = 0.1),
        legend.position = "right",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.key.width = ggplot2::unit(1, "mm"),
        legend.spacing.x = ggplot2::unit(0.5, "mm"),
        legend.margin= ggplot2::margin(0,0,0,0),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

final_plot <- ggpubr::ggarrange(feature_plot_UGDH,
                                feature_plot_PHIP,
                                feature_plot_hist,
                                ncol = 3,
                                nrow = 1) +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggpubr::annotate_figure(final_plot, top = ggpubr::text_grob("Proteomics",
                                                            color = "black", face = "bold", size = 8))

# ggplot2::ggsave(here::here("doc/figures/figure_2/proteomics_feature_plots/feature_plots_proteomics.png"),
#                 units = "mm",
#                 width = 120,
#                 height = 35)

# Feature plot ANXA 11 ----------------------------------------------------

feature_plot_ANXA11 <- Seurat::FeaturePlot(seurat_proteome,
                                           features = c("ANXA11"),
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
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent', colour = "transparent") #transparent legend panel
        # legend.box.margin = ggplot2::element_blank()
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_PHIP <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("PHIP"),
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
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent', colour = "transparent") #transparent legend panel
        # legend.box.margin = ggplot2::element_blank()
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_TPPP <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("TPPP"),
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
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent', colour = "transparent") #transparent legend panel
        # legend.box.margin = ggplot2::element_blank()
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_hist <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("HIST1H2AB"),
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
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent', colour = "transparent") #transparent legend panel
        # legend.box.margin = ggplot2::element_blank()
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")















