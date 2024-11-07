
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
    dplyr::mutate(
        PC2 = PC2 * -1,
        PC1 = PC1 * -1
    ) |>
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
                height = 60,
                width = 195)

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

# ggplot2::ggsave(here::here("doc/figures/figure_2/figure_2F.png"),
#                 units = "mm",
#                 width = 120,
#                 height = 35)


################################################################################################################################################
###############################################      FIGURE 2G-H   ##########################################################
################################################################################################################################################

data_proteomics <-read.csv(here::here("data/figure_2/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "X") |>
    tibble::column_to_rownames("Gene.name")

pca_object <- prcomp(t(data_proteomics), scale = T)

PC_drivers_proteomics <- pca_object$rotation |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    tibble::rownames_to_column("Gene.name")


PC_drivers_transcriptomics <-
    filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@reductions$pca@feature.loadings |>
    as.data.frame() |>
    dplyr::select(PC_1, PC_2) |>
    tibble::rownames_to_column("Gene.name")

selection_vector <- PC_drivers_proteomics |>
    dplyr::pull(Gene.name)

PC_drivers <- PC_drivers_transcriptomics |>
    dplyr::filter(Gene.name %in% selection_vector) |>
    dplyr::inner_join(PC_drivers_proteomics) |>
    tibble::column_to_rownames("Gene.name") |>
    limma::normalizeBetweenArrays() |>
    as.data.frame() |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::rename(PC1_transcriptomics = PC_1,
                  PC2_transcriptomics = PC_2,
                  PC1_proteomics = PC1,
                  PC2_proteomics = PC2) |>
    dplyr::mutate(PC2_transcriptomics = PC2_transcriptomics * -1) %>%
    dplyr::mutate(PC1_proteomics = -PC1_proteomics)

Genes <- PC_drivers$Gene.name

PC_drivers |>
    dplyr::mutate(label = dplyr::case_when(
        Gene.name %in% c(
            "RPL38",
            "RPL31",
            "ANXA11",
            "RPL27",
            "RPS9",
            "PSMA1",
            "ME1",
            "TNNC2",
            "SAR1B",
            "RPL18",
            "TRIM72"
        ) ~ Gene.name,
        TRUE ~ ""
    )) |>
    ggplot2::ggplot(ggplot2::aes(x = PC1_transcriptomics,
                                 y = PC1_proteomics,
                                 label = label)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.65) +
    ggrepel::geom_label_repel(size = 1.8,
                              label.padding=0.1,
                              max.overlaps = Inf,
                              min.segment.length=0.1,
                              segment.size=0.2,
                              force = 46) +
    ggplot2::theme_minimal() +
    ggplot2::annotate(geom = "text",
                      x = - 0.07,
                      y = 0.06,
                      label = paste("r = -0.027"),
                      size = 2.5) +
    ggplot2::annotate(geom = "text",
                      x = - 0.07,
                      y = 0.05,
                      label = paste("p = 0.315"),
                      size = 2.5) +
    ggplot2::xlab("PC1 Transcriptomics") +
    ggplot2::ylab("PC1 Proteomics") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7, face = "bold")
    )

ggsave(here::here("doc/figures/figure_2/figure_2H.png"),
       width = 60,
       height = 60,
       units="mm")



# Correlation PC 2 --------------------------------------------------------

PC_drivers |>
    dplyr::mutate(PC2_proteomics = -PC2_proteomics) %>%
    dplyr::mutate(label = dplyr::case_when(
        Gene.name %in% c(
            "TNNT1",
            "MYH7",
            "ATP2A2",
            "MYH7B",
            "DHRS7C",
            "MYLPF",
            "MYH2",
            "MYH1",
            "ACTA1"
        ) ~ Gene.name,
        TRUE ~ ""
    )) |>
    ggplot2::ggplot(ggplot2::aes(x = PC2_transcriptomics,
                                 y = PC2_proteomics,
                                 label = label)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.65) +
    ggrepel::geom_label_repel(size = 1.8,
                              label.padding=0.1,
                              max.overlaps = Inf,
                              min.segment.length=0.1,
                              segment.size=0.2,
                              force = 45) +
    ggplot2::theme_minimal() +
    ggplot2::annotate(geom = "text",
                      x = - 0.05,
                      y = 0.06,
                      label = paste("r = 0.663"),
                      size = 2.5) +
    ggplot2::annotate(geom = "text",
                      x = - 0.05,
                      y = 0.05,
                      label = paste("p = 5.81e-173"),
                      size = 2.5) +
    ggplot2::xlab("PC2 Transcriptomics") +
    ggplot2::ylab("PC2 Proteomics") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7, face = "bold")
    )

ggsave(here::here("doc/figures/figure_2/figure_2G.png"),
       width = 60,
       height = 60,
       units="mm")









