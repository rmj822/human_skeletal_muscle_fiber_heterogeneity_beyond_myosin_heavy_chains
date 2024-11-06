
# Load required packages --------------------------------------------------
library(tidyverse)
library(ggrastr)
library(Seurat)
library(SeuratData)
library(Matrix)

# Load in gene annotation file --------------------------------------------
annotations <- read.csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

#################################################################################################################################################
#########################################      LOAD TRANSCRIPTOME AND PROTEOME DATA  ############################################################
#################################################################################################################################################

# Transcriptome -----------------------------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata") # 925 fibers for 7418 genes

# Transcriptome metadata
transcriptome_metadata <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data

# Proteome -----------------------------------------------------------
data_proteomics <- read.csv("~/single_fiber_heterogeneity/data/data_proteomics_filtered.csv") # 974 fibers for 2983 proteins
data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = "Gene.name") |>
    tibble::column_to_rownames("Protein")

# Proteome metadata
proteomics_metadata <- read.csv("~/single_fiber_heterogeneity/data/metadata_proteomics.csv")


proteomics_metadata <- proteomics_metadata |>  dplyr::rename("fiberID" = X)
rownames(proteomics_metadata) <- proteomics_metadata$fiberID


#################################################################################################################################################
#########################################      FILTER TRANSCRIPTOME FOR OVERLAPPING GENES/PROTEINS  ###########################
#################################################################################################################################################

# Filter transcriptome data
counts_transcriptome <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")
counts_transcriptome <- as.data.frame(counts_transcriptome)

counts_transcriptome$Gene <- rownames(counts_transcriptome)
counts_transcriptome <- counts_transcriptome %>%
    dplyr::filter(Gene %in% rownames(data_proteomics)) %>%
    dplyr::select(-Gene) # 925 for 2243 genes

#################################################################################################################################################
###############################################      TRANSCRIPTOMICS SEURAT DATA PROCESSING  ####################################################
#################################################################################################################################################

# Create seurat object
seurat_transcriptome <- Seurat::CreateSeuratObject(counts_transcriptome)

# Add metadata
all(colnames(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) %in% colnames(seurat_transcriptome))
all(colnames(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) == colnames(seurat_transcriptome))

seurat_transcriptome$subject <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest$subject
seurat_transcriptome$fiber <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest$fiber
seurat_transcriptome$condition <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest$condition
seurat_transcriptome$fiber_type_MYH <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest$fiber_type_MYH
seurat_transcriptome$fiber_type_MYH_hybrids <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest$fiber_type_MYH_hybrids
seurat_transcriptome$final_cluster <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest$final_cluster

# Normalize data with NormalizeData
seurat_transcriptome <- NormalizeData(seurat_transcriptome, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
seurat_transcriptome <- FindVariableFeatures(seurat_transcriptome, selection.method = "vst", nfeatures = 2000)

# Data scaling (mean expression = 0, variance across cells = 1)
genes_scale_transcriptome <- rownames(seurat_transcriptome)
seurat_transcriptome <- ScaleData(seurat_transcriptome, features = genes_scale_transcriptome, verbose = TRUE)

# Perform PCA reduction
seurat_transcriptome <- RunPCA(seurat_transcriptome, features = VariableFeatures(object = seurat_transcriptome), verbose = TRUE)

PCA_transcriptome_PC1_PC2 <- DimPlot(seurat_transcriptome, reduction = "pca", group.by = "fiber_type_MYH_hybrids") # No separation of fiber type by PC1, limited separation by PC2
PCA_transcriptome_PC1_PC3 <- DimPlot(seurat_transcriptome, reduction = "pca", dims = c(1,3), group.by = "fiber_type_MYH_hybrids") # Good separation by PC3

# Determine dimensionality of dataset
ElbowPlot(seurat_transcriptome) # 10 PCs?

# Clustering
seurat_transcriptome <- FindNeighbors(seurat_transcriptome, dims = 1:10, verbose = TRUE)
seurat_transcriptome <- FindClusters(seurat_transcriptome, resolution = c(0.2, 0.4, 0.6, 0.8), verbose = TRUE)

# Run UMAP
seurat_transcriptome <- Seurat::RunUMAP(seurat_transcriptome, dims = 1:10)

# Run T-SNE
seurat_transcriptome <- RunTSNE(seurat_transcriptome, dims = 1:10)

# Check UMAPs
UMAP_transcriptome_res0.2 <- DimPlot(seurat_transcriptome, reduction = "umap", group.by = "RNA_snn_res.0.2")
UMAP_transcriptome_res0.4 <- DimPlot(seurat_transcriptome, reduction = "umap", group.by = "RNA_snn_res.0.4")
UMAP_transcriptome_res0.6 <- DimPlot(seurat_transcriptome, reduction = "umap", group.by = "RNA_snn_res.0.6")
UMAP_transcriptome_res0.8 <- DimPlot(seurat_transcriptome, reduction = "umap", group.by = "RNA_snn_res.0.8")
UMAP_transcriptome_fibertype <- DimPlot(seurat_transcriptome, reduction = "umap", group.by = "fiber_type_MYH_hybrids")
UMAP_transcriptome_origcluster <- DimPlot(seurat_transcriptome, reduction = "umap", group.by = "final_cluster")


################################################################################################################################################
###########################################      PERFORM PCA WITH PRCOMP      ##################################################################
################################################################################################################################################

# Extract scaled data
transcriptomics_scaledata <- GetAssayData(object = seurat_transcriptome, assay = "RNA", slot = "scale.data")

# Extract non log-transformed counts
transcriptomics_counts <- GetAssayData(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, assay = "SCT", slot = "counts")

pca_object <- prcomp(t(transcriptomics_scaledata),  center = F, scale. = F)

data_pca <- pca_object$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2, PC3)

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
################################################      BY GO RIBOSOMAL TERM      ################################################################
################################################################################################################################################

vector_ribosome <- meta_go2gene[grep("cytosolic ribosome",
                                     meta_go2gene$GO),
                                "SYMBOL"] |>
    unique()

ribosome_filtered <- transcriptomics_counts |>
    as.data.frame() %>%
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(Gene.name %in% vector_ribosome) |>
    dplyr::select(!Gene.name)

total_intensity <- colSums(transcriptomics_counts, na.rm = T)

rel_abundance_ribosomes <- (colSums(ribosome_filtered,
                                    na.rm = T) / total_intensity) * 100

data_pca_ribosomes <- as.data.frame(rel_abundance_ribosomes) |>
    tibble::rownames_to_column("fiberID") |>
    dplyr::inner_join(data_pca)

data_ribosomal_proteins <- factoextra::get_pca_var(pca_object)$coor |>
    as.data.frame() |>
    dplyr::select(Dim.1, Dim.2, Dim.3) |>
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


custom_ribo_plot_with_legend <- data_pca_ribosomes |>
    ggplot2::ggplot(ggplot2::aes(
        x = PC1,
        y = PC3,
        color = rel_abundance_ribosomes
    )) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_viridis_c("Rel. \nAbundance", option = "plasma") +
    ggplot2::geom_segment(
        data = data_ribosomal_proteins,
        ggplot2::aes(x = x.centroid,
                     y = y.centroid,
                     xend = Dim.1,
                     yend = Dim.3),
        color = "black",
        arrow = grid::arrow(type = "closed",
                            length = unit(1, "mm"))
    ) +
    ggrepel::geom_label_repel(data = data_ribosomal_proteins,
                              ggplot2::aes(
                                  x = Dim.1,
                                  y = Dim.3,
                                  label = label),
                              size = 2,
                              color = "black",
                              label.padding = 0.1,
                              min.segment.length = 0.1,
                              segment.size = 0.2,
                              force = 20,
                              max.overlaps = Inf) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA colored by Cytosolic Ribosome GO term") +
    ggplot2::xlab("PC1 (XX %)") +
    ggplot2::ylab("PC3 (XX %)") +
    ggplot2::theme(
        text = element_text(face="bold", colour="black", size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 3),
        legend.position = "bottom",
        legend.key.height = ggplot2::unit(2, "mm")
    ) +
    scale_x_continuous(trans = "reverse",
                       breaks=c(-50, 0, 50, 100),
                       labels=c("50", "0", "-50", "-100")
    )

ggsave(custom_ribo_plot_with_legend,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_5/PCA_transcriptomics_ribosomes_overlapping_genes.png",
       width = 90,
       height = 60,
       units="mm")




