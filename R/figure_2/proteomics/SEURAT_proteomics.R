
# Loading data and creating Seurat object: --------------------------------

source(here::here("R/figure_1/MYH_curves.r"))

# Unfiltered data
data_proteomics <- vroom::vroom(here::here("data/data_proteomics_filtered.csv"))

# Filtered and processed data
data_pca <- vroom::vroom(here::here("data/data_pca_proteomics.csv")) |>
    dplyr::rename("Genes" = 1) |>
    tibble::column_to_rownames("Genes")

metadata <- vroom::vroom(
    "C:/Users/jns822/Desktop/Scripts/Heterofiber/data-raw/metadata.txt"
) |>
    dplyr::filter(fiberID %in% colnames(data_proteomics)) |>
    dplyr::filter(!duplicated(fiberID)) |>
    dplyr::inner_join(
        data_fiber_type |>
            dplyr::rename("fiberID" = "fiber_ID") |>
            dplyr::select(fiberID, fiber_type)
    ) |>
    dplyr::arrange(desc(fiberID))

seurat_object <- Seurat::CreateSeuratObject(counts = data_pca,
                                            project = "proteomics_seurat",
                                            assay = "LFQ",
                                            meta.data = metadata |>
                                                tibble::column_to_rownames("fiberID"))

# Seurat processing of filtered data --------------------------------------

seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst")

seurat_object <- Seurat::ScaleData(seurat_object)

seurat_object <- Seurat::RunPCA(seurat_object,
                                features = Seurat::VariableFeatures(object = seurat_object))

seurat_object <-
    Seurat::FindNeighbors(seurat_object,
                          dims = 1:6
    )

seurat_object <-
    Seurat::FindClusters(seurat_object,
                         resolution = 0.5
    )

seurat_object <- Seurat::RunUMAP(seurat_object,
                                 dims = 1:5
)

seurat_object <- Seurat::RunTSNE(seurat_object,
                                 dims = 1:6)

Seurat::DimPlot(seurat_object, reduction = "umap")

# Unfiltered data ---------------------------------------------------------
load("C:/Users/jns822/Desktop/Scripts/single_fiber_heterogeneity/data/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata")

#Extracting genes identified in transcriptomics data:
data_transcriptome <- filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@assays$SCT@data@Dimnames

transcriptomics_genes <- unlist(data_transcriptome[1])

#Filtering proteomics data with the genes identified in transcriptomics:
data_integration <- data_proteomics |>
    dplyr::select(!1) |>
    dplyr::filter(Gene.name %in% transcriptomics_genes) |>
    tibble::column_to_rownames("Gene.name") |>
    log2() |>
    limma::normalizeQuantiles() |>
    as.matrix() |>
    limma::removeBatchEffect( #Batch correction with limma
        batch = metadata$MS_batch,
        batch2 = metadata$digestion_batch
    ) |>
    PhosR::tImpute()


# Seurat workflow with integration data -----------------------------------

seurat_object <- Seurat::CreateSeuratObject(data_integration,
                                            project = "Proteomics integration",
                                            assay = "LFQ",
                                            meta.data = metadata |>
                                                tibble::column_to_rownames("fiberID"))

seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst")

seurat_object <- Seurat::ScaleData(seurat_object)

seurat_object <- Seurat::RunPCA(seurat_object,
                                features = Seurat::VariableFeatures(object = seurat_object))

Seurat::DimPlot(seurat_object, reduction = "pca", group.by = "fiber_type")

seurat_object <-
    Seurat::FindNeighbors(seurat_object,
                          dims = 1:10
    )

seurat_object <-
    Seurat::FindClusters(seurat_object,
                         resolution = 0.5
    )

seurat_object <- Seurat::RunUMAP(seurat_object,
                                 dims = 1:50)

seurat_object <- Seurat::RunTSNE(seurat_object,
                                 dims = 1:50)

Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "fiber_type")


# Exporting dataframes for Seurat integration -----------------------------

proteomics_stringent <- data_pca |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes %in% transcriptomics_genes)

readr::write_csv(proteomics_stringent,
                 here::here("data/seurat_data_integration/proteomics_stringent.csv"))

proteomics_non_filtered <- data_proteomics |>
    dplyr::select(!1) |>
    dplyr::filter(Gene.name %in% transcriptomics_genes) |>
    tibble::column_to_rownames("Gene.name") |>
    log2() |>
    limma::normalizeQuantiles() |>
    as.matrix() |>
    limma::removeBatchEffect( #Batch correction with limma
        batch = metadata$MS_batch,
        batch2 = metadata$digestion_batch
    ) |>
    PhosR::tImpute() |>
    as.data.frame() |>
    tibble::rownames_to_column("Genes")

readr::write_csv(proteomics_non_filtered,
                 here::here("data/seurat_data_integration/proteomics_non_stringent.csv"))
