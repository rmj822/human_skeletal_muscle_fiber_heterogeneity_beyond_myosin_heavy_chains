
# Loading data, ranking and extracting PCA drivers ------------------------

proteomics_data <- vroom::vroom(
    here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_proteomics_filtered.csv"),
    col_select = !c(1)
) |>
    as.data.frame() |>
    tibble::column_to_rownames("Gene.name")

pca_object <- vroom::vroom(here::here("data/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "...1") |>
    tibble::column_to_rownames("Gene.name") |>
    t() |>
    prcomp(scale = TRUE)

data_pca <- pca_object$rotation |>
    as.data.frame() |>
    dplyr::select(PC1, PC2)

selecting_hits <- round(nrow(data_pca) * 5 /100, 0)

PC1_enrichment_positive <- data_pca |>
    dplyr::select(PC1) |>
    dplyr::arrange(desc(PC1)) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::slice_head(n = selecting_hits) |>
    dplyr::pull(PC1, name = Gene.name)

PC1_enrichment_negative <- data_pca |>
    dplyr::select(PC1) |>
    dplyr::arrange(PC1) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::slice_head(n = selecting_hits) |>
    dplyr::pull(PC1, name = Gene.name)

PC2_enrichment_positive <- data_pca |>
    dplyr::select(PC2) |>
    dplyr::arrange(desc(PC2)) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::slice_head(n = selecting_hits) |>
    dplyr::pull(PC2, name = Gene.name)

PC2_enrichment_negative <- data_pca |>
    dplyr::select(PC2) |>
    dplyr::arrange(PC2) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::slice_head(n = selecting_hits) |>
    dplyr::pull(PC2, name = Gene.name)


# Overrepresentation analysis PC2 ---------------------------------------------------------

# Positive hits:

universe <- rownames(proteomics_data)
library(org.Hs.eg.db)

GSEA_PC2_positive <- clusterProfiler::enrichGO(
names(PC2_enrichment_positive),
ont = "ALL",
OrgDb = org.Hs.eg.db,
universe,
keyType = "SYMBOL",
pvalueCutoff = 0.05,
pAdjustMethod = "BH",
readable = TRUE
)

GSEA_PC2_positive <- clusterProfiler::simplify(
    GSEA_PC2_positive
)

clusterProfiler::dotplot(GSEA_PC2_positive,
                         font.size = 8,
                         showCategory = 5) +
    ggplot2::theme_bw()

# negative hits:

GSEA_PC2_negative <- clusterProfiler::enrichGO(
    names(PC2_enrichment_negative),
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    universe,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

GSEA_PC2_negative <- clusterProfiler::simplify(
    GSEA_PC2_negative
)

clusterProfiler::dotplot(GSEA_PC2_negative,
                         font.size = 8,
                         showCategory = 5) +
    ggplot2::theme_bw()


# Overrepresentation analysis PC1 ---------------------------------------------------------

# Positive hits:

GSEA_PC1_positive <- clusterProfiler::enrichGO(
    names(PC1_enrichment_positive),
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    universe,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

GSEA_PC1_positive <- clusterProfiler::simplify(
    GSEA_PC1_positive
)

clusterProfiler::dotplot(GSEA_PC1_positive,
                         font.size = 8,
                         showCategory = 5) +
    ggplot2::theme_bw()

# negative hits:

GSEA_PC1_negative <- clusterProfiler::enrichGO(
    names(PC1_enrichment_negative),
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    universe,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

GSEA_PC1_negative <- clusterProfiler::simplify(
    GSEA_PC1_negative
)

clusterProfiler::dotplot(GSEA_PC1_negative,
                         font.size = 8,
                         showCategory = 5) +
    ggplot2::theme_bw()

results_PC2_positive <- GSEA_PC2_positive@result
results_PC2_negative <- GSEA_PC2_negative@result
results_PC1_positive <- GSEA_PC1_positive@result
results_PC1_negative <- GSEA_PC1_negative@result

# Using over representation analysis independent of direction -------------

contrib__PC1 <- factoextra::get_pca_var(pca_object)$contrib |>
    as.data.frame() |>
    dplyr::select(Dim.1) |>
    dplyr::arrange(desc(Dim.1)) |>
    dplyr::slice_head(n = selecting_hits) |>
    tibble::rownames_to_column("genes") |>
    dplyr::pull(genes)

GSEA_PC1_contrib <- clusterProfiler::enrichGO(
    contrib__PC1,
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    universe,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

GSEA_PC1_contrib <- clusterProfiler::simplify(
    GSEA_PC1_contrib
)

clusterProfiler::dotplot(GSEA_PC1_contrib,
                         font.size = 8,
                         showCategory = 5) +
    ggplot2::theme_bw()

# PC2

contrib__PC2 <- factoextra::get_pca_var(pca_object)$contrib |>
    as.data.frame() |>
    dplyr::select(Dim.2) |>
    dplyr::arrange(desc(Dim.2)) |>
    dplyr::slice_head(n = selecting_hits) |>
    tibble::rownames_to_column("genes") |>
    dplyr::pull(genes)

GSEA_PC2_contrib <- clusterProfiler::enrichGO(
    contrib__PC2,
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    universe,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

GSEA_PC2_contrib <- clusterProfiler::simplify(
    GSEA_PC2_contrib
)

clusterProfiler::dotplot(GSEA_PC2_contrib,
                         font.size = 8,
                         showCategory = 5) +
    ggplot2::theme_bw()

#STILL NOT SURE IF WE WILL USE THE FOLLOWING:

# Using gprofiler2 for PC2 enrichment -------------------------------------

PC2_enrichment_positive <- as.data.frame(PC2_enrichment) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(PC2_enrichment > 0) |>
    dplyr::arrange(desc(PC2_enrichment)) |>
    dplyr::pull(Gene.name)


PC2_enrichment_positive <- gprofiler2::gost(
    PC2_enrichment_positive,
    organism = "hsapiens",
    exclude_iea = TRUE,
    custom_bg = rownames(proteomics_data),
    user_threshold = 0.05,
    correction_method = "fdr",
    ordered_query = TRUE
)

gprofiler2::gostplot(PC2_enrichment_positive)

# PC2 negative:

PC2_enrichment_negative <- as.data.frame(PC2_enrichment) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(PC2_enrichment < 0) |>
    dplyr::arrange(PC2_enrichment) |>
    dplyr::pull(Gene.name)


PC2_enrichment_negative <- gprofiler2::gost(
    PC2_enrichment_negative,
    organism = "hsapiens",
    exclude_iea = TRUE,
    custom_bg = rownames(proteomics_data),
    user_threshold = 0.05,
    correction_method = "fdr",
    ordered_query = TRUE
)

gprofiler2::gostplot(PC2_enrichment_negative)

# Using gprofiler2 for PC1 enrichment -------------------------------------

PC1_enrichment_positive <- as.data.frame(PC1_enrichment) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(PC1_enrichment > 0) |>
    dplyr::arrange(desc(PC1_enrichment)) |>
    dplyr::pull(Gene.name)


PC1_enrichment_positive <- gprofiler2::gost(
    PC1_enrichment_positive,
    organism = "hsapiens",
    exclude_iea = TRUE,
    custom_bg = rownames(proteomics_data),
    user_threshold = 0.05,
    correction_method = "fdr",
    ordered_query = TRUE
)

gprofiler2::gostplot(PC1_enrichment_positive)

# PC2 negative:

PC1_enrichment_negative <- as.data.frame(PC1_enrichment) |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::filter(PC1_enrichment < 0) |>
    dplyr::arrange(PC1_enrichment) |>
    dplyr::pull(Gene.name)


PC1_enrichment_negative <- gprofiler2::gost(
    PC1_enrichment_negative,
    organism = "hsapiens",
    exclude_iea = TRUE,
    custom_bg = rownames(proteomics_data),
    user_threshold = 0.05,
    correction_method = "fdr",
    ordered_query = TRUE
)

gprofiler2::gostplot(PC1_enrichment_negative)

# GSEA ranked by contribution ---------------------------------------------

data_pca <- vroom::vroom(here::here("data/data_pca_proteomics.csv")) |>
    dplyr::rename("Gene.name" = "...1") |>
    tibble::column_to_rownames("Gene.name") |>
    t() |>
    prcomp(scale = TRUE)

contrib_PC1 <- factoextra::get_pca_var(data_pca)$contrib |>
    as.data.frame() |>
    dplyr::select(Dim.1) |>
    tibble::rownames_to_column("Gene.names") |>
    dplyr::arrange(desc(Dim.1)) |>
    dplyr::pull(Dim.1, name = Gene.names)

GSEA_PC1 <- clusterProfiler::gseGO(
    contrib_PC1,
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    scoreType = "pos",
    pAdjustMethod = "fdr"
)

GSEA_PC1 <- clusterProfiler::simplify(
    GSEA_PC1
)

clusterProfiler::dotplot(GSEA_PC1,
                         split = ".sign",
                         font.size = 8,
                         showCategory = 10) +
    ggplot2::theme_bw()

# For PC2:

contrib_PC2 <- factoextra::get_pca_var(data_pca)$contrib |>
    as.data.frame() |>
    dplyr::select(Dim.2) |>
    tibble::rownames_to_column("Gene.names") |>
    dplyr::arrange(desc(Dim.2))

translating <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = contrib_PC2$Gene.names,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
)

contrib_PC2 <- translating |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::inner_join(contrib_PC2 |>
                          dplyr::rename("SYMBOL" = "Gene.names")) |>
    dplyr::pull(Dim.2, name = ENTREZID)

GSEA_PC2 <- clusterProfiler::gseGO(
    contrib_PC2,
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    scoreType = "pos",
    pAdjustMethod = "BH"
)

GSEA_PC2 <- clusterProfiler::simplify(
    GSEA_PC2
)

clusterProfiler::dotplot(GSEA_PC2,
                         split = ".sign",
                         font.size = 8,
                         showCategory = 10) +
    ggplot2::theme_bw()
