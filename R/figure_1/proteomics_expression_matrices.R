

# supplementary table S25. Metadata heterofiber -------------------------------------------------

sup_table_metadata <-
    vroom::vroom(here::here("data/metadata_proteomics.csv")) |>
    dplyr::select(subject,
                  fiber_number,
                  plate_number,
                  plate_position,
                  digestion_batch,
                  MS_batch,
                  date_isolation) |>
    dplyr::mutate(subject = dplyr::case_when(
        subject == "FOR2" ~ "P1",
        subject == "FOR4" ~ "P2",
        subject == "FOR9" ~ "P3",
        subject == "FOR10" ~ "P4",
        subject == "FOR11" ~ "P5",
        TRUE ~ "unknown"
    )) |>
    dplyr::mutate(participant = subject,
                  number = fiber_number) |>
    tidyr::unite(
        fiberID,
        participant,
        number,
        sep = "_fibNumber"
    ) |>
    dplyr::select(fiberID, everything())

write.csv(
    sup_table_metadata,
    here::here("doc/supplementary_tables/proteomics/metadata_heterofiber.csv")
)


# sup table 27. metadata MD -----------------------------------------------

metadata_MD <- vroom::vroom(here::here("data/metadata_MD_w_fiber_type_w_anonim.csv")) |>
    dplyr::mutate("collection_order" = as.character(collection_order)) |>
    tidyr::unite(
        col = "new_fiber_ID",
        anonimized_subject,
        collection_order,
        sep = "_"
    ) |>
    dplyr::select(
new_fiber_ID,
Plate_position,
Plate_number,
condition,
date_isolation
    ) |>
    dplyr::mutate(fiber_ID = new_fiber_ID,
                  condition = dplyr::case_when(
                      condition == "troponin" ~ "TNNT1_nemaline_myopaty",
                      condition == "actin" ~ "ACTA1_nemaline_myopaty",
                      TRUE ~ "control"
                  )) |>
    tidyr::separate(new_fiber_ID,
                    into = c(
                        "subject",
                        "fiber_number"
                    ),
                    sep = "_") |>
    dplyr::select(
        fiber_ID,
        subject,
        fiber_number,
        condition,
        Plate_number,
        Plate_position,
        date_isolation
    )

write.csv(
    metadata_MD,
    here::here("doc/supplementary_tables/proteomics/metadata_muscle_disease.csv")
)

# Supplementary table 2. Raw expression --------------------------

sup_table_proteomics_raw <- vroom::vroom(here::here(
  "data-raw/data_raw_supp_table.csv"
)) |>
  dplyr::select(!1) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("fiberID") |>
  tidyr::separate(
    fiberID,
    into = c(
      "subject",
      "sample"
    ),
    sep = "_"
  ) |>
  dplyr::mutate(subject = dplyr::case_when(
    subject == "FOR2" ~ "P1",
    subject == "FOR4" ~ "P2",
    subject == "FOR9" ~ "P3",
    subject == "FOR10" ~ "P4",
    subject == "FOR11" ~ "P5",
    TRUE ~ "unknown"
  )) |>
  tidyr::unite(
    fiberID,
    subject,
    sample,
    sep = "_"
  ) |>
  tibble::column_to_rownames("fiberID") |>
  t() |>
  as.data.frame() |>
  dplyr::rename("Gene_name" = "unknown_NA") |>
  dplyr::select(
    Gene_name,
    everything()
  ) |>
  dplyr::mutate(dplyr::across(
    .cols = 2:1039,
    as.numeric
  ))

write.csv(
  sup_table_proteomics_raw,
  here::here("doc/supplementary_tables/proteomics/proteomics_raw_expression.csv")
)


# Supplementary table 5. Filtered values ----------------------------------

sup_table_proteomics_filtered <- read.csv(here::here("data/data_pca_proteomics.csv")) |>
  dplyr::rename("Gene_name" = "X") |>
  tibble::column_to_rownames("Gene_name") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("fiberID") |>
  tidyr::separate(
    fiberID,
    into = c(
      "subject",
      "sample"
    ),
    sep = "_"
  ) |>
  dplyr::mutate(subject = dplyr::case_when(
    subject == "FOR2" ~ "P1",
    subject == "FOR4" ~ "P2",
    subject == "FOR9" ~ "P3",
    subject == "FOR10" ~ "P4",
    subject == "FOR11" ~ "P5",
    TRUE ~ "unknown"
  )) |>
  tidyr::unite(
    fiberID,
    subject,
    sample,
    sep = "_"
  ) |>
  tibble::column_to_rownames("fiberID") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Gene_name")

write.csv(
  sup_table_proteomics_filtered,
  here::here("doc/supplementary_tables/proteomics/proteomics_filtered_expression.csv")
)


# Supplementary table 8. PCA loadings -------------------------------------

load(here::here("R/figure_2/PCA_proteomics.r"))

pca_loadings <- pca_object$rotation |>
  as.data.frame() |>
  dplyr::select(PC1, PC2) |>
  tibble::rownames_to_column("Gene_name")

write.csv(
  pca_loadings,
  here::here("doc/supplementary_tables/proteomics/proteomics_pca_loadings.csv")
)


# GSEA_results_PCA --------------------------------------------------------

results_PC1_negative <- vroom::vroom(
  here::here("data/GSEA_PCA_proteomics/results_PC1_negative.csv")
) |>
  tibble::column_to_rownames("...1") |>
  dplyr::mutate("enrichment" = "PC1_negative") |>
  dplyr::select(enrichment, everything())

results_PC2_negative <- vroom::vroom(
  here::here("data/GSEA_PCA_proteomics/results_PC2_negative.csv")
) |>
  tibble::column_to_rownames("...1") |>
  dplyr::mutate("enrichment" = "PC2_negative") |>
  dplyr::select(enrichment, everything())

results_PC1_positive <- vroom::vroom(
  here::here("data/GSEA_PCA_proteomics/results_PC1_positive.csv")
) |>
  tibble::column_to_rownames("...1") |>
  dplyr::mutate("enrichment" = "PC1_positive") |>
  dplyr::select(enrichment, everything())

results_PC2_positive <- vroom::vroom(
  here::here("data/GSEA_PCA_proteomics/results_PC2_positive.csv")
) |>
  tibble::column_to_rownames("...1") |>
  dplyr::mutate("enrichment" = "PC2_positive") |>
  dplyr::select(enrichment, everything())

proteomics_GSEA_results <- dplyr::bind_rows(
  results_PC1_negative,
  results_PC1_positive,
  results_PC2_negative,
  results_PC2_positive
) |>
  tibble::remove_rownames()

write.csv(
  proteomics_GSEA_results,
  here::here("doc/supplementary_tables/proteomics/proteomics_GSEA_results.csv")
)


# Supplementary table 12. Pseudobulk data ----------------------------------

data_pseudobulk <- vroom::vroom(here::here("data/data_proteomics_pseudobulk_fastvslow.csv")) |>
  dplyr::rename("Gene_name" = "...1")

write.csv(
  data_pseudobulk,
  here::here("doc/supplementary_tables/proteomics/proteomics_data_pseudobulk_slowvfast.csv")
)


# Supplementary table 14. DE matrix ----------------------------------------

DE_results <- vroom::vroom(here::here("data/DE_analysis_slow_vs_fast_proteomics.csv")) |>
  dplyr::rename("Gene_name" = "...1") |>
    dplyr::mutate(signifficant = dplyr::case_when(
        adj.P.Val < 0.05 ~ "YES",
        TRUE ~ "NO"
    ))

write.csv(
  DE_results,
  here::here("doc/supplementary_tables/proteomics/proteomics_DE_data_slowvfast.csv")
)


# Supplementary table 6. GSEA slow versus fast ----------------------------


# Supplementary table 7. Ribosomal proteins -------------------------------

ribosomal_clusters <- vroom::vroom("data/ribosomal_clusters.csv") |>
    dplyr::mutate("Gene_name" = Gene.name)

ribosomal_proteins <- vroom::vroom(here::here("data/data_ribosomes_proteomics.csv")) |>
  dplyr::rename("Gene_name" = "...1") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("fiberID") |>
  tidyr::separate(
    fiberID,
    into = c(
      "subject",
      "sample"
    ),
    sep = "_"
  ) |>
  dplyr::mutate(subject = dplyr::case_when(
    subject == "FOR2" ~ "P1",
    subject == "FOR4" ~ "P2",
    subject == "FOR9" ~ "P3",
    subject == "FOR10" ~ "P4",
    subject == "FOR11" ~ "P5",
    TRUE ~ "Gene"
  )) |>
  tidyr::unite(
    fiberID,
    subject,
    sample,
    sep = "_"
  ) |>
  tibble::column_to_rownames("fiberID") |>
  t() |>
  as.data.frame() |>
    dplyr::inner_join(ribosomal_clusters |>
                          dplyr::select(Gene_name,
                                        ribosomal_clusters)) |>
    dplyr::select(ribosomal_clusters, everything()) |>
  tibble::remove_rownames()

write.csv(
  ribosomal_proteins,
  here::here("doc/supplementary_tables/proteomics/proteomics_ribosomal_data.csv")
)


# Supplementary table 8. Raw expression MD --------------------------------


MD_raw_expression <- vroom::vroom("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_muscle_disease.csv") |>
  tibble::column_to_rownames("...1")

metadata_MD <- vroom::vroom(here::here("data/metadata_MD_w_fiber_type_w_anonim.csv")) |>
  dplyr::mutate("collection_order" = as.character(collection_order)) |>
  tidyr::unite(
    col = "new_sample_ID",
    anonimized_subject,
    collection_order,
    sep = "_"
  )

MD_raw_expression <- MD_raw_expression |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("fiber_ID") |>
  dplyr::inner_join(metadata_MD |>
    dplyr::select(
      fiber_ID,
      new_sample_ID
    )) |>
  dplyr::select(!fiber_ID) |>
  tibble::column_to_rownames("new_sample_ID") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Gene_name")

write.csv(
  MD_raw_expression,
  here::here("doc/supplementary_tables/MD_raw_expression.csv")
)


# Supplementary table 9. Filtered expression MD ---------------------------

MD_filtered_expression <- vroom::vroom(here::here("data/MD_data_PCA.csv")) |>
  tibble::column_to_rownames("...1") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("fiber_ID") |>
  dplyr::inner_join(metadata_MD |>
    dplyr::select(
      fiber_ID,
      new_sample_ID
    )) |>
  dplyr::select(!fiber_ID) |>
  tibble::column_to_rownames("new_sample_ID") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Gene_name")

write.csv(
  MD_filtered_expression,
  here::here("doc/supplementary_tables/MD_filtered_expression.csv")
)


# Supplementary table 10. Pseudobulk MD -----------------------------------


MD_pseudobulk_expression <- vroom::vroom(here::here("data/data_MD_pseudobulk.csv")) |>
  dplyr::rename("Gene_name" = "Genes")

write.csv(
  MD_pseudobulk_expression,
  here::here("doc/supplementary_tables/MD_pseudobulk.csv")
)


# Supplementary table 21. MD differential expression analysis --------

MD_DE_results_ACTA1_vs_control <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls.csv")) |>
    dplyr::rename(
        LogFC_ACTA1_vs_control = logFC,
        P.Value_ACTA1_vs_control = P.Value,
        adj.P.Val_ACTA1_vs_control = adj.P.Val,
        xiao_ACTA1_vs_control = xiao,
        signifficant_ACTA1_vs_control = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_ACTA1_vs_control,
                  P.Value_ACTA1_vs_control,
                  adj.P.Val_ACTA1_vs_control,
                  xiao_ACTA1_vs_control,
                  signifficant_ACTA1_vs_control)

MD_DE_results_TNNT1_vs_control <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls.csv")) |>
    dplyr::rename(
        LogFC_TNNT1_vs_control = logFC,
        P.Value_TNNT1_vs_control = P.Value,
        adj.P.Val_TNNT1_vs_control = adj.P.Val,
        xiao_TNNT1_vs_control = xiao,
        signifficant_TNNT1_vs_control = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_TNNT1_vs_control,
                  P.Value_TNNT1_vs_control,
                  adj.P.Val_TNNT1_vs_control,
                  xiao_TNNT1_vs_control,
                  signifficant_TNNT1_vs_control)

MD_DE_results_ACTA1_vs_TNNT1 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin.csv")) |>
    dplyr::rename(
        LogFC_ACTA1_vs_TNNT1 = logFC,
        P.Value_ACTA1_vs_TNNT1 = P.Value,
        adj.P.Val_ACTA1_vs_TNNT1 = adj.P.Val,
        xiao_ACTA1_vs_TNNT1 = xiao,
        signifficant_ACTA1_vs_TNNT1 = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_ACTA1_vs_TNNT1,
                  P.Value_ACTA1_vs_TNNT1,
                  adj.P.Val_ACTA1_vs_TNNT1,
                  xiao_ACTA1_vs_TNNT1,
                  signifficant_ACTA1_vs_TNNT1)

MD_DE_results <- MD_DE_results_ACTA1_vs_control |>
    dplyr::inner_join(
        MD_DE_results_TNNT1_vs_control
    ) |>
    dplyr::inner_join(
        MD_DE_results_ACTA1_vs_TNNT1
    )

write.csv(
  MD_DE_results,
  here::here("doc/supplementary_tables/proteomics/MD_DE_results.csv")
)


# Supplementary table 23. MD DE analysis fast and slow -----------------------

# Type 1:

MD_DE_results_ACTA1_vs_control_type_1 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls_type_1.csv")) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::rename(
        LogFC_ACTA1_vs_control = logFC,
        P.Value_ACTA1_vs_control = P.Value,
        adj.P.Val_ACTA1_vs_control = adj.P.Val,
        xiao_ACTA1_vs_control = xiao,
        signifficant_ACTA1_vs_control = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_ACTA1_vs_control,
                  P.Value_ACTA1_vs_control,
                  adj.P.Val_ACTA1_vs_control,
                  xiao_ACTA1_vs_control,
                  signifficant_ACTA1_vs_control)

MD_DE_results_TNNT1_vs_control_type_1 <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls_type_1.csv")) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::rename(
        LogFC_TNNT1_vs_control = logFC,
        P.Value_TNNT1_vs_control = P.Value,
        adj.P.Val_TNNT1_vs_control = adj.P.Val,
        xiao_TNNT1_vs_control = xiao,
        signifficant_TNNT1_vs_control = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_TNNT1_vs_control,
                  P.Value_TNNT1_vs_control,
                  adj.P.Val_TNNT1_vs_control,
                  xiao_TNNT1_vs_control,
                  signifficant_TNNT1_vs_control)


MD_DE_results_ACTA1_vs_TNNT1_type_1 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin_type_1.csv")) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::rename(
        LogFC_ACTA1_vs_TNNT1 = logFC,
        P.Value_ACTA1_vs_TNNT1 = P.Value,
        adj.P.Val_ACTA1_vs_TNNT1 = adj.P.Val,
        xiao_ACTA1_vs_TNNT1 = xiao,
        signifficant_ACTA1_vs_TNNT1 = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_ACTA1_vs_TNNT1,
                  P.Value_ACTA1_vs_TNNT1,
                  adj.P.Val_ACTA1_vs_TNNT1,
                  xiao_ACTA1_vs_TNNT1,
                  signifficant_ACTA1_vs_TNNT1)


# Type 2A:

MD_DE_results_ACTA1_vs_control_type_2 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls_type_2.csv")) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::rename(
        LogFC_ACTA1_vs_control = logFC,
        P.Value_ACTA1_vs_control = P.Value,
        adj.P.Val_ACTA1_vs_control = adj.P.Val,
        xiao_ACTA1_vs_control = xiao,
        signifficant_ACTA1_vs_control = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_ACTA1_vs_control,
                  P.Value_ACTA1_vs_control,
                  adj.P.Val_ACTA1_vs_control,
                  xiao_ACTA1_vs_control,
                  signifficant_ACTA1_vs_control)

MD_DE_results_TNNT1_vs_control_type_2 <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls_type_2.csv")) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::rename(
        LogFC_TNNT1_vs_control = logFC,
        P.Value_TNNT1_vs_control = P.Value,
        adj.P.Val_TNNT1_vs_control = adj.P.Val,
        xiao_TNNT1_vs_control = xiao,
        signifficant_TNNT1_vs_control = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_TNNT1_vs_control,
                  P.Value_TNNT1_vs_control,
                  adj.P.Val_TNNT1_vs_control,
                  xiao_TNNT1_vs_control,
                  signifficant_TNNT1_vs_control)


MD_DE_results_ACTA1_vs_TNNT1_type_2 <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin_type_2.csv")) |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not signifficant"
    )) |>
    dplyr::rename(
        LogFC_ACTA1_vs_TNNT1 = logFC,
        P.Value_ACTA1_vs_TNNT1 = P.Value,
        adj.P.Val_ACTA1_vs_TNNT1 = adj.P.Val,
        xiao_ACTA1_vs_TNNT1 = xiao,
        signifficant_ACTA1_vs_TNNT1 = signifficant
    ) |>
    dplyr::select(Genes,
                  LogFC_ACTA1_vs_TNNT1,
                  P.Value_ACTA1_vs_TNNT1,
                  adj.P.Val_ACTA1_vs_TNNT1,
                  xiao_ACTA1_vs_TNNT1,
                  signifficant_ACTA1_vs_TNNT1)


MD_DE_results_fast <- MD_DE_results_ACTA1_vs_control_type_2 |>
    dplyr::inner_join(
        MD_DE_results_TNNT1_vs_control_type_2
    ) |>
    dplyr::inner_join(MD_DE_results_ACTA1_vs_TNNT1_type_2) |>
    dplyr::mutate("fiber_type" = "type_2A")


MD_DE_results_fiber_type <- MD_DE_results_slow |>
    dplyr::rows_append(
        MD_DE_results_fast
    )


write.csv(
  MD_DE_results_fiber_type,
  here::here("doc/supplementary_tables/proteomics/MD_DE_results_fiber_type.csv")
)

# Supplementary table 22. MD GSEA fiber type ------------------------------

# Enrichment actin vs controls:
result_GSEA_actin_control <- vroom::vroom(
  here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_control.csv")
) |>
  dplyr::select(!1) |>
  dplyr::mutate("enrichment" = dplyr::case_when(
    enrichmentScore < 0 ~ "controls",
    enrichmentScore > 0 ~ "actin",
    TRUE ~ "unknown"
  ))

# Enrichment actin vs controls in type 1:
result_GSEA_actin_control_type_1 <- vroom::vroom(
  here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_control_type_1.csv")
) |>
  dplyr::select(!1) |>
  dplyr::mutate("enrichment" = dplyr::case_when(
    enrichmentScore < 0 ~ "controls_type_1",
    enrichmentScore > 0 ~ "actin_type_1",
    TRUE ~ "unknown"
  ))

# Enrichment actin vs controls in type 2:
result_GSEA_actin_control_type_2 <- vroom::vroom(
  here::here("data/MD_GSEA_fiber_type/result_GSEA_actin_control_type_2.csv")
) |>
  dplyr::select(!1) |>
  dplyr::mutate("enrichment" = dplyr::case_when(
    enrichmentScore < 0 ~ "controls_type_2",
    enrichmentScore > 0 ~ "actin_type_2",
    TRUE ~ "unknown"
  ))

# Enrichment troponin vs controls:
result_GSEA_troponin_control <- vroom::vroom(
  here::here("data/MD_GSEA_fiber_type/result_GSEA_troponin_control.csv")
) |>
  dplyr::select(!1) |>
  dplyr::mutate("enrichment" = dplyr::case_when(
    enrichmentScore < 0 ~ "controls",
    enrichmentScore > 0 ~ "troponin",
    TRUE ~ "unknown"
  ))

# Enrichment troponin vs controls  in type 1:
result_GSEA_troponin_control_type_1 <- vroom::vroom(
  here::here("data/MD_GSEA_fiber_type/result_GSEA_troponin_control_type_1.csv")
) |>
  dplyr::select(!1) |>
  dplyr::mutate("enrichment" = dplyr::case_when(
    enrichmentScore < 0 ~ "controls_type_1",
    enrichmentScore > 0 ~ "troponin_type_1",
    TRUE ~ "unknown"
  ))

# Enrichment troponin vs controls  in type 2:
result_GSEA_troponin_control_type_2 <- vroom::vroom(
  here::here("data/MD_GSEA_fiber_type/result_GSEA_troponin_control_type_2.csv")
) |>
  dplyr::select(!1) |>
  dplyr::mutate("enrichment" = dplyr::case_when(
    enrichmentScore < 0 ~ "controls_type_2",
    enrichmentScore > 0 ~ "troponin_type_2",
    TRUE ~ "unknown"
  ))

combined_GSEA_results <- result_GSEA_actin_control |>
  dplyr::bind_rows(
    result_GSEA_actin_control_type_1,
    result_GSEA_actin_control_type_2,
    result_GSEA_troponin_control,
    result_GSEA_troponin_control_type_1,
    result_GSEA_troponin_control_type_2
  ) |>
  dplyr::select(enrichment, everything())

write.csv(
  combined_GSEA_results,
  here::here("doc/supplementary_tables/MD_GSEA_results_fiber_type.csv")
)
