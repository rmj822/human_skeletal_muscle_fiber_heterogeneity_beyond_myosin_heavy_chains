################################################################################################################################################
#################################################     Panel B  ##############################################################
################################################################################################################################################

# Loading data ------------------------------------------------------------

# Also simplifying sample names

bulk_data <- vroom::vroom(here::here("data-raw/nemaline_myopathy_bulk.pg_matrix.tsv")) |>
    dplyr::filter(!duplicated(Genes)) |>
    dplyr::filter(!Genes == "") |>
    tibble::column_to_rownames("Genes") |>
    dplyr::select(dplyr::contains("Roger")) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_name") |>
    dplyr::mutate(sample_name = gsub(
        pattern = ".*96_",
        replacement = "",
        sample_name
    )) |>
    dplyr::mutate(sample_name = gsub(
        pattern = "_S.*",
        replacement = "",
        sample_name
    )) |>
    tibble::column_to_rownames("sample_name") |>
    t() |>
    as.data.frame()

# Create metadata:

metadata <- data.frame(
    "sample_id" = c(
        "control_1",
        "control_2",
        "control_3",
        "troponin_1",
        "troponin_2",
        "troponin_3"
    ),
    "grouping" = c(
        "control",
        "control",
        "control",
        "troponin",
        "troponin",
        "troponin"
    ),
    "subject" = c(
        "control_1",
        "control_2",
        "control_3",
        "troponin_1",
        "troponin_2",
        "troponin_3"
    )
)

# Data processing ---------------------------------------------------------

data_processed <- bulk_data |>
    log2() |>
    PhosR::selectGrps(
        grps = metadata$grouping,
        percent = 0.7
    ) |>
    limma::normalizeBetweenArrays(method = "quantile") |>
    as.data.frame()


# PCA ---------------------------------------------------------------------

data_pca <- data_processed |>
    PhosR::tImpute() |>
    t() |>
    prcomp()

data_pca$x |>
    as.data.frame() |>
    dplyr::select(PC1, PC2) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::inner_join(metadata) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = grouping
        )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::geom_point(
        size = 1.5,
        alpha = 0.65
    ) +
    ggplot2::scale_color_manual(
        "Grouping",
        values = c(
            "#969594",
            "#d662c4"
        )
    ) +
    ggrepel::geom_label_repel(
        data = data_pca$x |>
            as.data.frame() |>
            dplyr::select(PC1, PC2) |>
            tibble::rownames_to_column("sample_id") |>
            dplyr::inner_join(metadata),
        mapping = ggplot2::aes(
            x = PC1,
            y = PC2,
            fill = grouping,
            label = sample_id
        ),
        color = "black",
        size = 1.8,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 10
    ) +
    ggplot2::scale_fill_manual(
        "Grouping",
        values = c(
            ggplot2::alpha("#969594", alpha = 0.35),
            ggplot2::alpha("#d662c4", alpha = 0.35)
        )
    ) +
    ggplot2::xlab("PC 1 (53.5%)") +
    ggplot2::ylab("PC 2 (25.2%)") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6),
        legend.position = "none"
    )

ggplot2::ggsave(here::here("doc/figures/figure_6_S4/figure_6_S4B.png"),
                units = "mm",
                height = 60,
                width = 60)

################################################################################################################################################
#################################################     Panel C  ##############################################################
################################################################################################################################################

# Limma analysis: ---------------------------------------------------------

design_matrix <- model.matrix(
    ~ 0 + metadata$grouping,
    data_processed
)

colnames(design_matrix) <- c(
    "control",
    "troponin"
)

fit <- limma::lmFit(
    data_processed,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "troponin-control" = troponin - control,
    levels = design_matrix
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp <- limma::eBayes(tmp)

DE_results <- limma::topTable(tmp,
                              sort.by = "P",
                              n = Inf
)

DE_results <- DE_results |>
    dplyr::mutate(xiao = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

DE_results <- DE_results |>
    dplyr::mutate("significant" = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC > 0 ~ "Upregulated (π < 0.05)",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated (FDR < 0.05)",
        xiao <= 0.05 & logFC < 0 ~ "Downregulated (π < 0.05)",
        TRUE ~ "not significant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate("names" = dplyr::case_when(
        Genes %in% c(
            "RCC1",
            "ASPN",
            "BGN",
            "COL6A2",
            "SOD3",
            "ANXA1",
            "ANXA5",
            "TNNT1",
            "CKMT2",
            "COX5A",
            "ATP5PF",
            "NDUFA5"
        ) ~ Genes,
        TRUE ~ ""
    ))

DE_results |>
    ggplot2::ggplot(ggplot2::aes(
        x = logFC,
        y = -log10(P.Value),
        color = significant,
        names = Genes
    )) +
    ggplot2::geom_point(
        size = 0.5,
        shape = ifelse(
            DE_results$significant == "Upregulated (π < 0.05)",
            1,
            ifelse(DE_results$significant == "Downregulated (π < 0.05)", 1, 16)
        ),
        alpha = ifelse(DE_results$significant == "not significant", 0.3, 0.8)
    ) +
    ggplot2::scale_color_manual(values = c("#969594",
                                           "#969594",
                                           "lightgrey",
                                           "#d662c4",
                                           "#d662c4")) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("TNNT1-NM vs control") +
    ggplot2::xlab("log2FC (TNNT1-NM - control)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(
            size = 5,
            colour = "black",
            face = "bold"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           face = "bold")
    ) +
    ggrepel::geom_label_repel(
        data = DE_results |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = significant,
            label = names
        ),
        color = "black",
        size = 1.5,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 40
    ) +
    ggplot2::scale_fill_manual(values = c("#d5d4d4",
                                          "#eab0e1"))

ggplot2::ggsave(here::here("doc/figures/figure_6_S4/figure_6_S4C.png"),
                units = "mm",
                height = 60,
                width = 60)


################################################################################################################################################
#################################################     Panel D  ##############################################################
################################################################################################################################################

# Ranked GSEA -------------------------------------------------------------

library(org.Hs.eg.db)

GSEA_proteins <- DE_results |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::filter(!logFC == "") |>
    dplyr::pull(logFC, name = Genes)

GSEA_results <- clusterProfiler::gseGO(
    GSEA_proteins,
    keyType = "SYMBOL",
    ont = "ALL",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
)

GSEA_main_effect_training <- clusterProfiler::simplify(
    GSEA_results
)

dotplot_title <- c(
    `activated` = "Upregulated",
    `suppressed` = "Downregulated"
)

results <- GSEA_results@result

results |>
    dplyr::mutate(
        GeneRatio =stringr::str_count(core_enrichment, "\\w+")/setSize,
        regulated = dplyr::case_when(
            NES > 0 ~ "Upregulated \nin TNNT1-NM",
            TRUE ~ "Downregulated \nin TNNT1-NM"
        )
    ) |>
    dplyr::filter(Description %in% c(
        "collagen-containing extracellular matrix",
        "extracellular matrix",
        "cell-matrix adhesion",
        "inflammatory response",
        "mitochondrion",
        "oxidative phosphorylation",
        "mitochondrial ATP synthesis coupled electron transport"
    )) |>
    dplyr::mutate(Description = dplyr::case_when(
        Description == "mitochondrial ATP synthesis coupled electron transport" ~ "mitochondrial ATP synthesis \ncoupled electron transport",
        TRUE ~ Description
    )) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = GeneRatio,
            y = Description,
            color = p.adjust,
            size = abs(NES)
        )
    ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~regulated) +
    ggplot2::scale_color_viridis_c("Adj.P.Val", option = "plasma") +
    ggplot2::scale_size_continuous("Normalized\nenrichment\nscore",
                                   limits = c(2, 3.9),
                                   breaks = c(2, 3, 4)) +
    ggplot2::theme(
        text = ggplot2::element_text(size = 6),
        legend.key.size = ggplot2::unit(3, "mm"),
        strip.text = ggplot2::element_text(size = 6, face = "bold"),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5, face = "bold"),
        axis.title.y = ggplot2::element_blank()
    )

ggplot2::ggsave(here::here("doc/figures/figure_6_S4/figure_6_S4D.png"),
                units = "mm",
                height = 60,
                width = 120)
