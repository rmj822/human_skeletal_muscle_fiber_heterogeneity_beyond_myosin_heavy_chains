################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)

################################################################################################################################################
########################################################       FIGURE 5A    ###################################################################
################################################################################################################################################

# Load Transcriptomics results
transcriptomics <- read.csv(here::here("data/figure_5/slow_vs_fast.csv"), header = T) %>%
    drop_na(GENEID) %>%
    drop_na(padj)

# Extract significant genes per fiber factor
fast <- transcriptomics %>% dplyr::filter(padj < 0.05 & log2FoldChange < 0)
slow <- transcriptomics %>% dplyr::filter(padj < 0.05 & log2FoldChange > 0)


# Extract only non-coding DEGs
fast_noncoding <- fast %>% dplyr::filter(GENEBIOTYPE != "protein_coding") %>% dplyr::arrange(log2FoldChange)
slow_noncoding <- slow %>% dplyr::filter(GENEBIOTYPE != "protein_coding") %>% dplyr::arrange(log2FoldChange)

# Plot slow

slow_noncoding_plot <- slow_noncoding %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::slice(1:10)
plot_slow <- ggplot(slow_noncoding_plot,
                    aes(x = fct_reorder(GENE,log2FoldChange),
                        y = log2FoldChange,
                        alpha = log2FoldChange)
) +
    geom_bar(stat="identity",position = "dodge", fill = "#440154FF", color = "#440154FF", size = 0.35,alpha = 0.65) +
    coord_flip() +
    scale_alpha_continuous(name = "", range = c(0.75, 1)) +
    theme_classic() +
    ggtitle("Slow") +
    ylab("Fold change (Log2)") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme(
        text = ggplot2::element_text(face = "bold",size = 6, colour = "black"),
        plot.title = element_text(hjust = 0.5, face="bold", size=6),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"),
        legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        plot.margin = margin(0,0,0,0)
    )

# Plot fast
fast_noncoding_plot <- fast_noncoding %>% dplyr::arrange(log2FoldChange) %>% dplyr::slice(1:10)
plot_fast <- ggplot(fast_noncoding_plot,
                    aes(x = fct_reorder(GENE,-log2FoldChange),
                        y = -log2FoldChange,
                        alpha = -log2FoldChange)
) +
    geom_bar(stat="identity",position = "dodge", fill = "#618F70", color = "#618F70", size = 0.35, alpha = 0.65) +
    coord_flip() +
    scale_alpha_continuous(name = "", range = c(0.75, 1)) +
    theme_classic() +
    ggtitle("Fast") +
    ylab("Fold change (Log2)") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme(
        text = ggplot2::element_text(face = "bold",size = 6, colour = "black"),
        plot.title = element_text(hjust = 0.5, face="bold", size=6),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"),
        legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        plot.margin = margin(0,0,0,0)
    )

# Combine into one figure ------------------------------------------------------------------
Combined_plot_nc_type <- ggarrange(plot_slow,
                                   plot_fast,
                                   nrow=1)

annotate_figure(Combined_plot_nc_type, top = text_grob("Non-coding RNA", color = "black", face = "bold", size = 7))

ggplot2::ggsave(here::here("doc/figures/figure_5/figure_5A.png"),
                units = "mm",
                height = 55,
                width = 130)

################################################################################################################################################
########################################################       FIGURE 5C    ###################################################################
################################################################################################################################################


# Enter data from RNAscope quantification ---------------------------------
data <- data.frame(
    "gene" = rep(
        c("LINC01405",
          "RP11-255P5.3"
        ), 6
    ),
    "fiber_type" = c("slow",
                     "slow",
                     "slow",
                     "slow",
                     "slow",
                     "slow",
                     "fast",
                     "fast",
                     "fast",
                     "fast",
                     "fast",
                     "fast"),
    "replicate" = rep(c(
        "1",
        "1",
        "2",
        "2",
        "3",
        "3"
    ), 2),
    "dots_mm2" =
        c(
            28954.65,
            1547.0970,
            19968.54,
            1163.1080,
            14351.29,
            968.8703,
            856.0667,
            4431.203,
            799.0687,
            4666.297,
            304.7845,
            8430.472
        )
)

t_test_Linc <- t.test(
    x = data |>
        dplyr::filter(gene == "LINC01405") |>
        dplyr::filter(fiber_type == "fast") |>
        dplyr::pull(dots_mm2),
    y = data |>
        dplyr::filter(gene == "LINC01405") |>
        dplyr::filter(fiber_type == "slow") |>
        dplyr::pull(dots_mm2),
    var.equal = TRUE
)

t_test_RP <- t.test(
    x = data |>
        dplyr::filter(gene == "RP11-255P5.3") |>
        dplyr::filter(fiber_type == "fast") |>
        dplyr::pull(dots_mm2),
    y = data |>
        dplyr::filter(gene == "RP11-255P5.3") |>
        dplyr::filter(fiber_type == "slow") |>
        dplyr::pull(dots_mm2),
    var.equal = TRUE
)

box_plot_linc <- data |>
    dplyr::mutate(fiber_type = factor(fiber_type, levels = c("slow", "fast"))) |>
    dplyr::filter(gene == "LINC01405") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = fiber_type,
            y = dots_mm2
        )
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(fill = fiber_type), alpha = 0.85
    ) +
    ggplot2::geom_line(
        ggplot2::aes(group = replicate)
    ) +
    ggplot2::geom_point(ggplot2::aes(x = fiber_type, y = dots_mm2, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 31000, label = paste(round(
        t_test_Linc$p.value,
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    # ggplot2::scale_color_manual(values=c("#440154FF",
    #                                     "#5DC863FF")) +
    ggplot2::ggtitle("LINC01405") +
    ggplot2::ylab(bquote('Dots/mm' ^2)) +
    ggplot2::xlab("") +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        # axis.text.x = ggplot2::element_text(size = 7),
        axis.title.x = ggplot2::element_blank(),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=6.5),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5)
    ) +
    ggplot2::ylim(0, 32000) +
    ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "mm"))

box_plot_rp <- data |>
    dplyr::mutate(fiber_type = factor(fiber_type, levels = c("slow", "fast"))) |>
    dplyr::filter(gene == "RP11-255P5.3") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = fiber_type,
            y = dots_mm2
        )
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(fill = fiber_type), alpha = 0.85
    ) +
    ggplot2::geom_line(
        ggplot2::aes(group = replicate)
    ) +
    ggplot2::geom_point(ggplot2::aes(x = fiber_type, y = dots_mm2, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 10500, label = paste(round(
        t_test_RP$p.value,
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    # ggplot2::scale_color_manual(values=c("#440154FF",
    #                                     "#5DC863FF")) +
    ggplot2::ggtitle("RP11-255P5.3") +
    ggplot2::ylab(bquote('Dots/mm' ^2)) +
    ggplot2::xlab("") +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        # axis.text.x = ggplot2::element_text(size = 7),
        axis.title.x = ggplot2::element_blank(),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=6.5),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5)
    ) +
    ggplot2::ylim(0, 32000) +
    ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "mm"))


patchwork::wrap_plots(box_plot_linc,
                      box_plot_rp +
                          ggplot2::theme(
                              axis.line.y = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank(),
                              axis.title.y = ggplot2::element_blank()
                          ), guides = "collect") +
    ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "mm"))

ggplot2::ggsave(
    here::here("doc/figures/figure_5/figure_5C.png"),
    units = "mm",
    height = 50,
    width = 60
)

################################################################################################################################################
########################################################       FIGURE 5E    ###################################################################
################################################################################################################################################

# 1000 fiber proteome -----------------------------------------------------

data_raw <- vroom::vroom(here::here("data-raw/20230930_112903_SMF_directDIA_lncRNAs_DS_13092023_Report_proteins.tsv"))

metadata <- vroom::vroom(here::here("data/metadata_proteomics_seurat_clusters.csv")) |>
    dplyr::select(!...1)

contaminants <- data_raw |>
    dplyr::filter(grepl("KRT", PG.Genes)) |>
    dplyr::pull(PG.Genes)

data_raw <- data_raw |>
    dplyr::filter(!PG.Genes %in% contaminants)

old_column_names <- colnames(data_raw)

new_column_names <- gsub(
    pattern = ".*DIA_S",
    replacement = "",
    old_column_names
)

data_wrangled <- data_raw

new_column_names <- gsub(
    pattern = "_.*",
    replacement = "",
    new_column_names
)

BiocGenerics::colnames(data_wrangled) <- new_column_names

Gene.name <- data_wrangled$PG.Genes

selection_vector <- as.character(metadata$randomized_number)

data_wrangled <- data_wrangled |>
    dplyr::select(PG.ProteinGroups,
                  selection_vector) |>
    tibble::column_to_rownames("PG.ProteinGroups") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("randomized_number") |>
    dplyr::mutate(randomized_number = as.numeric(randomized_number))


data_wrangled <- data_wrangled |>
    dplyr::inner_join(
        metadata |>
            # dplyr::mutate(randomized_number = as.numeric(randomized_number)) |>
            dplyr::select(fiberID,
                          randomized_number)
    ) |>
    tibble::column_to_rownames("fiberID") |>
    dplyr::select(!randomized_number) |>
    t() |>
    as.data.frame()

# Filtering genes with less than 70% valid --------------------------------

#' Filtering missing values from rows
#'
#' @param .data dataset to filter
#' @param percentage_accepted_missing % of accepted missing values
#'
#' @return a dataframe
#' @export
#'
#' @examples
filtering_rows_Na <- function(.data, percentage_accepted_missing) {
    row_keep_vector <- .data |>
        is.na() |>
        rowSums()

    row_keep_vector <- row_keep_vector / ncol(.data)

    row_keep_vector <- row_keep_vector <= percentage_accepted_missing

    data_filtered <- .data |>
        tibble::add_column(row_keep_vector) |>
        dplyr::filter(row_keep_vector == T) |>
        dplyr::select(!starts_with("row"))

    return(data_filtered)
}

data_row_filtered <- data_wrangled |>
    filtering_rows_Na(0.3) |>
    tibble::rownames_to_column("protein_groups")

proteomics_data <- data_wrangled |>
    log2() |>
    as.data.frame()

# Doing pseudobulk on fast and slow fibers --------------------------------

metadata <- metadata |>
    dplyr::rename(
        "fiber_type_seurat" = "seurat_clusters"
    )


pseudobulk_maker <- function(.data, metadata, subject_id, grouping, colname) {
    selection_vector <- metadata |>
        dplyr::filter(fiber_type_seurat == grouping) |>
        dplyr::filter(subject == subject_id) |>
        dplyr::pull("fiberID")

    pseudobulk_median_calculation <- .data |>
        as.data.frame() |>
        dplyr::select(selection_vector) |>
        t() |>
        as.data.frame() |>
        dplyr::mutate(dplyr::across(.cols = everything(), median, na.rm = T)) |>
        unique() |>
        t()

    colnames(pseudobulk_median_calculation) <- colname
    return(pseudobulk_median_calculation)
}

data_pseudobulk <- data.frame(
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P1",
        grouping = "slow",
        colname = "P1_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P2",
        grouping = "slow",
        colname = "P2_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P3",
        grouping = "slow",
        colname = "P3_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P4",
        grouping = "slow",
        colname = "P4_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P5",
        grouping = "slow",
        colname = "P5_slow"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P1",
        grouping = "fast",
        colname = "P1_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P2",
        grouping = "fast",
        colname = "P2_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P3",
        grouping = "fast",
        colname = "P3_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P4",
        grouping = "fast",
        colname = "P4_fast"
    ),
    pseudobulk_maker(
        .data = as.data.frame(proteomics_data),
        metadata = metadata,
        subject_id = "P5",
        grouping = "fast",
        colname = "P5_fast"
    )
)

data_grouping_pseudobulk <- data.frame(
    "sample_ID" = c(
        "P1_slow",
        "P2_slow",
        "P3_slow",
        "P4_slow",
        "P5_slow",
        "P1_fast",
        "P2_fast",
        "P3_fast",
        "P4_fast",
        "P5_fast"
    ),
    "fiber_type" = c(
        "slow",
        "slow",
        "slow",
        "slow",
        "slow",
        "fast",
        "fast",
        "fast",
        "fast",
        "fast"
    )
)

data_limma <- data_pseudobulk |>
    limma::normalizeBetweenArrays(method = "quantile")

vector_factor_fiber_type <- factor(data_grouping_pseudobulk$fiber_type,
                                   levels = c(
                                       "slow",
                                       "fast"
                                   )
)

design_matrix <-
    model.matrix(
        ~ 0 + vector_factor_fiber_type,
        data_grouping_pseudobulk
    )

colnames(design_matrix) <- c(
    "slow",
    "fast"
)

fit <- limma::lmFit(
    data_limma,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "slow_vs_fast" = slow - fast,
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

genes_protein_groups <- data_raw |>
    dplyr::select(PG.ProteinGroups,
                  PG.Genes) |>
    dplyr::filter(PG.ProteinGroups %in% rownames(DE_results))

DE_results <- DE_results |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::inner_join(genes_protein_groups) |>
    dplyr::mutate(
        Genes = dplyr::case_when(
            grepl("ENS", PG.ProteinGroups) == TRUE ~ PG.ProteinGroups,
            TRUE ~ PG.Genes
        )
    ) |>
    dplyr::mutate(Genes = gsub(
        pattern = ";.*",
        replacement = "",
        Genes
    )) |>
    dplyr::select(!PG.ProteinGroups) |>
    dplyr::select(!PG.Genes)

data_volcano <- DE_results |>
    dplyr::mutate("signifficant" = dplyr::case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "enriched in slow",
        adj.P.Val < 0.05 & logFC < 0 ~ "enriched in fast",
        TRUE ~ "not signifficant"
    ))

DE_lnc <-
    data_volcano |>
    dplyr::filter(
        grepl("ENS", Genes)
    )

# slow LncRNA -------------------------------------------------------------

slow_lnc_proteomics <-
    data_volcano |>
    dplyr::filter(
        signifficant == "enriched in slow",
        grepl("ENS", Genes)
    ) |>
    dplyr::mutate(Genes = gsub(
        pattern = "_.*",
        replacement = "",
        Genes
    )) |>
    dplyr::pull(Genes) |>
    unique()

slow_lnc_transcriptomics <-
    vroom::vroom(here::here("data/figure_5/noncoding_DEG_slow_vs_fast/noncoding_DEG_slow.csv")) |>
    dplyr::pull(GENEID)


euler_object <- eulerr::euler(
    list(
        "LncRNA \nproteomics" = slow_lnc_proteomics,
        "LncRNA \ntranscriptomics" = slow_lnc_transcriptomics
    ),
    quantities = list(type = c("percent", "counts"))
)

euler_diagram <- ggplotify::as.ggplot(
    plot(
        euler_object,
        labels = list(cex = 0.75),
        quantities = list(
            type = c("counts"),
            font = 1,
            cex = 0.75
        ),
        fill = c(
            ggplot2::alpha("#440154FF", 0.85),
            ggplot2::alpha("#440154FF", 0.5)
        ),
        alpha = 0.65,
        edges = list(
            col = c(
                ggplot2::alpha("#440154FF", 1),
                ggplot2::alpha("#440154FF", 1)
            ),
            lex = 2
        ),
    )
)

euler_diagram

# Fast LncRNA -------------------------------------------------------------

fast_lnc_proteomics <-
    data_volcano |>
    dplyr::filter(
        signifficant == "enriched in fast",
        grepl("ENS", Genes)
    ) |>
    dplyr::mutate(Genes = gsub(
        pattern = "_.*",
        replacement = "",
        Genes
    )) |>
    dplyr::pull(Genes) |>
    unique()

fast_lnc_transcriptomics <-
    vroom::vroom(here::here("data/figure_5/noncoding_DEG_slow_vs_fast/noncoding_DEG_fast.csv")) |>
    dplyr::pull(GENEID)


euler_object <- eulerr::euler(
    list(
        "LncRNA \nproteomics" = fast_lnc_proteomics,
        "LncRNA \ntranscriptomics" = fast_lnc_transcriptomics
    ),
    quantities = list(type = c("percent", "counts"))
)

euler_diagram <- ggplotify::as.ggplot(
    plot(
        euler_object,
        labels = list(cex = 0.75),
        quantities = list(
            type = c("counts"),
            font = 1,
            cex = 0.75
        ),
        fill = c(
            ggplot2::alpha("#5DC863FF", 0.85),
            ggplot2::alpha("#5DC863FF", 0.5)
        ),
        alpha = 0.65,
        edges = list(
            col = c(
                ggplot2::alpha("#5DC863FF", 1),
                ggplot2::alpha("#5DC863FF", 1)
            ),
            lex = 2
        ),
    )
)

euler_diagram

# Identifying candidates:

lncs_complete <- unique(c(
    slow_lnc_proteomics,
    slow_lnc_transcriptomics,
    fast_lnc_proteomics,
    fast_lnc_transcriptomics
)) |>
    as.data.frame() |>
    dplyr::rename("lnc_list" = 1)


# For it to be a candidate, it needs to be shared between proteomics and transcriptomics
# And be common for the same fiber type exclusively:

candidates_slow <- slow_lnc_transcriptomics |>
    as.data.frame() |>
    dplyr::filter(slow_lnc_transcriptomics %in% slow_lnc_proteomics) |>
    dplyr::filter(!slow_lnc_transcriptomics %in% fast_lnc_proteomics) |>
    dplyr::rename("slow_candidates" = 1)


# Same for fast candidates:

candidates_fast <- fast_lnc_transcriptomics |>
    as.data.frame() |>
    dplyr::filter(fast_lnc_transcriptomics %in% fast_lnc_proteomics) |>
    dplyr::filter(!fast_lnc_transcriptomics %in% slow_lnc_proteomics) |>
    dplyr::rename("fast_candidates" = 1)


# Translate entries:

library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

# Make annotation dataframe with ENSEMBL ID and BIOTYPE

candidates_slow <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = candidates_slow$slow_candidates,
                                         columns = c("GENENAME","GENEID","GENEBIOTYPE"),
                                         keytype = "GENEID")

candidates_fast <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = candidates_fast$fast_candidates,
                                         columns = c("GENENAME","GENEID","GENEBIOTYPE"),
                                         keytype = "GENEID")

# Create plots of stringent lncRNAs ---------------------------------------

DE_lnc <-
    data_volcano |>
    dplyr::filter(grepl("ENS", Genes)) |>
    dplyr::mutate(
        transcript_id = gsub(pattern = "_.*",
                             replacement = "",
                             Genes),
        extra_info = gsub(pattern = ".*_",
                          replacement = "",
                          Genes)
    )

gene_ids <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                  keys = DE_lnc$transcript_id,
                                  columns = c("GENENAME","GENEID","GENEBIOTYPE"),
                                  keytype = "GENEID") |>
    dplyr::rename(transcript_id = GENEID)

selected_transcripts <- c(
    "RP11-296E23.1_ORF467:22475:22134",
    "LINC00598_ORF399:29559:29852",
    "RP13-143G15.4_ORF796:16308:16213",
    "LINC01405_ORF310:17438:17355"
)

DE_lnc <- DE_lnc |>
    dplyr::inner_join(
        gene_ids
    ) |>
    tidyr::unite(
        col = transcript_id,
        GENENAME, extra_info,
        sep = "_"
    ) |>
    dplyr::select(!Genes) |>
    dplyr::mutate(names = dplyr::case_when(
        transcript_id %in% selected_transcripts ~ transcript_id,
        TRUE ~ ""
    )) |>
    dplyr::mutate(names = gsub(pattern = ":.*",
                               replacement = "",
                               names)
    )

selected_transcripts <- c(
    "RP11-296E23.1_ORF467:22475:22134",
    "LINC00598_ORF399:29559:29852",
    "RP13-143G15.4_ORF796:16308:16213",
    "LINC01405_ORF310:17438:17355"
)

DE_lnc$signifficant <- factor(DE_lnc$signifficant, levels = c(
    "enriched in slow",
    "not signifficant",
    "enriched in fast"
))

DE_lnc |>
    ggplot2::ggplot(ggplot2::aes(
        x = logFC,
        y = -log10(P.Value),
        color = signifficant,
        label = transcript_id
    )) +
    ggplot2::geom_point(
        size = 0.5,
        alpha = ifelse(DE_lnc$signifficant == "not signifficant", 1, 1)
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "#440154FF",
            "grey",
            "#5DC863FF"
        ),
        name = ""
    ) +
    ggrepel::geom_label_repel(
        data = DE_lnc |>
            dplyr::filter(!names == ""),
        mapping = ggplot2::aes(
            x = logFC,
            y = -log10(P.Value),
            fill = signifficant,
            label = names
        ),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#efedf5",
        "#e5f5e0"
    )) +
    ggplot2::ggtitle("Slow Vs Fast microproteins") +
    ggplot2::theme_classic() +
    ggplot2::theme(
        text = ggplot2::element_text(
            face = "bold",
            size = 12,
            colour = "black"
        ),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::xlab("log2FC (Slow - Fast)") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7),
        legend.position = "none"
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_5_S1/figure_5_S1C.pdf"),
#     height = 60,
#     width = 60,
#     units = "mm"
# )

# Boxplots Lnc01405 -------------------------------------------------------

data_boxplots <- data_limma |>
    as.data.frame() |>
    tibble::rownames_to_column("protein_groups") |>
    dplyr::filter(grepl("ENSG00000185847", protein_groups))

data_boxplots <- data_boxplots |>
    dplyr::mutate(gene_id = gsub(pattern = "_.*",
                                 replacement = "",
                                 data_boxplots$protein_groups
    ))

gene_ids <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                  keys = data_boxplots$gene_id,
                                  columns = c("GENENAME","GENEID","GENEBIOTYPE"),
                                  keytype = "GENEID") |>
    dplyr::rename(gene_id = GENEID)

data_boxplots <- data_boxplots |>
    dplyr::inner_join(gene_ids) |>
    dplyr::mutate(
        protein_groups = gsub(
            pattern = ";.*",
            replacement = "",
            protein_groups
        )
    ) |>
    dplyr::mutate(protein_groups = gsub(
        pattern = ".*_",
        replacement = "",
        protein_groups
    )) |>
    tidyr::unite(col = "transcript_id",
                 GENENAME, protein_groups,
                 sep = "_") |>
    dplyr::select(!GENEBIOTYPE) |>
    dplyr::select(!gene_id)

data_boxplots <- data_boxplots |>
    tidyr::pivot_longer(
        cols = 2:11,
        values_to = "LFQ_intensities",
        names_to = "sample_id"
    ) |>
    dplyr::mutate(fiber_type =
                      dplyr::case_when(
                          grepl("slow", sample_id) == TRUE ~ "slow",
                          TRUE ~ "fast"
                      )
    )


data_boxplots$fiber_type <- factor(data_boxplots$fiber_type, levels = c("slow", "fast"))

boxplot_LNC_fiber_type <- function(.data, LNC_name, height) {


    # data_boxplots <- data_boxplots |>
    # dplyr::mutate(
    #     sample_id = gsub(
    #         pattern = "_.*",
    #         replacement = "",
    #         sample_id
    #     )
    # )

    output <- .data |>
        dplyr::filter(transcript_id == LNC_name) |>
        ggplot2::ggplot(
            ggplot2::aes(
                x = fiber_type,
                y = LFQ_intensities
            )
        ) +
        ggplot2::geom_boxplot(
            ggplot2::aes(fill = fiber_type), alpha = 0.85
        ) +
        ggplot2::geom_point(ggplot2::aes(x = fiber_type, y = LFQ_intensities, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
        ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = height, label = paste("P.val = ", round(
            3.1,
            4
        )),
        label.size = 2) +
        ggplot2::scale_fill_manual(values=c("#440154FF",
                                            "#5DC863FF")) +
        ggplot2::ggtitle(LNC_name) +
        ggplot2::labs(
            y = "LFQ intensities",
            x = ""
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(size = 7),
            text = ggplot2::element_text(face="bold",
                                         colour="black",
                                         size=7),
            plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
            plot.margin = grid::unit(c(0,0,0,0), "mm")
        )

    return(output)
}

# LINC01405_ORF310:17438:17355

# boxplot_LNC_fiber_type(data_boxplots,
#                        LNC_name = "LINC01405_ORF408:17441:17358",
#                        height = 4.6)

t_test_result <- t.test(
    x = data_boxplots |>
        dplyr::mutate(
            sample_id = gsub(
                pattern = "_.*",
                replacement = "",
                sample_id
            )
        ) |>
        dplyr::filter(fiber_type == "slow") |>
        dplyr::filter(transcript_id == "LINC01405_ORF310:17438:17355") |>
        dplyr::pull(LFQ_intensities, name = sample_id),
    y = data_boxplots |>
        dplyr::mutate(
            sample_id = gsub(
                pattern = "_.*",
                replacement = "",
                sample_id
            )
        ) |>
        dplyr::filter(fiber_type == "fast") |>
        dplyr::filter(transcript_id == "LINC01405_ORF310:17438:17355") |>
        dplyr::pull(LFQ_intensities, name = sample_id)
)

# data_boxplots <- data_boxplots |>
#     dplyr::mutate(
#         sample_id = gsub(
#             pattern = "_.*",
#             replacement = "",
#             sample_id
#         )
#     )

data_boxplots |>
    dplyr::mutate(subject =
                      gsub(
                          pattern = "_.*",
                          replacement = "",
                          sample_id
                      )) |>
    dplyr::filter(transcript_id == "LINC01405_ORF310:17438:17355") |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = fiber_type,
            y = LFQ_intensities
        )
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(fill = fiber_type), alpha = 0.85
    ) +
    ggplot2::geom_line(
        ggplot2::aes(group = subject)
    ) +
    ggplot2::geom_point(ggplot2::aes(x = fiber_type, y = LFQ_intensities, fill = fiber_type), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 4.75, label = paste(round(
        DE_results |>
            dplyr::filter(grepl("ORF310:17438:17355",Genes)) |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("LINC01405_ORF310:17438:17355") +
    ggplot2::labs(
        y = "LFQ intensities",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 6),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=6),
        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
        plot.margin = grid::unit(c(0,0,0,0), "mm"))

ggplot2::ggsave(here::here("doc/figures/figure_5/figure_5E.pdf"),
                units = "mm",
                width = 60,
                height = 50)
