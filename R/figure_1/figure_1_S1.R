################################################################################################################################################
################################################       Panel B      ########################################################################
################################################################################################################################################

load(here::here("data/figure_2/filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest.Rdata"))

transcriptomics_data <- filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest@meta.data |>
    as.data.frame()

transcriptomics_data|>
    ggplot2::ggplot(
        ggplot2::aes(
            x = subject,
            y = nGene,
            fill = subject
        )
    ) +
    ggplot2::geom_violin(alpha = 0.75) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual(values = c(
        "#070B08",
        "#162119",
        "#26372B",
        "#354D3C",
        "#43624D",
        "#52795F",
        "#618F70",
        "#73A182",
        "#88AF94",
        "#9EBEA8",
        "#B5CCBB",
        "#C9DBCF",
        "#DFEAE2",
        "#F5F8F5"
    ),
    guide = "none") +
    ggplot2::ylim(1000, 7500) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1000 fiber transcriptome",
                  y = "Number of quantified genes") +
    ggplot2::theme(text = ggplot2::element_text(size = 9.5)) +
    ggplot2::theme(legend.position = "none") +
    scale_x_discrete(labels=c("1" = "T1",
                              "2" = "T2",
                              "3" = "T3",
                              "4" = "T4",
                              "5" = "T5",
                              "6" = "T6",
                              "7" = "T7",
                              "8" = "T8",
                              "9" = "T9",
                              "10" = "T10",
                              "11" = "T11",
                              "12" = "T12",
                              "13" = "T13",
                              "14" = "T14"))

ggsave(filename = "doc/figures/figure_1_S1/figure_1_S1B.png", width = 128, height = 60, units="mm")


colnames(transcriptomics_data) <- "number_of_proteins"

# Loading proteomics data -------------------------------------------------

proteomics_data <- vroom::vroom(
    here::here("data/data_proteomics_filtered.csv")
) |>
    dplyr::select(!1) |>
    as.data.frame() |>
    dplyr::filter(
        !Gene.name == "",
        !duplicated(Gene.name)
    ) |>
    tibble::column_to_rownames("Gene.name")

metadata <- vroom::vroom(
    here::here("data/metadata_proteomics_fiber_type.csv")
)

heterofiber_number_quantified_proteins <- colSums(!is.na(proteomics_data)) |>
    as.data.frame()

colnames(heterofiber_number_quantified_proteins) <- "number_of_proteins"

heterofiber_number_quantified_proteins <- heterofiber_number_quantified_proteins |>
    dplyr::mutate(study = rep("heterofiber",
                              ncol(proteomics_data))) |>
    tibble::rownames_to_column("fiberID")

heterofiber_number_quantified_proteins <- metadata |>
    dplyr::filter(fiberID %in% heterofiber_number_quantified_proteins$fiberID) |>
    dplyr::filter(!duplicated(fiberID)) |>
    dplyr::select(fiberID,
                  subject) |>
    dplyr::inner_join(heterofiber_number_quantified_proteins) |>
    tibble::column_to_rownames("fiberID")

data_plot <- heterofiber_number_quantified_proteins |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::filter(!number_of_proteins < 1000) |>
    dplyr::mutate(project = dplyr::case_when(
        study == "heterofiber" ~ "proteomics",
        TRUE ~ "muscle_disease"
    ))

################################################################################################################################################
########################################################      Panel C   ############################################################################
################################################################################################################################################


# Violins 1000 fiber proteome ---------------------------------------------
data_plot_heterofiber <- data_plot |>
    dplyr::filter(project == "proteomics")

data_plot_heterofiber$subject <- factor(data_plot_heterofiber$subject,
                                        levels = c("P1",
                                                   "P2",
                                                   "P3",
                                                   "P4",
                                                   "P5"))

data_plot_heterofiber |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = subject,
            y = number_of_proteins,
            fill = subject
        )
    ) +
    ggplot2::geom_violin(alpha = 0.75) +
    ggplot2::geom_point(position = ggplot2::position_jitter(seed = 1, width = 0.3),
                        size = 0.15) +
    ggplot2::scale_fill_manual(values = c(
        "#08519c",
        "#3182bd",
        "#6baed6",
        "#bdd7e7",
        "#eff3ff"
    ),
    guide = "none") +
    ggplot2::ylim(1000, 2500) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1000 fiber proteome",
                  y = "Number of quantified proteins") +
    ggplot2::theme(text = ggplot2::element_text(size = 9.5)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_discrete(labels = c("P1",
                                         "P2",
                                         "P3",
                                         "P4",
                                         "P5"))

ggplot2::ggsave(here::here("doc/figures/figure_1_S1/figure_1_S1C.png"),
                height = 60,
                width = 60,
                units = "mm")

################################################################################################################################################
########################################################      Panel D   ############################################################################
################################################################################################################################################

# CV ----------------------------------------------------------------------

proteomics_raw <- utils::read.table(
    here::here(
        "data-raw/loading_controls_proteomics.txt"
    ),
    header = TRUE,
    sep = "\t"
)

contaminants <- proteomics_raw |>
    dplyr::distinct(Gene.name, .keep_all = TRUE) |>
    dplyr::filter(!Gene.name == "") |>
    tibble::column_to_rownames("Gene.name") |>
    t() |>
    as.data.frame() |>
    dplyr::select(dplyr::contains("KRT"))

contaminants <- colnames(contaminants)

proteomics_raw <- proteomics_raw |>
    dplyr::filter(!Gene.name %in% contaminants)

Genes <- proteomics_raw |>
    dplyr::pull(Gene.name)


# Renaming columns and filtering for technical controls -------------------

old_column_names <- colnames(proteomics_raw)

new_column_names <- gsub(
    pattern = ".*DIA_",
    replacement = "",
    old_column_names
)

new_column_names <- gsub(
    pattern = "_.*",
    replacement = "",
    new_column_names
)

colnames(proteomics_raw) <- new_column_names

proteomics_controls <- proteomics_raw |>
    dplyr::select(dplyr::starts_with("C")) |>
    dplyr::bind_cols("Gene.name" = Genes) |>
    dplyr::select(!C11) |>
    dplyr::filter(!Gene.name == "",
                  !duplicated(Gene.name)) |>
    tibble::column_to_rownames("Gene.name")


# Calculating mean of intensities and CVs ---------------------------------


proteomics_controls <- proteomics_controls |>
    dplyr::mutate(coefficient_variation = apply(proteomics_controls,
                                                1,
                                                FUN = function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)) |>
    dplyr::mutate(mean_LFQ = rowMeans(proteomics_controls, na.rm = TRUE))


# Loading filtered dataset and annotating CV dataframe --------------------

proteomics_filtered <- vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1)

data_plot <- proteomics_controls |>
    tibble::rownames_to_column("Gene.name") |>
    dplyr::mutate(filtered = dplyr::case_when(
        Gene.name %in% proteomics_filtered$Gene.name ~ "kept",
        TRUE ~ "dropped"
    )) |>
    dplyr::mutate(perc_CVs = dplyr::case_when(
        coefficient_variation <= 30 & coefficient_variation >= 20 ~ "under_30",
        coefficient_variation <= 20 & coefficient_variation >= 10~ "under_20",
        coefficient_variation <= 10 ~ "under_10",
        TRUE ~ "high_perc"
    ))

table(data_plot$perc_CVs)


data_plot |>
    ggplot2::ggplot(
        ggplot2::aes(x = log2(mean_LFQ),
                     y = coefficient_variation,
                     color = perc_CVs)
    ) +
    ggplot2::geom_point(size = 0.5,
                        alpha = 0.5,
                        na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(
        values = c("#bdbdbd",
                   "#08519c",
                   "#3182bd",
                   "#6baed6")
    ) +
    ggplot2::geom_hline(
        yintercept = 10,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 24.5,
        y = 5,
        fontface = 2,
        label = "182",
        size = 2.25
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 23.5,
        y = 15,
        fontface = 2,
        label = "1153",
        size = 2.25
    ) +
    ggplot2::annotate(
        geom = "text",
        x = 22.5,
        y = 25,
        fontface = 2,
        label = "2048",
        size = 2.25
    ) +
    ggplot2::geom_hline(
        yintercept = 20,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = 30,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::ylab("Mean coefficient of variation (%)") +
    ggplot2::xlab("LFQ intensities (Log2)") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size = 8),
        legend.position = "none",
    )
# ggplot2::coord_cartesian(xlim = c(7.55, 23), clip = "off")

ggplot2::ggsave(filename = here::here("doc/figures/figure_1_S1/figure_1_S1D.png"),
                height = 60,
                width = 60,
                units = "mm")

################################################################################################################################################
########################################################      Panel E   ############################################################################
################################################################################################################################################

counts <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")

count_plot <- counts

# Calculate relative expression for each gene for each fiber
count_plot@x <- (count_plot@x/rep.int(colSums(count_plot), diff(count_plot@p))) * 100

# Take average across fibers
mean_expression_filtered <- rowMeans(count_plot)
mean_expression_filtered <- as.data.frame(mean_expression_filtered)

# Rename mean expression column
mean_expression_filtered <- mean_expression_filtered %>% dplyr::rename(mean_expression = mean_expression_filtered)

# Add column with gene names
mean_expression_filtered$gene <- rownames(mean_expression_filtered)

# Create rank order
mean_expression_filtered <- mean_expression_filtered %>%
    dplyr::arrange(desc(mean_expression))
mean_expression_filtered$order <- seq_len(nrow(mean_expression_filtered))


ggplot() +

    # Add all genes
    geom_point(data = mean_expression_filtered %>% filter(mean_expression != 0), aes(x=order, y=log(mean_expression,10)), colour="#B7DFB3", size=0.25) +

    # Add rectangle for filtered genes
    annotate("rect", xmin = -Inf, xmax = 7418, ymin = -Inf, ymax = Inf,
             alpha = .15, fill = "#4F7A5D") +

    # Add only filtered genes
    geom_point(data = mean_expression_filtered,
               aes(x=order, y=log(mean_expression,10)),
               colour="#4F7A5D",
               size=0.25) +

    # Add vertical line as cut-off of low abundant gene filtering
    geom_vline(xintercept = nrow(mean_expression_filtered)) +

    # Add horizontal lines to indicate % instead of log scale
    geom_hline(yintercept = log(20,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(5,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(1,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(0.1,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(0.01,10), linetype="dashed", linewidth=0.2) +

    # Add text for % expression
    annotate("text", x=26000, y=log(40,10), label= "20%", colour="black", fontface=2, size=2) +
    annotate("text", x=26000, y=log(9.4,10), label= "5%", colour="black", fontface=2, size=2) +
    annotate("text", x=26000, y=log(2.1,10), label= "1%", colour="black", fontface=2, size=2) +
    annotate("text", x=26000, y=log(0.21,10), label= "0.1%", colour="black", fontface=2, size=2) +
    annotate("text", x=26000, y=log(0.02,10), label= "0.01%", colour="black", fontface=2, size=2) +

    # Add labels interesting muscle cell markers
    geom_label_repel(data = mean_expression_filtered %>% dplyr::filter(gene == "ACTA1" | gene == "MYH2" | gene == "MYH7" | gene == "TNNT1" | gene == "TNNT3"),
                     mapping = aes(order, log(mean_expression,10), label = gene),
                     size = 1.5, max.overlaps = Inf, label.padding=0.1, min.segment.length=0.1, segment.size=0.2, max.time = 10, fill = "white", force=80) +

    # Add dotsMuscle cell markers
    geom_point(data = mean_expression_filtered %>% dplyr::filter(gene == "ACTA1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_filtered %>% dplyr::filter(gene == "MYH2"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_filtered %>% dplyr::filter(gene == "MYH7"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_filtered %>% dplyr::filter(gene == "TNNT1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_filtered %>% dplyr::filter(gene == "TNNT3"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Change design
    ylab("% total counts, 10log") +
    xlab("Gene rank") +
    ggtitle("Transcriptomics") +
    theme_classic() +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8)
    ) +
    scale_y_continuous(limits = c(-6.5,3), expand = c(0, 0))

ggsave(plot_expression_all_without_nonmuscle, filename = "doc/figures/figure_1_S1/figure_S1E.png", width = 60, height = 60, units="mm")



################################################################################################################################################
########################################################      Panel F   ############################################################################
################################################################################################################################################


# Dynamic range plot ------------------------------------------------------

sum_of_intensities <- colSums(proteomics_filtered |>
                                  tibble::column_to_rownames("Gene.name"),
                              na.rm = TRUE)

rel_abundance <- proteomics_filtered |>
    tibble::column_to_rownames("Gene.name")|>
    t() |>
    as.data.frame()

rel_abundance <- rel_abundance/sum_of_intensities * 100

mean_rel_abundance <- rel_abundance |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        mean,
        na.rm = TRUE
    )) |>
    dplyr::slice_head(n = 1) |>
    t() |>
    as.data.frame() |>
    log10()

colnames(mean_rel_abundance) <- "log10_rel_abundance"

mean_rel_abundance <- mean_rel_abundance |>
    dplyr::arrange(desc(log10_rel_abundance)) |>
    dplyr::mutate(
        order = 1:nrow(proteomics_filtered)
    )

Gene_names <- rownames(mean_rel_abundance)

mean_expression_all <- mean_rel_abundance |>
    tibble::rownames_to_column("gene")

ggplot2::ggplot() +

    # Add all genes
    ggplot2::geom_point(data = mean_expression_all,
                        ggplot2::aes(
                            x=order,
                            y=log10_rel_abundance),
                        colour = "#045a8d",
                        size = 0.25,
                        alpha = 0.5) +

    # Add horizontal lines to indicate % instead of log scale
    ggplot2::geom_hline(yintercept = log(20,10), linetype="dashed", linewidth=0.2) +
    ggplot2::geom_hline(yintercept = log(5,10), linetype="dashed", linewidth=0.2) +
    ggplot2::geom_hline(yintercept = log(1,10), linetype="dashed", linewidth=0.2) +
    ggplot2::geom_hline(yintercept = log(0.1,10), linetype="dashed", linewidth=0.2) +
    ggplot2::geom_hline(yintercept = log(0.01,10), linetype="dashed", linewidth=0.2) +

    # Add text for % expression
    ggplot2::annotate("text", x=2900, y=log(28,10), label= "20%", colour="black", fontface=2, size=2) +
    ggplot2::annotate("text", x=2900, y=log(7,10), label= "5%", colour="black", fontface=2, size=2) +
    ggplot2::annotate("text", x=2900, y=log(1.4,10), label= "1%", colour="black", fontface=2, size=2) +
    ggplot2::annotate("text", x=2900, y=log(0.14,10), label= "0.1%", colour="black", fontface=2, size=2) +
    ggplot2::annotate("text", x=2900, y=log(0.014,10), label= "0.01%", colour="black", fontface=2, size=2) +

    # Add custom labels:
    ggrepel::geom_label_repel(data = mean_expression_all |> dplyr::filter(gene %in% c("MYH7",
                                                                                   "MYH2",
                                                                                   "ACTA1",
                                                                                   "TNNT1",
                                                                                   "TNNT3")),
                     mapping = ggplot2::aes(order, log10_rel_abundance, label = gene),
                     size = 1.8, label.padding=0.1, max.overlaps = Inf, min.segment.length=0.1, segment.size=0.2, force = 10) +

# Change design
    ggplot2::ylab("% total intensities, 10log") +
    ggplot2::xlab("Protein rank") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Proteomics") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size = 6),
        plot.title = ggplot2::element_text(face = "bold", color = "black", size = 8, hjust = 0.5),
        strip.text = ggplot2::element_text(colour = "white"),
        strip.background = ggplot2::element_rect(fill="black"),
        legend.position = "none",
    )

ggplot2::ggsave(here::here("doc/figures/figure_1_S1/figure_1_S1F.png"),
       units = "mm",
       height = 60,
       width = 60)

################################################################################################################################################
########################################################      Panel G   ############################################################################
################################################################################################################################################

library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(rtracklayer)
library(GenomicFeatures)
library(ggrepel)
library(VennDiagram)

transcriptomics <- GetAssayData(object = filtered_normalized_fibertype_clustered_seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")

# Transcriptomics ---------------------------------------------------------------
transcriptomics_plot <- transcriptomics

# Calculate relative expression for each gene for each fiber
transcriptomics_plot@x <- (transcriptomics_plot@x/rep.int(colSums(transcriptomics_plot), diff(transcriptomics_plot@p))) * 100

# Take average across fibers
mean_expression_transcriptomics<- rowMeans(transcriptomics_plot)
mean_expression_transcriptomics <- as.data.frame(mean_expression_transcriptomics)

# Rename mean expression column
mean_expression_transcriptomics <- mean_expression_transcriptomics %>% dplyr::rename(mean_expression = mean_expression_transcriptomics)

# Add column with gene names
mean_expression_transcriptomics$gene <- rownames(mean_expression_transcriptomics)

# Create rank order
mean_expression_transcriptomics <- mean_expression_transcriptomics %>%
    dplyr::arrange(desc(mean_expression))
mean_expression_transcriptomics$order <- seq_len(nrow(mean_expression_transcriptomics))

# Proteomics ---------------------------------------------------------------
proteomics_plot <- proteomics

# Transpose df
proteomics_plot_transpose <- t(proteomics_plot)

# Change NA to 0
proteomics_plot_transpose[is.na(proteomics_plot_transpose)] <- 0

# Calculate relative abundance for each gene, separately for each fiber
proteomics_plot_transpose_relabundance <- t(apply(proteomics_plot_transpose, 1, function(x) (x/sum(x)*100)))

# Re-transpose relative abundances
proteomics_plot_relabundance <- t(proteomics_plot_transpose_relabundance)

# Take average across fibers
mean_expression_proteomics <- rowMeans(proteomics_plot_relabundance)
mean_expression_proteomics <- as.data.frame(mean_expression_proteomics)

# Rename mean expression column
mean_expression_proteomics <- mean_expression_proteomics %>% dplyr::rename(mean_expression = mean_expression_proteomics)

# Add column with gene names
mean_expression_proteomics$protein <- rownames(mean_expression_proteomics)

# Create rank order
mean_expression_proteomics <- mean_expression_proteomics %>%
    dplyr::arrange(desc(mean_expression))
mean_expression_proteomics$order <- seq_len(nrow(mean_expression_proteomics))

# Tx
transcriptomics_raw <- transcriptomics

mean_counts_transcriptomics <- rowMeans(transcriptomics_raw)
mean_counts_transcriptomics <- as.data.frame(mean_counts_transcriptomics)
mean_counts_transcriptomics$gene <- rownames(mean_counts_transcriptomics)


mean_counts_transcriptomics <- mean_counts_transcriptomics %>%
    dplyr::arrange(desc(mean_counts_transcriptomics))
mean_counts_transcriptomics$order_T <- seq_len(nrow(mean_counts_transcriptomics))

# Px
mean_counts_proteomics <- rowMeans(proteomics_plot, na.rm=TRUE)
mean_counts_proteomics <- as.data.frame(mean_counts_proteomics)
mean_counts_proteomics <- na.omit(mean_counts_proteomics)
mean_counts_proteomics$protein <- rownames(mean_counts_proteomics)

mean_counts_proteomics <- mean_counts_proteomics %>%
    dplyr::arrange(desc(mean_counts_proteomics))
mean_counts_proteomics$order_P <- seq_len(nrow(mean_counts_proteomics))

# Filter transcriptomics and proteomics to only keep genes detected in both datasets
mean_counts_transcriptomics_corr <- subset(mean_counts_transcriptomics, gene %in% mean_counts_proteomics$protein)
mean_counts_preoteomics_corr <- subset(mean_counts_proteomics, protein %in% mean_counts_transcriptomics$gene)

# Combine dataframes
combined_raw <- merge(mean_counts_transcriptomics_corr, mean_counts_preoteomics_corr, by.x="gene", by.y="protein")


################################################################################################################################################
####################################     Venn diagram detected genes transcriptomics vs proteomics      ###################################
################################################################################################################################################

Venn_detected_genes <- venn.diagram(
    # General
    filename=NULL,
    disable.logging=T,
    x = list(
        mean_expression_transcriptomics %>% dplyr::select(gene) %>% unlist(use.names=F),
        mean_expression_proteomics %>% dplyr::select(protein) %>% unlist(use.names=F)
    ),
    category.names = c("Transcriptomics" , "Proteomics"),
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.cex = 0.5,

    # Circles
    lwd = 2,
    col=c("#4F7A5D", "#045a8d"),
    fill = c(alpha("#4F7A5D",0.3), alpha('#045a8d',0.3)),

    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",

    # Names
    cat.cex = 0.75,
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-25,15),
    cat.dist = c(0.05, 0.1),
    cat.col = c("#4F7A5D", "#045a8d")
)

ggsave(Venn_detected_genes,
       filename = "doc/figures/figure_1_S1/figure_1_S1G.png",
       width = 60,
       height =60,
       units="mm")

################################################################################################################################################
####################################     Panel   H    ###################################
################################################################################################################################################

# Calculate Spearman correlation
cor.test(combined_raw$mean_counts_transcriptomics, combined_raw$mean_counts_proteomics,  method = "spearman") # r = 0.38
cor.test(combined_raw$mean_counts_transcriptomics, combined_raw$mean_counts_proteomics,  method = "pearson") # r = 0.52 (non-log)
cor <- cor.test(log(combined_raw$mean_counts_transcriptomics,10), log(combined_raw$mean_counts_proteomics,10),  method = "pearson") # r = 0.52 (log)

round(cor$estimate, 3)
formatC(cor$p.value, format = "e", digits = 2)

# Create plot

TvsP_plot <- combined_raw %>%
    ggplot() +
    geom_point(aes(x=log(mean_counts_transcriptomics,10), y=log(mean_counts_proteomics,10)), size=0.2) +
    geom_smooth(aes(x=log(mean_counts_transcriptomics,10), y=log(mean_counts_proteomics,10)), method='lm', formula= y~x, se=F, colour = "#28666E", size=0.75) +
    xlab("Transcriptomics (avg counts, log10)") +
    ylab("Proteomics (avg intensity, log10)") +
    annotate("text", x=3, y=3, label= "r = 0.519", colour="black", fontface=2, size=2.5) +
    annotate("text", x=3, y=2.75, label= "p < 0.001", colour="black", fontface=2, size=2.5) +
    theme_classic() +
    theme(
        text = element_text(face="bold", size=6, colour="black"),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        legend.position = "none"
    )

ggsave(TvsP_plot, filename = "doc/figures/figure_1_S1/figure_1_S1H.png", width = 60, height = 60, units="mm")

