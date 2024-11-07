
# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)

################################################################################################################################################
########################################################       FIGURE 5 S1A-B    ###################################################################
################################################################################################################################################

# No custom code


################################################################################################################################################
########################################################       FIGURE 5 S1C    ###################################################################
################################################################################################################################################

source(here::here("R/figure_5/figure_5.R"))

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

ggplot2::ggsave(
    here::here("doc/figures/figure_5_S1/figure_5_S1C.pdf"),
    height = 60,
    width = 60,
    units = "mm"
)


################################################################################################################################################
########################################################       FIGURE 5 S1D    ###################################################################
################################################################################################################################################

# LINC01405_ORF420:14512:14291

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
        dplyr::filter(transcript_id == "LINC01405_ORF420:14512:14291") |>
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
        dplyr::filter(transcript_id == "LINC01405_ORF420:14512:14291") |>
        dplyr::pull(LFQ_intensities, name = sample_id),
    paired = TRUE
)

# data_boxplots <- data_boxplots |>
#     dplyr::mutate(
#         sample_id = gsub(
#             pattern = "_.*",
#             replacement = "",
#             sample_id
#         )
#     )

LINC01405_ORF420 <- data_boxplots |>
    dplyr::mutate(subject =
                      gsub(
                          pattern = "_.*",
                          replacement = "",
                          sample_id
                      )) |>
    dplyr::filter(transcript_id == "LINC01405_ORF420:14512:14291") |>
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
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 4.9, label = paste(round(
        DE_results |>
            dplyr::filter(grepl("ORF420:14512:14291",Genes)) |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("LINC01405_ORF420:14512:14291") +
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
        plot.margin = grid::unit(c(1,1,0,0), "mm"))


# LINC01405_ORF403:18859:18767

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
        dplyr::filter(transcript_id == "LINC01405_ORF403:18859:18767") |>
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
        dplyr::filter(transcript_id == "LINC01405_ORF403:18859:18767") |>
        dplyr::pull(LFQ_intensities, name = sample_id),
    paired = TRUE
)

# data_boxplots <- data_boxplots |>
#     dplyr::mutate(
#         sample_id = gsub(
#             pattern = "_.*",
#             replacement = "",
#             sample_id
#         )
#     )

LINC01405_ORF403 <- data_boxplots |>
    dplyr::mutate(subject =
                      gsub(
                          pattern = "_.*",
                          replacement = "",
                          sample_id
                      )) |>
    dplyr::filter(transcript_id == "LINC01405_ORF403:18859:18767") |>
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
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 4.05, label = paste(round(
        DE_results |>
            dplyr::filter(grepl("ORF403:18859:18767",Genes)) |>
            dplyr::pull(adj.P.Val),
        4
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#440154FF",
                                        "#5DC863FF")) +
    ggplot2::ggtitle("LINC01405_ORF403:18859:18767") +
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

patchwork::wrap_plots(LINC01405_ORF420,
                      LINC01405_ORF403,
                      ncol = 2,
                      nrow = 1)

ggplot2::ggsave(here::here("doc/figures/figure_5_S1/figure_5_S1D.pdf"),
                units = "mm",
                width = 90,
                height = 60)
