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








