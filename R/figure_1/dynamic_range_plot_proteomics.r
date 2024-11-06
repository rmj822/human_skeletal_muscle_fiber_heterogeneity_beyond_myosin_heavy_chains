# Here I will create a dynamic range plot
# On the X axis I will show the proteins ranked by
# median LFQ intensity, and on the Y axis
# The log2 transformed LFQ intensities

proteomics_data <- vroom::vroom(
    here::here("C:/Users/jns822/Desktop/Scripts/Heterofiber/data/data_proteomics_filtered.csv"),
    col_select = !c(1)
) |>
    as.data.frame() |>
    dplyr::filter(
        !Gene.name == "",
        !duplicated(Gene.name)
    ) |>
    tibble::column_to_rownames("Gene.name")

sum_of_intensities <- colSums(proteomics_data,
                              na.rm = TRUE)

rel_abundance <- proteomics_data |>
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
        order = 1:nrow(proteomics_data)
    )

Gene_names <- rownames(mean_rel_abundance)

mean_expression_all <- mean_rel_abundance |>
    tibble::rownames_to_column("gene")

my_colors <- c("gray", "#583B2BFF", "#534C53FF", "#DAEAF6", "#AD8152FF", "#BBA78CFF")

nonmuscle <-mean_rel_abundance |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(gene == "PAX7" | gene == "NCAM1" | gene ==  "MRC1" | gene ==  "C1QA" | gene ==  "CDH5" | gene ==  "MYH11" | gene ==  "PDGFRA" | gene ==  "DCN" | gene == "")

nonmuscle$type <- c("FAP", "SC", "SMC", "EC")



# Load filtered dataset ---------------------------------------------------
library(ggplot2)
library(tidyverse)
library(ggrepel)

proteomics_filtered <- vroom::vroom(here::here("data/proteomics_working_data.csv")) |>
    dplyr::rename("gene" = "...1")

mean_expression_filtered <- mean_expression_all |>
    dplyr::filter(gene %in% proteomics_filtered$gene)

ggplot2::ggplot() +

    # Add all genes
    ggplot2::geom_point(data = mean_expression_all,
                        aes(x=order,
                            y=log10_rel_abundance),
                        colour = "#045a8d",
                        size = 0.25,
                        alpha = 0.5) +

    # Add horizontal lines to indicate % instead of log scale
    geom_hline(yintercept = log(20,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(5,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(1,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(0.1,10), linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = log(0.01,10), linetype="dashed", linewidth=0.2) +

    # Add text for % expression
    annotate("text", x=2900, y=log(28,10), label= "20%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(7,10), label= "5%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(1.4,10), label= "1%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(0.14,10), label= "0.1%", colour="black", fontface=2, size=2) +
    annotate("text", x=2900, y=log(0.014,10), label= "0.01%", colour="black", fontface=2, size=2) +

    # Add custom labels:
    geom_label_repel(data = mean_expression_filtered %>% dplyr::filter(gene %in% c("MYH7",
                                                                                   "MYH2",
                                                                                   "ACTA1",
                                                                                   "TNNT1",
                                                                                   "TNNT3")),
                     mapping = aes(order, log10_rel_abundance, label = gene),
                     size = 1.8, label.padding=0.1, max.overlaps = Inf, min.segment.length=0.1, segment.size=0.2, force = 10) +

    # Add labels non muscle cell markers
    # geom_label_repel(nonmuscle,
    #                  mapping = aes(order, log10_rel_abundance, label = gene, fill=type),
    #                  size = 1.8, max.overlaps = Inf, label.padding=0.1, min.segment.length=0.1, segment.size=0.2, force = 10) +
    # scale_fill_manual(values = c("#E8DFF5", "#FCE1E4", "#9CADCE", "#DAEAF6", "#FCF4DD")) +

    # # Add labels for legend of cell types
    # annotate("rect", xmin = 575, xmax = 2425, ymin = -3.4, ymax = -3.9, fill="#FCE1E4", colour="black") +
    # annotate("rect", xmin = 575, xmax = 2425, ymin = -4, ymax = -4.5, fill="#DAEAF6", colour="black") +
    # annotate("rect", xmin = 575, xmax = 2425, ymin = -4.6, ymax = -5.1, fill="#E8DFF5", colour="black") +
    # annotate("rect", xmin = 575, xmax = 2425, ymin = -5.2, ymax = -5.7, fill="#D6C0A9FF", colour="black") +
    # annotate("rect", xmin = 575, xmax = 2425, ymin = -5.8, ymax = -6.3, fill="#FCF4DD", colour="black") +
    #
    # annotate("text", x=1500, y=-3.65, label= "FAPs", colour="black", fontface=2, size=2.5) +
    # annotate("text", x=1500, y=-4.25, label= "Smooth muscle cells", colour="black", fontface=2, size=2) +
    # annotate("text", x=1500, y=-4.85, label= "Endothelial cells", colour="black", fontface=2, size=2) +
    # annotate("text", x=1500, y=-5.45, label= "Satellite cells", colour="black", fontface=2, size=2) +
    # annotate("text", x=1500, y=-6.05, label= "Macrophages", colour="black", fontface=2, size=2) +

    # Add dots non muscle cell markers

    # # Add PAX7 (Satellite cells)
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PAX7"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "NCAM1"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    # # Add MRC1 (Macrophages)
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MRC1"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    #
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "C1QA"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    #
    # # Add CDH5 (Endothelial cells)
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "CDH5"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    # # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PECAM1"),
    # #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    # # Add ACTA2 (Smooth muscle cells)
    # # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "ACTA2"),
    # #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MYH11"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    #
    # # Add PDGFRA (FAP cells)
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PDGFRA"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +
    #
    # geom_point(data = mean_expression_all %>% dplyr::filter(gene == "DCN"),
    #            aes(x=order, y=log10_rel_abundance), colour="black", size=0.25) +


    # Change design
    ylab("% total intensities, 10log") +
    xlab("Protein rank") +
    theme_classic() +
    ggtitle("Proteomics") +
    theme(
        text = element_text(face="bold", colour="black", size = 6),
        plot.title = element_text(face = "bold", color = "black", size = 8, hjust = 0.5),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        legend.position = "none",
    )

# ggsave(here::here("doc/figures/figure_1/dynamic_range_proteome.png"),
#        units = "mm",
#        height = 60,
#        width = 60)
