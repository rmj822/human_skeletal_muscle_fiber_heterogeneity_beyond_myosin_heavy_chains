################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Matrix)
library(Seurat)
library(viridis)
library(rstatix)
library(ggnewscale)
library(ComplexUpset)
library(viridis)
library(ggpubr)
library(ggdist)
library(janitor)
library(DropletUtils)
library(ggrepel)
library(Hmisc)
library(grid)

# Read fiber type file --------------------------------------------------------------
counts_ft <- read_csv(file=here::here("data/counts_ft.csv"))

################################################################################################################################################
############################################       MAIN FIGURE: MYH RATIO'S     ####@###########################################################
################################################################################################################################################

# Check number of fibers for each type
counts_ft %>% dplyr::filter(fiber_type_MYH == "Type 1") %>% nrow() # 324
counts_ft %>% dplyr::filter(fiber_type_MYH == "Hybrid 1/2X") %>% nrow() # 1

counts_ft %>% dplyr::filter(fiber_type_MYH == "Hybrid 1/2A") %>% nrow() # 41
counts_ft %>% dplyr::filter(fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction >= MYH1_fraction) %>% nrow() # 4

counts_ft %>% dplyr::filter(fiber_type_MYH == "Type 2A") %>% nrow() # 336

counts_ft %>% dplyr::filter(fiber_type_MYH == "Hybrid 2A/2X") %>% nrow() # 140
counts_ft %>% dplyr::filter(fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction < MYH1_fraction) %>% nrow() # 10

counts_ft %>% dplyr::filter(fiber_type_MYH == "Type 2X") %>% nrow() # 69

# Make list for order of figure
fiber_types <- list(
  type_1 <- counts_ft |>
    dplyr::filter(fiber_type_MYH == "Type 1" | fiber_type_MYH == "Hybrid 1/2X") |>
    dplyr::arrange(desc(MYH7_fraction)),
  Hybrid_1_2A <- counts_ft |>
    dplyr::filter(fiber_type_MYH == "Hybrid 1/2A" | fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction >= MYH1_fraction) |>
    dplyr::arrange(desc(MYH7_fraction)),
  type_2A <- counts_ft |>
    dplyr::filter(fiber_type_MYH == "Type 2A") |>
    dplyr::arrange(MYH2_fraction),
  Hybrid_2A_2X <- counts_ft |>
    dplyr::filter(fiber_type_MYH == "Hybrid 2A/2X" | fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction < MYH1_fraction) |>
    dplyr::arrange(MYH1_fraction),
  type_2X <- counts_ft |>
    dplyr::filter(fiber_type_MYH == "Type 2X") |>
    dplyr::arrange(MYH1_fraction)
)

# Give names to list
names(fiber_types) <- c(
  "type_1",
  "hybrid_1_2A",
  "type_2A",
  "hybrid_2A_2X",
  "type_2X"
)

# Make ordered matrix
order_matrix <- dplyr::bind_rows(
  fiber_types$type_1,
  fiber_types$hybrid_1_2A,
  fiber_types$type_2A,
  fiber_types$hybrid_2A_2X,
  fiber_types$type_2X
) |>
  tibble::add_column("sample_MYH" = seq_len(nrow(counts_ft)))

# Make df for rectangle colours
rectangle_colors <- data.frame(
  start = c(0, 324.5, 369.5, 705.5, 855.5),
  end = c(324.5, 369.5, 705.5, 855.5, 924),
  fiber_type = c("Type I", "Hybrid I/IIA", "Type IIA", "Hybrid IIA/IIX", "Type IIx")
)


# Create figure 128 x 60

text_MYH7v2 <- textGrob("MYH7 threshold", gp=gpar(fontsize=5))
text_MYH1v2 <- textGrob("MYH1 threshold", gp=gpar(fontsize=5))
text_MYH2v2 <- textGrob("MYH2 threshold", gp=gpar(fontsize=5))

plot_MYH_v2 <- order_matrix |>
  ggplot2::ggplot() +
  coord_cartesian(xlim = c(0,925), ylim = c(0,120), clip = "off") +
  ggplot2::geom_rect(
    data = rectangle_colors,
    ggplot2::aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf,
      fill = fiber_type
    ),
    alpha = 0.15
  ) +
  ggplot2::scale_fill_manual(values = c(
    "#3B528BFF",
               "#fdc325",
               "#440154FF",
               "#5DC863FF",
               "#D2631C"
  )) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = sample_MYH,
      y = MYH7_fraction
    ),
    colour = "#440154FF",
    size = 0.25
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = sample_MYH,
      y = MYH2_fraction
    ),
    colour = "#5DC863FF",
    size = 0.25
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = sample_MYH,
      y = MYH1_fraction
    ),
    colour = "#D2631C",
    size = 0.25
  ) +
  ggplot2::ylab("% MYH isoform expressed") +
  ggplot2::xlab("Sample (ranked by MYH expression)") +
  ggplot2::ggtitle("MYH-based fiber typing (transcriptomics)") +
  ggplot2::geom_vline(xintercept = 324.5, colour = "grey20", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 369.5, colour = "grey20", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 705.5, colour = "grey20", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 855.5, colour = "grey20", linewidth = 0.25) +
  ggplot2::geom_hline(
    yintercept = 6,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.25
  ) +
  ggplot2::geom_hline(
    yintercept = 13,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.25
  ) +
  ggplot2::geom_hline(
    yintercept = 10,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.25
  ) +
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
  ggplot2::annotate(
    "text",
    x = 162,
    y = 107,
    label = "Type 1 \n 324 (35.1%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 537,
    y = 107,
    label = "Type 2A \n 336 (36.3%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 780,
    y = 107,
    label = "Hybrid 2A/2X \n 150 (16.2%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 348.5,
    y = 107,
    label = "Hybrid \n1/2A \n45 \n(4.9%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 889.5,
    y = 107,
    label = "Type 2X \n 69 (7.5%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 35,
    y = 88,
    label = "MYH7",
    colour = "#440154FF",
    fontface = 2,
    size=2
  ) +
  ggplot2::annotate(
    "text",
    x = 450,
    y = 88,
    label = "MYH2",
    colour = "#5DC863FF",
    fontface = 2,
    size=2
  ) +
  ggplot2::annotate(
    "text",
    x = 900,
    y = 88,
    label = "MYH1",
    colour = "#D2631C",
    fontface = 2,
    size=2
  ) +
  annotation_custom(text_MYH7v2,xmin=60,xmax=60,ymin=-15,ymax=-16) +
  annotation_custom(text_MYH1v2,xmin=180,xmax=180,ymin=-15,ymax=-16) +
  annotation_custom(text_MYH2v2,xmin=60,xmax=60,ymin=15,ymax=20) +
  annotate("segment", x = 60, xend = 60, y = -13, yend = 5, size=0.3) +
  annotate("segment", x = 180, xend = 180, y = -13, yend = 9, size=0.3) +
  ggplot2::theme(text = ggplot2::element_text(size = 7)) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 7)) +
  scale_x_continuous(limits = c(0,925), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", "100"))

ggsave(plot_MYH_v2, filename = here::here("doc/figures/figure_1/MYH_curves_transcriptomics.png"), width = 128, height = 60, units="mm")


################################################################################################################################################
#######################################     NUMBER OF FIBERS PER FIBER TYPE     ################################################################
################################################################################################################################################

# MYH isoforms per subject ------------------------------------------------------------

data_plot_fibertype_participant <- counts_ft %>%
  dplyr::mutate(
    fiber_type_MYH = dplyr::case_when(
      fiber_type_MYH == "Type 1" ~ "Type 1",
      fiber_type_MYH == "Type 2A" ~ "Type 2A",
      fiber_type_MYH == "Type 2X" ~ "Type 2X",
      fiber_type_MYH == "Hybrid 1/2A" ~ "Hybrid 1/2A",
      fiber_type_MYH == "Hybrid 1/2X" ~ "Type 1",
      fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction >= MYH1_fraction ~ "Hybrid 1/2A",
      fiber_type_MYH == "Hybrid 2A/2X" ~ "Hybrid 2A/2X",
      fiber_type_MYH == "Hybrid 1/2A/2X" & MYH7_fraction < MYH1_fraction ~ "Hybrid 2A/2X",
      TRUE ~ "NA"
      ))

data_plot_fibertype_participant$fiber_type_MYH <- factor(data_plot_fibertype_participant$fiber_type_MYH, levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))

plot_fibertype_MYH_subject <- data_plot_fibertype_participant %>%
  ggplot(aes(
    x = subject,
    fill = fiber_type_MYH
  )) +
  geom_bar(na.rm = TRUE, position="fill", alpha = 0.85, colour=NA, size=0) +
  scale_fill_manual(values = c("#440154FF", "#3B528BFF", "#5DC863FF", "#fdc325", "#D2631C")) +
  labs(
    x = "Subject",
    y = "Percentage"
  ) +
  theme_classic() +
  ggtitle("MYH-based fiber typing by participant (transcriptomics)") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = c("0", "25", "50", "75", "100")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")) +
  theme(
    text = element_text(face="bold", colour="black", size=8),
    plot.title = element_text(size = 8,face = "bold", hjust = 0.5),
    legend.position = "none",
  )

ggsave(plot_fibertype_MYH_subject, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/fiber_type_subject_transcriptomics.png", width = 128, height = 60, units="mm")



################################################################################################################################################
#################################       FIBER TYPE DISTRIBUTION PER PARTICIPANT (MYH SLOW VS FAST)       #######################################
################################################################################################################################################

counts_ft <- counts_ft %>%
  mutate(biopsy = ifelse(fiber >= 1 & fiber <= 25 , 1,
                                       ifelse(fiber > 75 & fiber <= 100 , 2,
                                                            ifelse(fiber > 150 & fiber <= 175 , 3,   NA))))


# Distribution per subject ------------------------------------------------

# Order subjects from high to low slow expression
counts_ft$subject <- factor(counts_ft$subject, levels = c(4, 8, 12, 13, 10, 11, 9, 2, 7, 1, 3, 5, 6, 14))

ft_distribution_subject <- counts_ft %>%
  ggplot(aes(x=subject, fill=fiber_type_MYH_slowfast)) +
  geom_bar(position="stack") +
  scale_fill_manual(values=c("#BC4749", "#8CB3E8", "#134057")) +
  labs(
    x = "Subject",
    y = "N fibers",
    title = "Fiber type distribution"
  )  +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text = element_text(colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(ft_distribution_subject, filename = "7 Targeted fiber type markers/Fibers at rest/Inflection/Figures and metrics/distribution per subject.pdf", width = 8, height = 5, scale=1.5)


# Distribution per biopsy per subject ------------------------------------------------

# Order subjects from 1 to 14
counts_ft$subject <- factor(counts_ft$subject, levels = c(1:14))


ft_distribution_biopsy <- counts_ft %>%
  ggplot(aes(x=biopsy, fill=fiber_type_MYH_slowfast)) +
  geom_bar(position="stack") +
  scale_fill_manual(values=c("#BC4749", "#8CB3E8", "#134057")) +
  labs(
    x = "Biopsy",
    y = "N fibers",
    title = "Fiber type distribution per biopsy"
  )  +
  theme_classic() +
  scale_x_continuous(breaks=c(1,2,3), labels=c("Day 1", "Day 2", "Day 3"),
                     limits=c(0.2, 3.8)) +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text = element_text(colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)


ggsave(ft_distribution_biopsy, filename = "7 Targeted fiber type markers/Fibers at rest/Inflection/Figures and metrics/distribution per biopsy.pdf", width = 8, height = 5, scale=1.5)

################################################################################################################################################
########################################      DENSITY PLOTS MYH DISTRIBUTION       ####@########################################################
################################################################################################################################################


# MYH7 --------------------------------------------------------------------
density_MYH7 <- counts_ft %>%
  ggplot() +
  geom_density(aes(x=MYH7_fraction), color="#134057", fill="#134057", alpha=0.5) +
  ylab("Fraction") +
  ggtitle("Density plot: MYH7") +
  geom_vline(xintercept = 20) +
  geom_vline(xintercept = 80) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# MYH2 --------------------------------------------------------------------
density_MYH2 <- counts_ft %>%
  ggplot() +
  geom_density(aes(x=MYH2_fraction), color="#BC4749", fill="#BC4749", alpha=0.5) +
  ylab("Fraction") +
  ggtitle("Density plot: MYH2") +
  geom_vline(xintercept = 20) +
  geom_vline(xintercept = 80) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# MYH1 --------------------------------------------------------------------
density_MYH1 <- counts_ft %>%
  ggplot() +
  geom_density(aes(x=MYH1_fraction), color="#8CB3E8", fill="#8CB3E8", alpha=0.5) +
  ylab("Fraction") +
  ggtitle("Density plot: MYH1") +
  geom_vline(xintercept = 20) +
  geom_vline(xintercept = 80) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Combined MYH --------------------------------------------------------------------
density_combined <- counts_ft %>%
  ggplot() +
  geom_density(aes(x=MYH7_fraction), color="#134057", fill="#134057", alpha=0.5) +
  geom_density(aes(x=MYH2_fraction), color="#BC4749", fill="#BC4749", alpha=0.5) +
  geom_density(aes(x=MYH1_fraction), color="#8CB3E8", fill="#8CB3E8", alpha=0.5) +
  ylab("Fraction") +
  ggtitle("Density plot: combined") +
  geom_vline(xintercept = 20) +
  geom_vline(xintercept = 80) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


# Arranged plot -----------------------------------------------------------
density_plot <- ggpubr::ggarrange(density_MYH7, density_MYH2, density_MYH1, density_combined)
ggsave(density_plot, filename = "7 Targeted fiber type markers/Fibers at rest/Inflection/Figures and metrics/density MYH.pdf", width = 6, height = 5, scale=1.5)
