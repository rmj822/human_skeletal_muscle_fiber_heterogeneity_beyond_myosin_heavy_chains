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


# Get counts_transcriptomics --------------------------------------------------------------
counts_transcriptomics <- vroom::vroom(here::here("data-raw/filtered_counts_transcriptomics.csv")) |>
    tibble::column_to_rownames("...1")


# Load Proteomics data ----------------------------------------------------
proteomics <- vroom::vroom(here::here("data-raw/raw_proteomics.csv"))
proteomics <- proteomics %>%
    filter(!duplicated(Gene_name)) |>
    filter(!is.na(Gene_name)) |>
    tibble::column_to_rownames("Gene_name")

################################################################################################################################################
#############################################@####       EXTRACT GENES      ####################################################################
################################################################################################################################################


# Extract MYH1 ------------------------------------------------------------
MYH1 <- c("MYH1")
counts_MYH1 <- counts_transcriptomics[MYH1, ] |>
    t()
counts_MYH1 <- as.data.frame(as.matrix(counts_MYH1))


# Extract MYH2 ------------------------------------------------------------
MYH2 <- c("MYH2")
counts_MYH2 <- counts_transcriptomics[MYH2, ] |>
    t()
counts_MYH2 <- as.data.frame(as.matrix(counts_MYH2))

# Extract MYH7 ------------------------------------------------------------
MYH7 <- c("MYH7")
counts_MYH7 <- counts_transcriptomics[MYH7, ] |>
    t()
counts_MYH7 <- as.data.frame(as.matrix(counts_MYH7))


# Merge dataframes --------------------------------------------------------
counts_ft <- merge(counts_MYH1, counts_MYH2, by="row.names")
counts_ft <- merge(counts_ft, counts_MYH7 |>
                       tibble::rownames_to_column("Row.names"), by="Row.names")

# Delete Row.names column -------------------------------------------------
rownames(counts_ft) <- counts_ft$Row.names
counts_ft <- counts_ft %>% dplyr::select(MYH1, MYH2, MYH7)
# Remove any NA values
counts_ft <- counts_ft %>%
  tidyr::replace_na(list(
    "MYH7" = 0,
    "MYH2" = 0,
    "MYH1" = 0
  ))


################################################################################################################################################
#############################      CALCULATE PERCENTAGE UMI counts_transcriptomics FOR EACH GENE COMBINATION    ################################################
################################################################################################################################################

# MYH -------------------------------------------------------

# Sum of UMI counts_transcriptomics
counts_ft$sum_MYH <- counts_ft$MYH7 + counts_ft$MYH2 + counts_ft$MYH1

# MYH7
counts_ft$MYH7_fraction <- counts_ft$MYH7 / counts_ft$sum_MYH * 100
counts_ft <- counts_ft %>%
  dplyr::arrange(desc(MYH7_fraction)) %>%
  dplyr::mutate("order_MYH7" = seq_len(nrow(counts_ft)))

# MYH2
counts_ft$MYH2_fraction <- counts_ft$MYH2 / counts_ft$sum_MYH * 100
counts_ft <- counts_ft %>%
  dplyr::arrange(desc(MYH2_fraction)) %>%
  dplyr::mutate("order_MYH2" = seq_len(nrow(counts_ft)))

# MYH1
counts_ft$MYH1_fraction <- counts_ft$MYH1 / counts_ft$sum_MYH * 100
counts_ft <- counts_ft %>%
  dplyr::arrange(desc(MYH1_fraction)) %>%
  dplyr::mutate("order_MYH1" = seq_len(nrow(counts_ft)))


# MYH - only MYH7 and MYH2 -------------------------------------------------------

# Sum of UMI counts_transcriptomics
counts_ft$sum_MYH_slowfast <- counts_ft$MYH7 + counts_ft$MYH2

# MYH7
counts_ft$MYH7_fraction_slowfast <- counts_ft$MYH7 / counts_ft$sum_MYH_slowfast * 100
counts_ft <- counts_ft %>%
    dplyr::arrange(desc(MYH7_fraction_slowfast)) %>%
    dplyr::mutate("order_MYH7_slowfast" = seq_len(nrow(counts_ft)))

# MYH2
counts_ft$MYH2_fraction_slowfast <- counts_ft$MYH2 / counts_ft$sum_MYH_slowfast * 100
counts_ft <- counts_ft %>%
    dplyr::arrange(desc(MYH2_fraction_slowfast)) %>%
    dplyr::mutate("order_MYH2_slowfast" = seq_len(nrow(counts_ft)))


################################################################################################################################################
#############################      CALCULATE PERCENTAGE MYH isofroms proteomics    ################################################
################################################################################################################################################

# Selecting only MYH
proteomics_ft <- proteomics |>
  t() |>
  as.data.frame() |>
  dplyr::select("MYH7", "MYH2", "MYH1") |>
  dplyr::mutate("Fiber_ID" = colnames(proteomics)) |>
  dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
  dplyr::arrange(desc(MYH7)) |>
  dplyr::mutate("order_7" = seq_len(ncol(proteomics))) |>
  dplyr::arrange(desc(MYH2)) |>
  dplyr::mutate("order_2" = seq_len(ncol(proteomics)))

# QC figure
proteomics_ft |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = order_7,
      y = MYH7
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal()


# Sum of MYH intensities
proteomics_ft$sum_MYH <- proteomics_ft$MYH7 + proteomics_ft$MYH2 + proteomics_ft$MYH1

# MYH7
proteomics_ft$MYH7_fraction <- proteomics_ft$MYH7 / proteomics_ft$sum_MYH * 100
proteomics_ft <- proteomics_ft %>%
  dplyr::arrange(desc(MYH7_fraction)) %>%
  dplyr::mutate("order_MYH7" = seq_len(nrow(proteomics_ft)))

# MYH2
proteomics_ft$MYH2_fraction <- proteomics_ft$MYH2 / proteomics_ft$sum_MYH * 100
proteomics_ft <- proteomics_ft %>%
  dplyr::arrange(desc(MYH2_fraction)) %>%
  dplyr::mutate("order_MYH2" = seq_len(nrow(proteomics_ft)))

# MYH1
proteomics_ft$MYH1_fraction <- proteomics_ft$MYH1 / proteomics_ft$sum_MYH * 100
proteomics_ft <- proteomics_ft %>%
  dplyr::arrange(desc(MYH1_fraction)) %>%
  dplyr::mutate("order_MYH1" = seq_len(nrow(proteomics_ft)))


################################################################################################################################################
################################################      CHECK CURVES AS QC: Tx    ###################################################################
################################################################################################################################################


# MYH7 --------------------------------------------------------------------
counts_ft %>%
  dplyr::select(MYH7_fraction, order_MYH7) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH7,
                 MYH7_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH7")

# MYH2 --------------------------------------------------------------------
counts_ft %>%
  dplyr::select(MYH2_fraction, order_MYH2) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH2,
                 MYH2_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH2")


# MYH1 --------------------------------------------------------------------
counts_ft %>%
  dplyr::select(MYH1_fraction, order_MYH1) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH1,
                 MYH1_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH1")


################################################################################################################################################
################################################      FINDING BOTTOM KNEE     ##################################################################
################################################################################################################################################

# Function detects top knee, so we take inverse for bottom knee --> used as threshold

# MYH7 --------------------------------------------------------------------
MYH7_curve <- counts_ft %>%
  dplyr::select(MYH7_fraction) %>%
  dplyr::mutate(across(everything(), ~ 100 - .x)) %>%
  t()

MYH7_curve <- DropletUtils::barcodeRanks(MYH7_curve, lower = 10)

bottom_knee_MYH7 <- 100 - MYH7_curve@metadata$knee

MYH7_curve_plot <- counts_ft %>%
  dplyr::select(MYH7_fraction, order_MYH7) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH7,
                 MYH7_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH7") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH7, color = "red") +
  annotate("text", x=50, y=(bottom_knee_MYH7+3), label= "6%", colour="black", fontface=2)


# MYH2 --------------------------------------------------------------------
MYH2_curve <- counts_ft %>%
  dplyr::select(MYH2_fraction) %>%
  dplyr::mutate(across(everything(), ~ 100 - .x)) %>%
  t()

MYH2_curve <- DropletUtils::barcodeRanks(MYH2_curve, lower = 10)

bottom_knee_MYH2 <- 100 - MYH2_curve@metadata$knee

MYH2_curve_plot <- counts_ft %>%
  dplyr::select(MYH2_fraction, order_MYH2) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH2,
                 MYH2_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH2") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH2, color = "red") +
  annotate("text", x=50, y=(bottom_knee_MYH2+3), label= "13%", colour="black", fontface=2)


# MYH1 --------------------------------------------------------------------
MYH1_curve <- counts_ft %>%
  dplyr::select(MYH1_fraction) %>%
  dplyr::mutate(across(everything(), ~ 100 - .x)) %>%
  t()

MYH1_curve <- DropletUtils::barcodeRanks(MYH1_curve, lower = 10)

bottom_knee_MYH1 <- 100 - MYH1_curve@metadata$knee

MYH1_curve_plot <- counts_ft %>%
  dplyr::select(MYH1_fraction, order_MYH1) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH1,
                 MYH1_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH1") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH1, color = "red") +
  annotate("text", x=50, y=(bottom_knee_MYH1+3), label= "10%", colour="black", fontface=2)


################################################################################################################################################
##################################################     BOTTOM KNEE PROTEOMICS     ##################################################################
################################################################################################################################################

# Replace NA by 0
proteomics_ft <- proteomics_ft %>% tidyr::replace_na(list(
  "MYH7_fraction" = 0,
  "MYH2_fraction" = 0,
  "MYH1_fraction" = 0
))

# MYH7 --------------------------------------------------------------------
MYH7_curve_proteomics <- proteomics_ft %>%
  dplyr::select(MYH7_fraction) %>%
  dplyr::mutate(across(everything(), ~ 100 - .x)) %>%
  t()

MYH7_curve_proteomics <- DropletUtils::barcodeRanks(MYH7_curve_proteomics, lower = 10)

bottom_knee_MYH7_proteomics <- 100 - MYH7_curve_proteomics@metadata$knee

proteomics_ft %>%
  dplyr::select(MYH7_fraction, order_MYH7) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH7,
                 MYH7_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH7") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH7_proteomics, color = "red") +
  annotate("text", x=50, y=(bottom_knee_MYH7+3), label= "31%", colour="black", fontface=2)

# MYH2 --------------------------------------------------------------------
MYH2_curve_proteomics <- proteomics_ft %>%
  dplyr::select(MYH2_fraction) %>%
  dplyr::mutate(across(everything(), ~ 100 - .x)) %>%
  t()

MYH2_curve_proteomics <- DropletUtils::barcodeRanks(MYH2_curve_proteomics, lower = 10)

bottom_knee_MYH2_proteomics <- 100 - MYH2_curve_proteomics@metadata$knee

proteomics_ft %>%
  dplyr::select(MYH2_fraction, order_MYH2) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH2,
                 MYH2_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH2") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH2_proteomics, color = "red") +
  annotate("text", x=50, y=(bottom_knee_MYH2+3), label= "9%", colour="black", fontface=2)

# MYH1 --------------------------------------------------------------------
MYH1_curve_proteomics <- proteomics_ft %>%
  dplyr::select(MYH1_fraction) %>%
  dplyr::mutate(across(everything(), ~ 100 - .x)) %>%
  t()

MYH1_curve_proteomics <- DropletUtils::barcodeRanks(MYH1_curve_proteomics, lower = 10)

bottom_knee_MYH1_proteomics <- 100 - MYH1_curve_proteomics@metadata$knee

proteomics_ft %>%
  dplyr::select(MYH1_fraction, order_MYH1) %>%
  ggplot2::ggplot(
    ggplot2::aes(order_MYH1,
                 MYH1_fraction
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal() +
  ggtitle("MYH1") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH1_proteomics, color = "red") +
  annotate("text", x=50, y=(bottom_knee_MYH1+3), label= "8%", colour="black", fontface=2)

################################################################################################################################################
##################################################    COMBINED PLOTS Tx and Px     #############################################################
################################################################################################################################################


# MYH7 --------------------------------------------------------------------
MYH7_TvsP <- ggplot() +
  geom_point(
    data = counts_ft %>% dplyr::arrange(desc(MYH7_fraction)),
    aes(x = order_MYH7, y = MYH7_fraction),
    colour = "#4F7A5D",
    size = 0.4) +
  geom_point(
    data = proteomics_ft %>% dplyr::arrange(desc(MYH7_fraction)),
    aes(x = order_MYH7, y = MYH7_fraction),
    colour = "#045a8d",
    size = 0.4) +

  geom_hline(yintercept = bottom_knee_MYH7, color = "#4F7A5D", alpha=0.8) +
  geom_hline(yintercept = bottom_knee_MYH7_proteomics, color = "#045a8d", alpha=0.8) +

  annotate("text", x=-10, y=(bottom_knee_MYH7-10), label= "Transcriptomics: 6%", colour="#4F7A5D", fontface=2, hjust = 0, size=1.75) +
  annotate("text", x=450, y=(bottom_knee_MYH7_proteomics+4), label= "Proteomics: 31%", colour="#045a8d", fontface=2, hjust = 0, size=1.75) +

  ylab("MYH7 fraction (%)") +
  xlab("Fiber rank") +
  theme_classic() +
  theme(
    text = element_text(face="bold", colour="black", size=6),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
  )

ggsave(MYH7_TvsP, filename = here::here("doc/figures/figure_1/threshold_MYH7.png"), width = 40, height = 40, units="mm")


# MYH2 --------------------------------------------------------------------
MYH2_TvsP <- ggplot() +
  geom_point(
    data = counts_ft %>% dplyr::arrange(desc(MYH2_fraction)),
    aes(x = order_MYH2, y = MYH2_fraction),
    colour = "#4F7A5D",
    size = 0.4) +
  geom_point(
    data = proteomics_ft %>% dplyr::arrange(desc(MYH2_fraction)),
    aes(x = order_MYH2, y = MYH2_fraction),
    colour = "#045a8d",
    size = 0.4) +

  geom_hline(yintercept = bottom_knee_MYH2, color = "#4F7A5D", alpha=0.8) +
  geom_hline(yintercept = bottom_knee_MYH2_proteomics, color = "#045a8d", alpha=0.8) +

  annotate("text", x=0, y=(bottom_knee_MYH2-14), label= "Transcriptomics: 13%", colour="#4F7A5D", fontface=2, hjust = 0, size=1.75) +
  annotate("text", x=0, y=(bottom_knee_MYH2_proteomics-4), label= "Proteomics: 9%", colour="#045a8d", fontface=2, hjust = 0, size=1.75) +

  ylab("MYH2 fraction (%)") +
  xlab("Fiber rank") +
  theme_classic() +
  theme(
    text = element_text(face="bold", colour="black", size=6),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
  )

ggsave(MYH2_TvsP, filename = here::here("doc/figures/figure_1/threshold_MYH2.png"), width = 40, height = 40, units="mm")

# MYH1 --------------------------------------------------------------------
MYH1_TvsP <- ggplot() +
  geom_point(
    data = counts_ft %>% dplyr::arrange(desc(MYH1_fraction)),
    aes(x = order_MYH1, y = MYH1_fraction),
    colour = "#4F7A5D",
    size = 0.4) +
  geom_point(
    data = proteomics_ft %>% dplyr::arrange(desc(MYH1_fraction)),
    aes(x = order_MYH1, y = MYH1_fraction),
    colour = "#045a8d",
    size = 0.4) +

  geom_hline(yintercept = bottom_knee_MYH1, color = "#4F7A5D", alpha=0.8) +
  geom_hline(yintercept = bottom_knee_MYH1_proteomics, color = "#045a8d", alpha=0.8) +

  annotate("text", x=300, y=(bottom_knee_MYH1+10), label= "Transcriptomics: 10%", colour="#4F7A5D", fontface=2, hjust = 0, size=1.75) +
  annotate("text", x=300, y=(bottom_knee_MYH1_proteomics+6), label= "Proteomics: 8%", colour="#045a8d", fontface=2, hjust = 0, size=1.75) +

  ylab("MYH1 fraction (%)") +
  xlab("Fiber rank") +
  theme_classic() +
  theme(
    text = element_text(face="bold", colour="black", size=6),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
  )

ggsave(MYH1_TvsP, filename = here::here("doc/figures/figure_1/threshold_MYH1.png"), width = 40, height = 40, units="mm")


################################################################################################################################################
##################################################     ASSIGN FIBER TYPES     ##################################################################
################################################################################################################################################

# MYH --------------------------------------------------------------------
counts_ft <- counts_ft %>%
  dplyr::mutate(
    fiber_type_MYH = dplyr::case_when(
      MYH7_fraction >= bottom_knee_MYH7 &
        MYH2_fraction <= bottom_knee_MYH2 &
          MYH1_fraction <= bottom_knee_MYH1 ~ "Type 1",
      MYH7_fraction <= bottom_knee_MYH7 &
        MYH2_fraction >= bottom_knee_MYH2 &
          MYH1_fraction <= bottom_knee_MYH1 ~ "Type 2A",
      MYH7_fraction <= bottom_knee_MYH7 &
        MYH2_fraction <= bottom_knee_MYH2 &
          MYH1_fraction >= bottom_knee_MYH1 ~ "Type 2X",
      MYH7_fraction <= bottom_knee_MYH7 &
        MYH2_fraction >= bottom_knee_MYH2 &
          MYH1_fraction >= bottom_knee_MYH1 ~ "Hybrid 2A/2X",
      MYH7_fraction >= bottom_knee_MYH7 &
        MYH2_fraction >= bottom_knee_MYH2 &
          MYH1_fraction <= bottom_knee_MYH1 ~ "Hybrid 1/2A",
      MYH7_fraction >= bottom_knee_MYH7 &
        MYH2_fraction <= bottom_knee_MYH2 &
          MYH1_fraction >= bottom_knee_MYH1 ~ "Hybrid 1/2X",
      MYH7_fraction >= bottom_knee_MYH7 &
        MYH2_fraction >= bottom_knee_MYH2 &
          MYH1_fraction >= bottom_knee_MYH1~ "Hybrid 1/2A/2X",
      TRUE ~ "NA"
    )
  )

# Check if any NA values
summary(is.na(counts_ft$fiber_type_MYH))

# MYH Slow vs Fast--------------------------------------------------------------------
counts_ft <- counts_ft %>%
    dplyr::mutate(
        fiber_type_MYH_slowfast = dplyr::case_when(
            fiber_type_MYH == "Type 1" ~ "Slow",
            fiber_type_MYH == "Type 2A" ~ "Fast",
            fiber_type_MYH == "Type 2X" ~ "Fast",
            fiber_type_MYH == "Hybrid 2A/2X" ~ "Fast",
            fiber_type_MYH == "Hybrid 1/2A" ~ "Hybrid",
            fiber_type_MYH == "Hybrid 1/2X" ~ "Hybrid",
            fiber_type_MYH == "Hybrid 1/2A/2X" ~ "Hybrid",
            TRUE ~ "Fast"
        )
    )



# Check if any NA values
summary(is.na(counts_ft$fiber_type_MYH_slowfast))


# MYH-based: simpler version ----------------------------------------------

counts_ft <- counts_ft %>%
    dplyr::mutate(
        fiber_type_MYH_hybrids = dplyr::case_when(
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

counts_ft$fiber_type_MYH_hybrids <- factor(counts_ft$fiber_type_MYH_hybrids, levels = c("Type 1", "Hybrid 1/2A", "Type 2A", "Hybrid 2A/2X", "Type 2X"))


################################################################################################################################################
###################################################     SAVE DATA FRAME     #########>>#########################################################
################################################################################################################################################
counts_ft$Sample <- rownames(counts_ft)

# Add subject and fiber column --------------------------------------------
counts_ft$subject <- gsub("_.*", "", counts_ft$Sample)
counts_ft$fiber <- gsub(".*_", "", counts_ft$Sample)


write_csv(counts_ft, file="~/single_fiber_heterogeneity/data/counts_ft.csv")
