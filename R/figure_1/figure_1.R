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
##################################################    PANEL B, C and D     #############################################################
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


# Assign fiber types

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


# Fiber typing curves -----------------------------------------------------

# Proteomics

# Load your data.
# Im using a format where samples are columns and genes are rows.
# Im sending my MYH-filtered data as an example:

proteomics_MYH_filtered <- vroom::vroom(
    here::here("data-raw/raw_proteomics.csv")
) |>
    as.data.frame() |>
    dplyr::filter(!duplicated(Gene_name)) |>
    dplyr::filter(!is.na(Gene_name)) |>
    tibble::column_to_rownames("Gene_name")

# Sorting dataset with the genes that are used for fiber typing,
# plotting intensities and rank:
data_fiber_typing <- proteomics_MYH_filtered |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2", "MYH1") |>
    dplyr::mutate("Fiber_ID" = colnames(proteomics_MYH_filtered)) |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
    dplyr::arrange(desc(MYH7)) |>
    dplyr::mutate("order_7" = seq_len(ncol(proteomics_MYH_filtered))) |>
    dplyr::arrange(desc(MYH2)) |>
    dplyr::mutate("order_2" = seq_len(ncol(proteomics_MYH_filtered)))

data_fiber_typing |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = order_7,
            y = MYH7
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal()


# Ben and I double checked Marta Murgia's paper on SMF proteomics
# and saw that her curve figure is made using percentages
# of MYH intensities, instead of raw intensities.
# Here, we will do the same, calculate percentage
# of MYH composition and rank fibers before plotting them:

sum_of_intensities <- proteomics_MYH_filtered |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2", "MYH1") |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
    tidyr::replace_na(list(
        "MYH7" = 0,
        "MYH2" = 0,
        "MYH1" = 0
    )) |>
    rowSums()


perc_MYHs <- proteomics_MYH_filtered |>
    t() |>
    as.data.frame() |>
    dplyr::select("MYH7", "MYH2", "MYH1") |>
    dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
    tidyr::replace_na(list(
        "MYH7" = 0,
        "MYH2" = 0,
        "MYH1" = 0
    )) |>
    tibble::add_column(sum_of_intensities) |>
    dplyr::mutate(across(
        .cols = c("MYH7", "MYH2", "MYH1"),
        ~ .x / sum_of_intensities * 100
    )) |>
    as.data.frame() |>
    dplyr::arrange(desc(MYH7)) |>
    dplyr::mutate("order_7" = seq_len(ncol(proteomics_MYH_filtered))) |>
    dplyr::arrange(desc(MYH2)) |>
    dplyr::mutate("order_2" = seq_len(ncol(proteomics_MYH_filtered))) |>
    dplyr::arrange(desc(MYH1)) |>
    dplyr::mutate("order_1" = seq_len(ncol(proteomics_MYH_filtered))) |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::arrange(desc(fiber_ID)) |>
    dplyr::select(c(
        "MYH7",
        "MYH2",
        "MYH1",
        "order_7",
        "order_2",
        "order_1",
        "fiber_ID"
    )) |>
    tidyr::pivot_longer(
        cols = c("MYH7", "MYH2", "MYH1"),
        names_to = "MYHs",
        values_to = "values"
    )

write.csv(perc_MYHs,
          here::here("data/perc_MYH_proteomics.csv"))

# Individual curve for intensity of MYH7:

perc_MYHs |>
    dplyr::filter(MYHs == "MYH7") |>
    ggplot2::ggplot(
        ggplot2::aes(order_7,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::theme_minimal()

# MYH7, MYH2 and MYH1 curves all together ranked by MYH7:

perc_MYHs |>
    ggplot2::ggplot(ggplot2::aes(
        x = order_7,
        y = values
    )) +
    ggplot2::geom_point(ggplot2::aes(colour = MYHs)) +
    ggplot2::scale_color_manual(values = c(
        "#fdc325",
        "#5DC863FF",
        "#440154FF"
    )) +
    ggplot2::theme_minimal()


# Finding bottom knee of MYH 7 curve,
# that's the one we are using as a threshold.
# The function detects top knee so we do the inverse of the curve:

MYH_7_curve <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH7) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    t()

MYH_7_curve <- DropletUtils::barcodeRanks(MYH_7_curve, lower = 10)

bottom_knee_MYH7 <- 100 - MYH_7_curve@metadata$knee

perc_MYHs |>
    dplyr::filter(MYHs == "MYH7") |>
    ggplot2::ggplot(
        ggplot2::aes(order_7,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = "#440154FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH7, color = "red") +
    ggplot2::theme_minimal()

# Same with MYH2:

MYH_2_curve <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH2) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    t()

MYH_2_curve <- DropletUtils::barcodeRanks(MYH_2_curve, lower = 10)

bottom_knee_MYH2 <- 100 - MYH_2_curve@metadata$knee

perc_MYHs |>
    dplyr::filter(MYHs == "MYH2") |>
    ggplot2::ggplot(
        ggplot2::aes(order_2,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = "#5DC863FF") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH2, color = "red") +
    ggplot2::theme_minimal()

# Now MYH1

MYH_1_curve <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    tibble::column_to_rownames("fiber_ID") |>
    dplyr::select(MYH1) |>
    dplyr::mutate(across(everything(), ~ 100 - .x)) |>
    t()

MYH_1_curve <- DropletUtils::barcodeRanks(MYH_1_curve, lower = 10)

bottom_knee_MYH1 <- 100 - MYH_1_curve@metadata$knee

perc_MYHs |>
    dplyr::filter(MYHs == "MYH1") |>
    ggplot2::ggplot(
        ggplot2::aes(order_1,
                     values,
                     color = MYHs
        )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = "#fdc325") +
    ggplot2::geom_hline(yintercept = bottom_knee_MYH1, color = "red") +
    ggplot2::theme_minimal()

# This is the case when conditional that I have used to assign fiber types:

data_fiber_type <- perc_MYHs |>
    tidyr::pivot_wider(
        names_from = MYHs,
        values_from = values
    ) |>
    dplyr::mutate(
        fiber_type = dplyr::case_when(
            MYH7 >= bottom_knee_MYH7 &
                MYH2 <= bottom_knee_MYH2 ~ "Type 1",
            MYH7 >= bottom_knee_MYH7 &
                MYH2 >= bottom_knee_MYH2 ~ "Hybrid 1/2A",
            MYH1 >= bottom_knee_MYH1 &
                MYH2 >= bottom_knee_MYH2 ~ "Hybrid 2A/2X",
            MYH2 >= bottom_knee_MYH2 &
                MYH7 <= bottom_knee_MYH7 ~ "Type 2A",
            TRUE ~ "Hybrid 1/2A"
        )
    ) |>
    dplyr::arrange(desc(fiber_ID))

# Again, how the % of MYH look like themselves

data_fiber_type |>
    tidyr::pivot_longer(
        cols = c("MYH7", "MYH2", "MYH1"),
        names_to = "MYHs",
        values_to = "values"
    ) |>
    ggplot2::ggplot(ggplot2::aes(
        x = order_7,
        y = values
    )) +
    ggplot2::geom_point(ggplot2::aes(colour = MYHs)) +
    ggplot2::scale_color_manual(values = c(
        "#fdc325",
        "#5DC863FF",
        "#440154FF"
    )) +
    ggplot2::theme_minimal()

# When I colour by fiber type:

data_fiber_type |>
    tidyr::pivot_longer(
        cols = c("MYH7", "MYH2", "MYH1"),
        names_to = "MYHs",
        values_to = "values"
    ) |>
    ggplot2::ggplot(ggplot2::aes(
        x = order_7,
        y = values
    )) +
    ggplot2::geom_point(ggplot2::aes(colour = fiber_type)) +
    ggplot2::scale_color_manual(values = c(
        "#3B528BFF",
        "#fdc325",
        "#440154FF",
        "#5DC863FF"
    )) +
    ggplot2::theme_minimal()

# Overall number of fibers in each fiber type:

data_fiber_type |>
    ggplot2::ggplot(ggplot2::aes(
        x = fiber_type,
        fill = fiber_type
    )) +
    ggplot2::geom_bar(na.rm = TRUE) +
    ggplot2::scale_fill_manual("Fiber types",
                               values = c(
                                   "#3B528BFF",
                                   "#fdc325",
                                   "#440154FF",
                                   "#5DC863FF"
                               )
    ) +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Fiber types") +
    ggplot2::theme(plot.title = ggplot2::element_text(
        size = 20,
        face = "bold"
    )) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.title = ggplot2::element_text(
        size = 16,
        face = "bold"
    )) +
    ggplot2::theme(
        axis.title.x = ggplot2::element_text(vjust = -0.35),
        axis.title.y = ggplot2::element_text(vjust = 0.35)
    ) +
    ggplot2::labs(
        x = "Fiber Type",
        y = "Count"
    )

write.csv(data_fiber_type, here::here("data/data_fiber_type_saved.csv"))

# Fiber_typing_curves_MYHs ------------------------------------------------

# Here I'm basing my final plot of fiber typing curves on Thibaux script:

fiber_types <- list(
    type_1 <- data_fiber_type |>
        dplyr::filter(fiber_type == "Type 1") |>
        dplyr::arrange(desc(MYH7)),
    Hybrid_1_2A <- data_fiber_type |>
        dplyr::filter(fiber_type == "Hybrid 1/2A") |>
        dplyr::arrange(desc(MYH7)),
    type_2A <- data_fiber_type |>
        dplyr::filter(fiber_type == "Type 2A") |>
        dplyr::arrange(MYH2),
    Hybrid_2A_2X <- data_fiber_type |>
        dplyr::filter(fiber_type == "Hybrid 2A/2X") |>
        dplyr::arrange(MYH1)
)
names(fiber_types) <- c(
    "type_1",
    "hybrid_1_2A",
    "type_2A",
    "hybrid_2A_2X"
)

order_matrix <- dplyr::bind_rows(
    fiber_types$type_1,
    fiber_types$hybrid_1_2A,
    fiber_types$type_2A,
    fiber_types$hybrid_2A_2X
) |>
    tibble::add_column("sample_MYH" = seq_len(nrow(data_fiber_type)))

rectangle_colors <- data.frame(
    start = c(0, 331, 414, 787),
    end = c(330, 413, 786, 974),
    fiber_type = c("Type I", "Hybrid I/IIA", "Type IIA", "Hybrid IIA/IIX")
)

################################################################################################################################################
########################################################      PANEL F   ############################################################################
################################################################################################################################################

order_matrix |>
    ggplot2::ggplot() +
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
        "#5DC863FF"
    )) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_point(
        ggplot2::aes(
            x = sample_MYH,
            y = MYH7
        ),
        colour = "#440154FF",
        size = 0.25
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
            x = sample_MYH,
            y = MYH2
        ),
        colour = "#5DC863FF",
        size = 0.25
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
            x = sample_MYH,
            y = MYH1
        ),
        colour = "#fdc325",
        size = 0.25
    ) +
    ggplot2::ylab("% MYH isoform expressed") +
    ggplot2::xlab("Sample (ranked by MYH expression)") +
    ggplot2::ggtitle("MYH-based fiber typing (proteomics)") +
    ggplot2::geom_vline(xintercept = 330, colour = "grey30", size = 0.25) +
    ggplot2::geom_vline(xintercept = 413, colour = "grey30", size = 0.25) +
    ggplot2::geom_vline(xintercept = 786, colour = "grey30", size = 0.25) +
    ggplot2::geom_hline(
        yintercept = 31,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = 8.9,
        linetype = "dotted",
        colour = "black",
        size = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = 7.9,
        linetype = "dotted",
        colour = "black",
        size = 0.25
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
        x = 165,
        y = 115,
        label = "Type 1 \n 330 (33.9%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 599,
        y = 115,
        label = "Type 2A \n 373 (38.3%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 880,
        y = 115,
        label = "Hybrid 2A/2X \n 188 (19.3%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 371.5,
        y = 110,
        label = "Hybrid \n1/2A \n83 (8.5%)",
        colour = "black",
        fontface = 2,
        size = 1.7
    ) +
    ggplot2::annotate(
        "text",
        x = 35,
        y = 95,
        label = "MYH7",
        colour = "#440154FF",
        fontface = 2,
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 750,
        y = 85,
        label = "MYH2",
        colour = "#5DC863FF",
        fontface = 2,
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 950,
        y = 65,
        label = "MYH1",
        colour = "#fdc325",
        fontface = 2,
        size = 2
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 974), clip = "off") +
    ggplot2::annotate(
        "text",
        x = 75,
        y = 35,
        label = "MYH7 threshold",
        colour = "black",
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 75,
        y = 13,
        label = "MYH2 threshold",
        colour = "black",
        size = 2
    ) +
    ggplot2::annotate(
        "text",
        x = 75,
        y = 5.5,
        label = "MYH1 threshold",
        colour = "black",
        size = 2
    ) +
    ggplot2::theme(text = ggplot2::element_text(size = 7)) +
    ggplot2::scale_x_continuous(limits = c(0, 974),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0,125),
                                expand = c(0, 0),
                                breaks = c(0, 25, 50, 75, 100),
                                labels = c("0", "25", "50", "75", "100")) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 7))

ggplot2::ggsave(
    here::here("doc/figures/figure_1/fiber_type_curves_proteomics.png"),
    device = "png",
    width = 128,
    height = 60,
    units = "mm"
)


# Fiber typing curves transcriptomics -------------------------------------



# Transcriptomics

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


# UMAPS -------------------------------------------------------------------

# Filtered and processed data
data_proteomics <- read.csv(here::here("data/proteomics_pca_data.csv")) # 974 fibers for 1685 proteins

data_proteomics <- data_proteomics |>
    dplyr::rename("Protein" = 1) |>
    tibble::column_to_rownames("Protein")

metadata <- vroom::vroom(
    here::here("data/metadata_proteomics_fiber_type.csv")
) |>
    dplyr::select(!1)

seurat_proteome <- Seurat::CreateSeuratObject(counts = data_proteomics,
                                              meta.data = metadata)

# PCA ---------------------------------------------------------------------


# Find Variable features
seurat_proteome <- Seurat::FindVariableFeatures(seurat_proteome,
                                                selection.method = "vst")

# Scale data
seurat_proteome <- Seurat::ScaleData(seurat_proteome)

#  Run PCA------------------------------------------------
seurat_proteome <- Seurat::RunPCA(object = seurat_proteome,  features = Seurat::VariableFeatures(object = seurat_proteome))

# Explore PCA plots
PCA_PC1_2 <- Seurat::DimPlot(seurat_proteome, reduction = "pca", dims = c(1,2), group.by = "fiber_type") # No separation of fiber type by PC1, good separation by PC2
PCA_PC1_3 <- Seurat::DimPlot(seurat_proteome, reduction = "pca", dims = c(1,3), group.by = "fiber_type") # No separation at all by PC3

# Identify significant PCS
# Elbow plot: visualizes SD of each PC, search for where SD begins to plateau --------
elbow_plot <- Seurat::ElbowPlot(object = seurat_proteome,
                                ndims = 40) # Biggest difference after 6 PCs, plateau only after about 30 PCs

# Visual inspection: 40 PCs

# Quantitative determination elbow plateau (two metrics, choose lowest value)--------------------------------

# Metric 1: New PC only contributes 5%, and all other cumulatively contribute 90%
pct <- seurat_proteome[["pca"]]@stdev / sum(seurat_proteome[["pca"]]@stdev) * 100 # Calculate percent of variation for each PC
cumu <- cumsum(pct) # Calculate cumulative percents with each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 # PC 41

# Metric 2: PC where percent change to next PC is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
co2 # PC 23

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Conclusion: use first 6 PCs to cluster (number less important with new versions of Seurat, biggest difference is time to compute with more PCs)

# Graph-based clustering using K-nearest neighbor graph  ---------------------------------

# Determine the K-nearest neighbor graph (dims is the selected number of PCs from previous step)
seurat_proteome <- Seurat::FindNeighbors(object = seurat_proteome,  dims = 1:6)

# Determine the clusters for various resolutions (resolution between 0.4-1.4 is often best for scRNAseq --> determine which resolution is best for our dataset)
seurat_proteome <- Seurat::FindClusters(object = seurat_proteome, resolution = c(0.4))

# Dimensionality reduction

# Run UMAP ----------------------------------------------------------------
seurat_proteome <- Seurat::RunUMAP(seurat_proteome, dims = 1:6, seed.use = 42)

seurat_metadata <- seurat_proteome@meta.data |>
    dplyr::mutate(
        seurat_clusters = dplyr::case_when(
            RNA_snn_res.0.4 %in% c(1, 3) ~ "slow",
            TRUE ~ "fast"
        )
    )

write.csv(seurat_metadata,
          here::here("data/metadata_proteomics_seurat_clusters.csv"))

# Run T-SNE ----------------------------------------------------------------
seurat_proteome <- Seurat::RunTSNE(seurat_proteome, dims = 1:6)



################################################################################################################################################
#################################################     Panel H  ##############################################################
################################################################################################################################################


my_cols <- c("#440154FF", "#8CB3E8", "#5DC863FF", "#fdc325")

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 5,
    pt.size = 0.3,
    cols = ggplot2::alpha(my_cols, 1),
    group.by = "fiber_type",
    order = c("Hybrid 2A/2X",
              "Type 2A",
              "Hybrid 1/2A",
              "Type 1")) +
    # scale_color_viridis_d(option = "magma")
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle("UMAP Proteomics - by fiber type") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        axis.text = ggplot2::element_text(size=6),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.text= ggplot2::element_text(size=4),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position = "none"
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_1/UMAP_proteomics_fiber_type.png"),
#     height = 55,
#     width = 65,
#     units = "mm"
# )

################################################################################################################################################
########################################################      PANEL J   ############################################################################
################################################################################################################################################

# Feature plots -----------------------------------------------------------

feature_plot_MYH2 <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("MYH2"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot_MYH7 <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("MYH7"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")


feature_plot_MYH1 <- Seurat::FeaturePlot(seurat_proteome,
                                         features = c("MYH1"),
                                         pt.size = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=4),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size=5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(2, 'mm'),
        legend.spacing.x = ggplot2::unit(2, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
        plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = ggplot2::element_blank(), #remove major gridlines
        panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
        legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = ggplot2::element_rect(fill='transparent') #transparent legend panel
    ) +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")

feature_plot <- ggpubr::ggarrange(feature_plot_MYH7,
                                  feature_plot_MYH2,
                                  feature_plot_MYH1,
                                  ncol = 3,
                                  nrow = 1,
                                  legend = "none") +
    ggplot2::ggtitle("Proteomics") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 8,
                                           face = "bold"),
    )

ggplot2::ggsave(here::here("doc/figures/umaps_now_in_fig1/feature_plot_MYHs.png"),
                width = 100,
                height = 33,
                units = "mm")

# UMAP subject ------------------------------------------------------------

my_cols <- viridisLite::turbo(n = 5)

Seurat::DimPlot(
    seurat_proteome,
    reduction = "umap",
    label = FALSE,
    label.size = 5,
    pt.size = 0.5,
    cols = ggplot2::alpha(my_cols, 0.65),
    group.by = "subject") +
    ggplot2::ggtitle("Resolution 0.4") +
    # scale_color_viridis_d(option = "magma")
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle("UMAP Proteomics - by subject") +
    ggplot2::theme(
        text = ggplot2::element_text(face="bold", colour="black", size=6),
        axis.text = ggplot2::element_text(size=6),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.text= ggplot2::element_text(size=4),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position = "right"
    )

# ggplot2::ggsave(
#     here::here("doc/figures/figure_1/UMAP_proteomics_subject.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )
