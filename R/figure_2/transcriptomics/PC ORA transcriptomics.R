################################################################################################################################################
#################################################       USED LINKS      ########################################################################
################################################################################################################################################




################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggh4x)


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")


################################################################################################################################################
###########################################       EXTRACT GENES OF PCs      ####################################################################
################################################################################################################################################

# Load PC gene loading list
PC_genes <- read_csv("~/single_fiber_heterogeneity/data/transcriptomics_PC_loadings.csv")

# Rename first column
colnames(PC_genes)[1] = "Gene"

# 5% of genes
perc5 <- nrow(PC_genes) * 0.05

# Extract 5% top genes PC1 negative
PC1_neg <- PC_genes %>% arrange(PC_1) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC1 positive
PC1_pos <- PC_genes %>% arrange(desc(PC_1)) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC1 all
PC1_all <- PC_genes %>% arrange(desc(abs(PC_1))) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC2 negative
PC2_neg <- PC_genes %>% arrange(PC_2) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC2 positive
PC2_pos <- PC_genes %>% arrange(desc(PC_2)) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC2 all
PC2_all <- PC_genes %>% arrange(desc(abs(PC_2))) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC3 negative
PC3_neg <- PC_genes %>% arrange(PC_3) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC3 positive
PC3_pos <- PC_genes %>% arrange(desc(PC_3)) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC3 all
PC3_all <- PC_genes %>% arrange(desc(abs(PC_3))) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC4 negative
PC4_neg <- PC_genes %>% arrange(PC_4) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC4 positive
PC4_pos <- PC_genes %>% arrange(desc(PC_4)) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)

# Extract  5% top genes PC4 all
PC4_all <- PC_genes %>% arrange(desc(abs(PC_4))) %>% dplyr::select(Gene) %>% dplyr::slice(1:perc5)


################################################################################################################################################
#######################################       EXTRACT ALL BACKGROUND GENES      ################################################################
################################################################################################################################################

allgenes_SCT <- rownames(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@assays$SCT@counts)

################################################################################################################################################
##############################################       ORA PC1 - POS      #########################################################################
################################################################################################################################################

# Perform ORA with all background genes SCT -------------------------------------------------------------
ORA_PC1_pos <- enrichGO(gene = PC1_pos$Gene,
                       universe      = allgenes_SCT,
                       keyType       = "SYMBOL",
                       OrgDb         = org.Hs.eg.db,
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = FALSE)


# Simplify GO terms redundancy --------------------------------------------
ORA_PC1_pos_simplify <- simplify(ORA_PC1_pos, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
ORA_PC1_pos_simplify_df <- as.data.frame(ORA_PC1_pos_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
ORA_PC1_pos_simplify_df <- mutate(ORA_PC1_pos_simplify_df, foldEnrich =
                 (as.numeric(sub("/\\d+", "", ORA_PC1_pos_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_PC1_pos_simplify_df$GeneRatio))) /
                 (as.numeric(sub("/\\d+", "", ORA_PC1_pos_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_PC1_pos_simplify_df$BgRatio)))
               )

# Save results ------------------------------------------------------------
write_csv(ORA_PC1_pos_simplify_df, file="~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC1 pos.csv")

################################################################################################################################################
##############################################       ORA PC1 - NEG      #########################################################################
################################################################################################################################################

# Perform ORA with all background genes SCT -------------------------------------------------------------
ORA_PC1_neg <- enrichGO(gene = PC1_neg$Gene,
                        universe      = allgenes_SCT,
                        keyType       = "SYMBOL",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)


# Simplify GO terms redundancy --------------------------------------------
ORA_PC1_neg_simplify <- simplify(ORA_PC1_neg, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
ORA_PC1_neg_simplify_df <- as.data.frame(ORA_PC1_neg_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
ORA_PC1_neg_simplify_df <- mutate(ORA_PC1_neg_simplify_df, foldEnrich =
                                    (as.numeric(sub("/\\d+", "", ORA_PC1_neg_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_PC1_neg_simplify_df$GeneRatio))) /
                                    (as.numeric(sub("/\\d+", "", ORA_PC1_neg_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_PC1_neg_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(ORA_PC1_neg_simplify_df, file="~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC1 neg.csv")



################################################################################################################################################
##############################################       ORA PC2 - POS      #########################################################################
################################################################################################################################################

# Perform ORA with all background genes SCT -------------------------------------------------------------
ORA_PC2_pos <- enrichGO(gene = PC2_pos$Gene,
                        universe      = allgenes_SCT,
                        keyType       = "SYMBOL",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)


# Simplify GO terms redundancy --------------------------------------------
ORA_PC2_pos_simplify <- simplify(ORA_PC2_pos, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
ORA_PC2_pos_simplify_df <- as.data.frame(ORA_PC2_pos_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
ORA_PC2_pos_simplify_df <- mutate(ORA_PC2_pos_simplify_df, foldEnrich =
                                      (as.numeric(sub("/\\d+", "", ORA_PC2_pos_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_PC2_pos_simplify_df$GeneRatio))) /
                                      (as.numeric(sub("/\\d+", "", ORA_PC2_pos_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_PC2_pos_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(ORA_PC2_pos_simplify_df, file="~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC2 pos.csv")

################################################################################################################################################
##############################################       ORA PC2 - NEG      #########################################################################
################################################################################################################################################

# Perform ORA with all background genes SCT -------------------------------------------------------------
ORA_PC2_neg <- enrichGO(gene = PC2_neg$Gene,
                        universe      = allgenes_SCT,
                        keyType       = "SYMBOL",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)


# Simplify GO terms redundancy --------------------------------------------
ORA_PC2_neg_simplify <- simplify(ORA_PC2_neg, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
ORA_PC2_neg_simplify_df <- as.data.frame(ORA_PC2_neg_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
ORA_PC2_neg_simplify_df <- mutate(ORA_PC2_neg_simplify_df, foldEnrich =
                                      (as.numeric(sub("/\\d+", "", ORA_PC2_neg_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", ORA_PC2_neg_simplify_df$GeneRatio))) /
                                      (as.numeric(sub("/\\d+", "", ORA_PC2_neg_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", ORA_PC2_neg_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(ORA_PC2_neg_simplify_df, file="~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC2 neg.csv")


################################################################################################################################################
###############################################      DOT PLOT Tx AND Px COMBINED      ##########################################################
################################################################################################################################################

Tx_PC1_pos <- read.csv("~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC1 pos.csv")
Tx_PC1_neg <- read.csv("~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC1 neg.csv")
Tx_PC2_pos <- read.csv("~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC2 pos.csv")
Tx_PC2_neg <- read.csv("~/single_fiber_heterogeneity/data/PCA_transcriptomics/PC2 neg.csv")

Px_PC1_pos <- read.csv("~/single_fiber_heterogeneity/data/GSEA_PCA_proteomics/results_PC1_positive.csv")
Px_PC1_neg <- read.csv("~/single_fiber_heterogeneity/data/GSEA_PCA_proteomics/results_PC1_negative.csv")
Px_PC2_pos <- read.csv("~/single_fiber_heterogeneity/data/GSEA_PCA_proteomics/results_PC2_positive.csv")
Px_PC2_neg <- read.csv("~/single_fiber_heterogeneity/data/GSEA_PCA_proteomics/results_PC2_negative.csv")

# Calculate foldEnrich for proteomics
Px_PC1_pos <- mutate(Px_PC1_pos, foldEnrich =
                                    (as.numeric(sub("/\\d+", "", Px_PC1_pos$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC1_pos$GeneRatio))) /
                                    (as.numeric(sub("/\\d+", "", Px_PC1_pos$BgRatio)) / as.numeric(sub(".*/", "", Px_PC1_pos$BgRatio))))

Px_PC1_neg <- mutate(Px_PC1_neg, foldEnrich =
                       (as.numeric(sub("/\\d+", "", Px_PC1_neg$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC1_neg$GeneRatio))) /
                       (as.numeric(sub("/\\d+", "", Px_PC1_neg$BgRatio)) / as.numeric(sub(".*/", "", Px_PC1_neg$BgRatio))))

Px_PC2_pos <- mutate(Px_PC2_pos, foldEnrich =
                       (as.numeric(sub("/\\d+", "", Px_PC2_pos$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC2_pos$GeneRatio))) /
                       (as.numeric(sub("/\\d+", "", Px_PC2_pos$BgRatio)) / as.numeric(sub(".*/", "", Px_PC2_pos$BgRatio))))

Px_PC2_neg <- mutate(Px_PC2_neg, foldEnrich =
                       (as.numeric(sub("/\\d+", "", Px_PC2_neg$GeneRatio)) / as.numeric(sub(".*/", "", Px_PC2_neg$GeneRatio))) /
                       (as.numeric(sub("/\\d+", "", Px_PC2_neg$BgRatio)) / as.numeric(sub(".*/", "", Px_PC2_neg$BgRatio))))


# Select categories for PC1 positive - Tx and Px
Px_PC1_pos_plot <- Px_PC1_pos %>%
  dplyr::filter(
    Description == "cytosolic ribosome" |
      Description == "cytoplasmic translation" |
      Description == "structural constituent of ribosome" |
      Description == "mitochondrial respirasome"
  )

Tx_PC1_pos_plot <- Tx_PC1_pos %>%
  dplyr::filter(
    Description == "cytosolic ribosome" |
      Description == "cytoplasmic translation" |
      Description == "structural constituent of ribosome" |
      Description == "mitochondrial respirasome"
  )

# Select categories for PC1 negative - Tx and Px
Px_PC1_neg_plot <- Px_PC1_neg %>%
  dplyr::filter(
    Description == "cell junction" |
      Description == "costamere" |
      Description == "enzyme binding"
  )

Tx_PC1_neg_plot <- Tx_PC1_neg %>%
  dplyr::filter(
    Description == "GTPase regulator activity" |
      Description == "cell-cell adhesion" |
      Description == "cell-cell junction"
  )

# Select categories for PC2 positive - Tx and Px
Px_PC2_pos_plot <- Px_PC2_pos %>%
  dplyr::filter(
    Description == "fatty-acyl-CoA binding" |
      Description == "oxidoreductase activity" |
      Description == "contractile fiber"
  )

Tx_PC2_pos_plot <- Tx_PC2_neg %>%
  dplyr::filter(
    Description == "fatty-acyl-CoA binding" |
      Description == "oxidoreductase activity" |
      Description == "contractile fiber"
  )

# Select categories for PC2 negative - Tx and Px
Px_PC2_neg_plot <- Px_PC2_neg %>%
  dplyr::filter(
    Description == "carbohydrate binding" |
      Description == "nucleotide phosphorylation" |
      Description == "contractile fiber"
  )

Tx_PC2_neg_plot <- Tx_PC2_pos %>%
  dplyr::filter(
    Description == "carbohydrate binding" |
      Description == "nucleotide phosphorylation" |
      Description == "contractile fiber"
  )

# Add required columns for plotting
Px_PC1_pos_plot$PCDir <- rep("PC1 - Pos", nrow(Px_PC1_pos_plot))
Px_PC1_pos_plot$graphdir <- rep(3, nrow(Px_PC1_pos_plot))
Px_PC1_pos_plot$data <- rep("Px", nrow(Px_PC1_pos_plot))

Tx_PC1_pos_plot$PCDir <- rep("PC1 - Pos", nrow(Tx_PC1_pos_plot))
Tx_PC1_pos_plot$graphdir <- rep(4, nrow(Tx_PC1_pos_plot))
Tx_PC1_pos_plot$data <- rep("Tx", nrow(Tx_PC1_pos_plot))

Px_PC1_neg_plot$PCDir <- rep("PC1 - Neg", nrow(Px_PC1_neg_plot))
Px_PC1_neg_plot$graphdir <- rep(1, nrow(Px_PC1_neg_plot))
Px_PC1_neg_plot$data <- rep("Px", nrow(Px_PC1_neg_plot))

Tx_PC1_neg_plot$PCDir <- rep("PC1 - Neg", nrow(Tx_PC1_neg_plot))
Tx_PC1_neg_plot$graphdir <- rep(2, nrow(Tx_PC1_neg_plot))
Tx_PC1_neg_plot$data <- rep("Tx", nrow(Tx_PC1_neg_plot))

Px_PC2_pos_plot$PCDir <- rep("PC2 - Pos", nrow(Px_PC2_pos_plot))
Px_PC2_pos_plot$graphdir <- rep(7, nrow(Px_PC2_pos_plot))
Px_PC2_pos_plot$data <- rep("Px", nrow(Px_PC2_pos_plot))

Tx_PC2_pos_plot$PCDir <- rep("PC2 - Pos", nrow(Tx_PC2_pos_plot))
Tx_PC2_pos_plot$graphdir <- rep(8, nrow(Tx_PC2_pos_plot))
Tx_PC2_pos_plot$data <- rep("Tx", nrow(Tx_PC2_pos_plot))

Px_PC2_neg_plot$PCDir <- rep("PC2 - Neg", nrow(Px_PC2_neg_plot))
Px_PC2_neg_plot$graphdir <- rep(5, nrow(Px_PC2_neg_plot))
Px_PC2_neg_plot$data <- rep("Px", nrow(Px_PC2_neg_plot))

Tx_PC2_neg_plot$PCDir <- rep("PC2 - Neg", nrow(Tx_PC2_neg_plot))
Tx_PC2_neg_plot$graphdir <- rep(6, nrow(Tx_PC2_neg_plot))
Tx_PC2_neg_plot$data <- rep("Tx", nrow(Tx_PC2_neg_plot))

# Combine all into one dataframe
PC_plot <- bind_rows(Tx_PC2_pos_plot, Px_PC2_pos_plot, Tx_PC2_neg_plot, Px_PC2_neg_plot, Tx_PC1_pos_plot, Px_PC1_pos_plot, Tx_PC1_neg_plot, Px_PC1_neg_plot)
PC_plot$graph <- rep(1, nrow(PC_plot))
PC_plot$PCDir = factor(PC_plot$PCDir, levels=c("PC2 - Pos", "PC2 - Neg", "PC1 - Pos", "PC1 - Neg"))
PC_plot$data = factor(PC_plot$data, levels=c("Tx", "Px"))

# Create dot plot
dotplot_combined <- ggplot(PC_plot, aes(graph, fct_reorder(Description, graphdir))) +
  geom_point(aes(color= data, alpha=foldEnrich, size = -log10(p.adjust))) +
  scale_color_manual(values = c("#354D3C", "#045a8d", "#354D3C", "#045a8d", "#354D3C","#045a8d", "#354D3C", "#045a8d")) +
  scale_alpha_continuous(
      name = "Fold enrichment",
      limits = c(2.0004, 12.574),
      breaks = c(2.0004, 12.574),
      range = c(0.2, 1)) +
  scale_size_continuous(
      name = "-log10(p.adjust)",
      limits = c(1.430, 33.73),
      breaks = c(1.430, 33.73),
      range=c(1.5, 5)) +
  ggh4x::facet_nested_wrap(~PCDir + data , nrow = 1, ncol=8) +
  theme_classic() +
  scale_x_continuous(limits = c(0.9,1.1), expand = c(0, 0)) +
  theme(
    text = ggplot2::element_text(face = "bold",size = 9, colour = "black"),
    axis.title.y= element_blank(),
    axis.title.x= element_blank(),
    axis.text.x= element_blank(),
    axis.line.x =element_blank(),
    axis.ticks.x =element_blank(),
    #axis.text.y = element_text(size=8),
    plot.title = element_text(hjust = 0.5, face="bold", vjust=-7),
    legend.position = "none",
    strip.background = element_rect(colour=NA, fill="white"),
    strip.text = element_text(colour = "white")
  ) +
  coord_fixed(ratio = 0.2)

ggsave(dotplot_combined,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/Combined_dot_plot.png",
       width = 128,
       height = 128,
       units="mm")

# small version:

dotplot_combined_small <- ggplot(PC_plot, aes(graph, fct_reorder(Description, graphdir))) +
    geom_point(aes(color= data, alpha=foldEnrich, size = -log10(p.adjust))) +
    scale_color_manual(values = c("#354D3C", "#045a8d", "#354D3C", "#045a8d", "#354D3C","#045a8d", "#354D3C", "#045a8d")) +
    scale_alpha_continuous(
        name = "Fold enrichment",
        limits = c(2.0004, 12.574),
        breaks = c(2.0004, 12.574),
        range = c(0.2, 1)) +
    scale_size_continuous(
        name = "-log10(p.adjust)",
        limits = c(1.430, 33.73),
        breaks = c(1.430, 33.73),
        range=c(1, 3.5)) +
    ggh4x::facet_nested_wrap(~PCDir + data , nrow = 1, ncol=8) +
    theme_classic() +
    scale_x_continuous(limits = c(0.9,1.1), expand = c(0, 0)) +
    theme(
        text = ggplot2::element_text(face = "bold",size = 8, colour = "black"),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x= element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x =element_blank(),
        #axis.text.y = element_text(size=8),
        # plot.title = element_text(hjust = 0.5, face="bold", vjust=-7),
        plot.title = element_text(hjust = 1, face="bold"),
        legend.position = "none",
        strip.background = element_rect(colour=NA, fill="white"),
        strip.text = element_text(colour = "white")
    )
# coord_fixed(ratio = 0.2)

ggsave(dotplot_combined_small,
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_4/Combined_dot_plot_small.png",
       width = 80,
       height = 100,
       units="mm")

# Export data frame to make legend in Illustrator
write.csv(PC_plot, file = "~/single_fiber_heterogeneity/data/combined_dotplot/PCA_enrichment_data_plot.csv")
