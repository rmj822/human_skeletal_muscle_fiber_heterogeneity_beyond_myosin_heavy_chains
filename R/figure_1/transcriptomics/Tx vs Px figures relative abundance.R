################################################################################################################################################
################################################       PREPARATION      ####@###################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(rtracklayer)
library(GenomicFeatures)
library(ggrepel)
library(VennDiagram)

BiocManager::install("GenomicFeatures")


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("7 Targeted fiber type markers/Fibers at rest/Inflection/filtered_normalized_fibertype_seurat_wo_MSTRG_rest.RData")

# Load Proteomics data
proteomics <- read.csv("~/single_fiber_heterogeneity/data/data_proteomics_filtered.csv")

################################################################################################################################################
#################################################      EXTRACT FILTERED COUNTS      ######################################################################
################################################################################################################################################

transcriptomics <- GetAssayData(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")

################################################################################################################################################
#################################################      CALCULATE RELATIVE EXPRESSION ALL GENES      ######################################################################
################################################################################################################################################

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

# Add genes as rownames
rownames(proteomics_plot) <- proteomics_plot$Gene.name

# Remove unnecessary columns
proteomics_plot <- proteomics_plot %>% dplyr::select(-X) %>% dplyr::select(-Gene.name)

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

################################################################################################################################################
#################################################      CREATE DF WITH RAW INTENSITIES/COUNTS      ##############################################
################################################################################################################################################

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
       filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/detected-genes_transcriptomics-vs-proteomics.png",
       width = 60,
       height =60,
       units="mm")



################################################################################################################################################
####################################     CREATE PLOT WITH RAW INTENSITY/COUNTS      ###################################
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

ggsave(TvsP_plot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/correlation_transcriptomics_vs_proteomics.png", width = 60, height = 60, units="mm")

