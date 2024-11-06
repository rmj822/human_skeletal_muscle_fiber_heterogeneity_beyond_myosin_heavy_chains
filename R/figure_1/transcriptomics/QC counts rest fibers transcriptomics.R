
################################################################################################################################################
#################################################       USED LINKS      ####@###################################################################
################################################################################################################################################

# https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_quality_control_analysis.html
# https://github.com/hbctraining/scRNA-seq/tree/master/schedule

################################################################################################################################################
################################################       PREPARATION      ####@###################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Matrix)
library(Seurat)
library(viridis)
library(rstatix)
library(readxl)

# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")


################################################################################################################################################
#############################################       LOAD SEURAT OBEJCT      ####@###############################################################
################################################################################################################################################
load("4 Seurat object creation/seurat_wo_MSTRG.RData")


################################################################################################################################################
########################################       FILTER ONLY FIBERS AT REST      #################################################################
################################################################################################################################################

seurat_wo_MSTRG_rest <- subset(seurat_wo_MSTRG, subset = time == "Pre") # 1044 fibers

################################################################################################################################################
########################################      EXTRACT META DATA HETEROFIBER FOR EGA UPLOAD      ################################################
################################################################################################################################################

EGA <- seurat_wo_MSTRG_rest@meta.data
EGA$sample <- rownames(EGA)
EGA <- EGA %>% dplyr::select(sample, subject, fiber, test_day, time, condition)

write.csv(EGA,
          file = "~/single_fiber_heterogeneity/data/EGA_upload/meta_JVS2101_heterofiber.csv")

################################################################################################################################################
########################################      EXTRACT META DATA EXERFIBER FOR EGA UPLOAD      ################################################
################################################################################################################################################

seurat_wo_MSTRG_exercise <- subset(seurat_wo_MSTRG, subset = time != "Pre") #  2096 fibers

EGA_exercise <- seurat_wo_MSTRG_exercise@meta.data
EGA_exercise$sample <- rownames(EGA_exercise)
EGA_exercise <- EGA_exercise %>% dplyr::select(sample, subject, fiber, test_day, time, condition)

write.csv(EGA_exercise,
          file = "~/single_fiber_heterogeneity/data/EGA_upload/meta_JVS2101_exerfiber.csv")

################################################################################################################################################
##############################################       EXTRACT METADATA      #####@###############################################################
################################################################################################################################################

# Make factors for variables ----------------------------------------------
seurat_wo_MSTRG_rest@meta.data$subject <- factor(seurat_wo_MSTRG_rest@meta.data$subject, levels = c(1:14))

seurat_wo_MSTRG_rest@meta.data$test_day <- factor(seurat_wo_MSTRG_rest@meta.data$test_day, levels = c("T1", "T2", "T3"))

seurat_wo_MSTRG_rest@meta.data$condition <- factor(seurat_wo_MSTRG_rest@meta.data$condition, levels = c("P", "H1", "H2"))


################################################################################################################################################
#############################################       ASSESSING QC METRICS      ####@#############################################################
################################################################################################################################################

# Fiber counts per subject -------------------------------------------------------------
pdf("5 QC/Fibers at rest/Counts/Initial metrics/N fibers.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=subject, fill=subject)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_fill_viridis(discrete=T, option="magma") +
  labs(
    x = "Subject",
    y = "N fibers",
    title = "Number of fibers per subject"
  )  +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text = element_text(colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

 # Conclusion: 75 fibers per subject, except for some due to missing samples (1 for each sub6 and sub8) or due to >100 genes filtering

# UMI counts per fiber ----------------------------------------------------

  # scRNAseq: 500-1000 UMIs per cell is low, but useable. Above 1000 better. 500 is cut-off line

# All subjects together
pdf("5 QC/Fibers at rest/Counts/Initial metrics/UMI distribution all.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

  # Each subject seperately
pdf("5 QC/Fibers at rest/Counts/Initial metrics/UMI distribution individual.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill="grey"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
facet_wrap(~ subject, nrow=5)
dev.off()

    # Conclusion: All subjects have similar profile (sub7 possibly worst). N fibers below 500 UMI counts is limited.

# Genes detected per fiber ------------------------------------------------

# All subjects together
pdf("5 QC/Fibers at rest/Counts/Initial metrics/Detected genes distribution all.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nGene, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Each subject seperately
pdf("5 QC/Fibers at rest/Counts/Initial metrics/Detected genes distribution individual.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nGene, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill="grey"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()

    # Conclusion: All subjects have similar profile. N fibers below 200 genes is very limited.


# UMIs vs genes detected --------------------------------------------------

    # Poor quality fibers have low UMIs and low genes = lower left in graph. Good linearity, with similar slope for each subject, expected
    # Low complexity fiber in bottom right: high number of UMI, but low genes. Delete, possible contamination from other cells

# All subjects together
pdf("5 QC/Fibers at rest/Counts/Initial metrics/UMI vs genes all.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, y=nGene, colour= subject)) +
  scale_colour_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  ylab("log10 nGene") +
  xlab("log10 nUMI") +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Each subject seperately
pdf("5 QC/Fibers at rest/Counts/Initial metrics/UMI vs genes individual.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, y=nGene, colour= subject)) +
  scale_colour_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  ylab("log10 nGene") +
  xlab("log10 nUMI") +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()

    # Conclusion: Good linearity, and very similar between subjects. Some fibers in lower-left quadrant

# Complexity based on UMI counts -----------------------------------------------------------------

    # Gives info on gene complexity of RNA content of sample (higher = better).
    # Lower value indicates not saturating seq for any given gene. General good value = 0.80, we will use 0.75

# All subjects together
pdf("5 QC/Fibers at rest/Counts/Initial metrics/UMI complexity all.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 0.75) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Each subject seperately
pdf("5 QC/Fibers at rest/Counts/Initial metrics/UMI complexity individual.pdf")
seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 0.75) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill="grey"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()

    # Conclusion: Relatively high number below 0.80 (maybe because of high counts for sarcomeric genes?). Almost none below 0.75


################################################################################################################################################
###########################################       FILTERING LOW ABUNDANT GENES     ####@########################################################
################################################################################################################################################

# Extract counts ----------------------------------------------------------
counts.seurat <- GetAssayData(object = seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")

# Calculate % of samples that express each gene --------
genes.percent.expression <- rowMeans(counts.seurat>0 ) * 100

# Keep only genes that are expressed in 30% of samples--------
genes.filter <- names(genes.percent.expression[genes.percent.expression>30])   # 7418 genes retained

# Filter counts ----------------
counts.sub <- counts.seurat[genes.filter,]

# Reassign to filtered Seurat object --------------------------------------
filtered_seurat_wo_MSTRG_rest <- CreateSeuratObject(counts.sub)

# Check if order of metadata and counts is identical
all(rownames(seurat_wo_MSTRG_rest@meta.data) %in% rownames(filtered_seurat_wo_MSTRG_rest@meta.data))
all(rownames(seurat_wo_MSTRG_rest@meta.data) == rownames(filtered_seurat_wo_MSTRG_rest@meta.data))

all(colnames(seurat_wo_MSTRG_rest) %in% colnames(filtered_seurat_wo_MSTRG_rest))
all(colnames(seurat_wo_MSTRG_rest) == colnames(filtered_seurat_wo_MSTRG_rest))

# Add metadata ----------------------------------------------

filtered_seurat_wo_MSTRG_rest@meta.data$log10GenesPerUMI <- log10(filtered_seurat_wo_MSTRG_rest$nFeature_RNA) / log10(filtered_seurat_wo_MSTRG_rest$nCount_RNA)

filtered_seurat_wo_MSTRG_rest$mitoRatio <- Seurat::PercentageFeatureSet(object = filtered_seurat_wo_MSTRG_rest, pattern = "^MT-")
filtered_seurat_wo_MSTRG_rest$mitoRatio <- filtered_seurat_wo_MSTRG_rest@meta.data$mitoRatio / 100

filtered_seurat_wo_MSTRG_rest@meta.data  <-filtered_seurat_wo_MSTRG_rest@meta.data %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

filtered_seurat_wo_MSTRG_rest@meta.data$subject <- seurat_wo_MSTRG_rest@meta.data$subject
filtered_seurat_wo_MSTRG_rest@meta.data$fiber <- seurat_wo_MSTRG_rest@meta.data$fiber
filtered_seurat_wo_MSTRG_rest@meta.data$test_day <- seurat_wo_MSTRG_rest@meta.data$test_day
filtered_seurat_wo_MSTRG_rest@meta.data$time <- seurat_wo_MSTRG_rest@meta.data$time
filtered_seurat_wo_MSTRG_rest@meta.data$condition <- seurat_wo_MSTRG_rest@meta.data$condition

filtered_seurat_wo_MSTRG_rest@meta.data

# Factorize metadata ----------------------------------------------
filtered_seurat_wo_MSTRG_rest@meta.data$subject <- factor(filtered_seurat_wo_MSTRG_rest@meta.data$subject, levels = c(1:14))
filtered_seurat_wo_MSTRG_rest@meta.data$fiber <- factor(filtered_seurat_wo_MSTRG_rest@meta.data$fiber, levels = c(1:225))
filtered_seurat_wo_MSTRG_rest@meta.data$test_day <- factor(filtered_seurat_wo_MSTRG_rest@meta.data$test_day, levels = c("T1", "T2", "T3"))
filtered_seurat_wo_MSTRG_rest@meta.data$condition <- factor(filtered_seurat_wo_MSTRG_rest@meta.data$condition, levels = c("P", "H1", "H2"))


################################################################################################################################################
#########################################       FILTERING LOW QUALITY SAMPLES      ####@########################################################
################################################################################################################################################

# Threshold: minimum 1000 UMIs, 1000 genes, complexity >0.75 and mitoRatio < 0.20

filtered_seurat_wo_MSTRG_rest <- subset(x = filtered_seurat_wo_MSTRG_rest,
                          subset= (nUMI >= 1000) &
                            (nGene >= 1000) &
                            (log10GenesPerUMI > 0.75) &
                            (mitoRatio < 0.20))


        # From 1044 samples to 926 (344 samples deleted)


################################################################################################################################################
#############################################       REASSESSING QC METRICS      ####@###########################################################
################################################################################################################################################

# Average genes per sample -------------------------------------------------------------
stats_avg_genes_subject <- filtered_seurat_wo_MSTRG_rest@meta.data %>%
  group_by(subject) %>%
  get_summary_stats(nGene, type = "mean_sd")

stats_avg_genes <- filtered_seurat_wo_MSTRG_rest@meta.data %>%
  get_summary_stats(nGene, type = "mean_sd")

write_csv(stats_avg_genes_subject, file="5 QC/Fibers at rest/Counts/After filtering/Average genes per subject.csv")
write_csv(stats_avg_genes, file="5 QC/Fibers at rest/Counts/After filtering/Average genes all.csv")

# Fiber counts per subject -------------------------------------------------------------
pdf("5 QC/Fibers at rest/Counts/After filtering/N fibers.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=subject, fill=subject)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_fill_viridis(discrete=T, option="magma") +
  labs(
    x = "Subject",
    y = "N fibers",
    title = "Number of fibers per subject"
  )  +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text = element_text(colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()


# UMI counts per fiber ----------------------------------------------------

        # scRNAseq: 500-1000 UMIs per cell is low, but useable. Above 1000 better. 500 is cut-off line

# All subjects together
pdf("5 QC/Fibers at rest/Counts/After filtering/UMI distribution all.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Each subject seperately
pdf("5 QC/Fibers at rest/Counts/After filtering/UMI distribution individual.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill="grey"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()

    # Conclusion: All subjects now > 1000 UMI counts

# Genes detected per fiber ------------------------------------------------

# All subjects together
Detected_genes_all <- filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nGene, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.4) +
  geom_density() +
  scale_x_log10(expand = c(0,0)) +
  ylab("Fiber density") +
  xlab("Detected genes (log10)") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=8, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(Detected_genes_all, filename = "5 QC/Fibers at rest/Counts/After filtering/Detected genes all.png", width = 128, height = 60, units="mm")

# Violins per subject

load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata") # Load final data for the figure plot

genes_per_subject <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data %>%
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

ggsave(genes_per_subject, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/detected_genes_per_subject_transcriptomics.png", width = 128, height = 60, units="mm")


# Each subject seperately: density plot
pdf("5 QC/Fibers at rest/Counts/After filtering/Detected genes distribution individual.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nGene, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  scale_x_log10() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill="grey"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()
      # Conclusion: All subjects now > 500 genes.


# UMIs vs genes detected --------------------------------------------------

      # Poor quality fibers have low UMIs and low genes = lower left in graph. Good linearity, with similar slope for each subject, expected
      # Low complexity fiber in bottom right: high number of UMI, but low genes. Delete, possible contamination from other cells

# All subjects together
pdf("5 QC/Fibers at rest/Counts/After filtering/UMI vs genes all.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, y=nGene, colour= subject)) +
  scale_colour_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  ylab("log10 nGene") +
  xlab("log10 nUMI") +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Each subject seperately
pdf("5 QC/Fibers at rest/Counts/After filtering/UMI vs genes individual.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=nUMI, y=nGene, colour= subject)) +
  scale_colour_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  ylab("log10 nGene") +
  xlab("log10 nUMI") +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 1000) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()

        # Conclusion: Good linearity, and very similar between subjects

# Complexity based on UMI counts -----------------------------------------------------------------

    # Gives info on gene complexity of sample (higher = better). Lower value indicates not saturating seq for any given gene. General good value = 0.80

# All subjects together
pdf("5 QC/Fibers at rest/Counts/After filtering/UMI complexity all.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 0.75) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Each subject seperately
pdf("5 QC/Fibers at rest/Counts/After filtering/UMI complexity individual.pdf")
filtered_seurat_wo_MSTRG_rest@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, fill= subject)) +
  scale_fill_viridis(discrete=T, option="magma", alpha=0.7) +
  geom_density() +
  ylab("log10 fiber density") +
  geom_vline(xintercept = 0.75) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill="grey"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ subject, nrow=5)
dev.off()

      # Conclusion: Relatively high number below 0.80 (maybe because of high counts for sarcomeric genes?)


################################################################################################################################################
###########################################       SAVE FILTER SEURAT OBJECT      ####@##########################################################
################################################################################################################################################

# Save Seurat object ------------------------------------------------------
save(filtered_seurat_wo_MSTRG_rest, file="5 QC/Fibers at rest/seurat_filtered_wo_MSTRG_rest.RData")



# Load final Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data %>%
    get_summary_stats(nGene, type = "mean_sd")


