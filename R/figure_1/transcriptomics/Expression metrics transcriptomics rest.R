################################################################################################################################################
################################################       PREPARATION      ####@###################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
#library(scater)
library(rtracklayer)
library(GenomicFeatures)
library(ggrepel)

BiocManager::install("GenomicFeatures")


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("7 Targeted fiber type markers/Fibers at rest/Inflection/filtered_normalized_fibertype_seurat_wo_MSTRG_rest.RData")

# Load bulk data
setwd("/Users/thibauxvds/Library/CloudStorage/OneDrive-UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Bulk muscle RNAseq/Own data")

counts_bulk <- read.table("2 Data preparation/counts.txt", header = TRUE, sep = "\t", dec = ".")
meta_bulk <- read.table("2 Data preparation/meta.txt", header = TRUE, sep = "\t", dec = ".")

txdb <- makeTxDbFromGFF("2 Data preparation/TotalRNA_polyA_transcriptome_filtered.gtf",format="gtf")

setwd("~/Library/CloudStorage/OneDrive-UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

################################################################################################################################################
#################################################      EXTRACT FILTERED COUNTS      ######################################################################
################################################################################################################################################

counts <- GetAssayData(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")

################################################################################################################################################
#########################################       EXTRACT NON-FILTERED COUNTS      ####@########################################################
################################################################################################################################################

# Goal: only do sample filtering, but not gene filtering to keep potential non-muscle genes

# Load non-filtered Seurat object -----------------------------------------
load("4 Seurat object creation/seurat_wo_MSTRG.RData")

# Only use fibers at rest (not exercise) -----------------------------------------
seurat_nonfiltered <- subset(seurat_wo_MSTRG, subset = time == "Pre") # 1045 fibers

# Filter for only fibers that passed criteria in original QC -----------------------------------------

# list the cells (=colnames) that you want to keep:
keep <- colnames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest)

# Perform filtering:
seurat_nonfiltered <- seurat_nonfiltered[,colnames(seurat_nonfiltered) %in% keep]

# Check that number of fibers is the same
all(colnames(seurat_nonfiltered) %in% colnames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest))
all(colnames(seurat_nonfiltered) == colnames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest))

# Extract raw counts -----------------------------------------
counts_nonfiltered <- GetAssayData(object = seurat_nonfiltered, assay = "RNA", slot = "counts")

################################################################################################################################################
#########################################       BULK DATA DATA WRANGLING      ####@########################################################
################################################################################################################################################

# Remove MSTRG artefact genes from counts ------------------------------------------------------------
counts_bulk <- counts_bulk[!grepl('MSTRG',rownames(counts_bulk)),]

# Filter genes  ------------------------------------------------------------

# Calculate % of samples that express each gene
genes.percent.expression <- rowMeans(counts_bulk>9 )*100  # Minimum of 10 counts to be assumed expressed with non-UMI RNAseq

# Keep only genes that are expressed in 30% of samples
genes.filter <- names(genes.percent.expression[genes.percent.expression>30])   # 13910 genes retained

# Filter counts
counts_bulk_filter <- counts_bulk[genes.filter,]

# TPM-normalize to take into account gene length  ------------------------------------------------------------

# Get exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")

# For each gene: reduce all exons to a set of non-overlapping exons, calculate lengths (widhts) and sum these
# Still only estimate, since no information of expression of which isoform per gene (or % of different isoforms)
exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))

# Filter for only expressed genes
exonic.gene.sizes.filter <- exonic.gene.sizes %>% dplyr::filter(rownames(exonic.gene.sizes) %in% rownames(counts_bulk_filter))

# Get exonic.gene.sizes in same order as counts
all(rownames(exonic.gene.sizes.filter) %in% rownames(counts_bulk_filter))
all(rownames(exonic.gene.sizes.filter) == rownames(counts_bulk_filter))

exonic.gene.sizes.filter %>% head()

# Calculate TPM matrix
TPM <- counts_bulk_filter / exonic.gene.sizes.filter$`sum(width(GenomicRanges::reduce(exons.list.per.gene)))`
tpm.mat <- t( t(TPM) * 1e6 / colSums(TPM) )

################################################################################################################################################
#################################################      ASSESS RIBO GENES      ######################################################################
################################################################################################################################################

# Single fibers -----------------------------------------------------------
rb.genes <- rownames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest)[grep("^RP[SL]",rownames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest))]
percent.ribo <- colSums(counts[rb.genes,])/Matrix::colSums(counts)*100
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- AddMetaData(filtered_normalized_fibertype_seurat_wo_MSTRG_rest, percent.ribo, col.name = "percent.ribo")

min(percent.ribo)
mean(percent.ribo)
max(percent.ribo)

pdf(file="11 Expression metrics/Fibers at rest//Ribo/VlnPlot ribo percentage.pdf")
VlnPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest, features = "percent.ribo", pt.size = 0.1) + NoLegend()
dev.off()

# Bulk muscle -----------------------------------------------------------
rb.genes_bulk <- rownames(tpm.mat)[grep("^RP[SL]",rownames(tpm.mat))]
percent.ribo_bulk <- colSums(tpm.mat[rb.genes_bulk,])/Matrix::colSums(tpm.mat)*100
percent.ribo_bulk <- as.data.frame(percent.ribo_bulk)

all(rownames(meta_bulk) %in% rownames(percent.ribo_bulk))
all(rownames(meta_bulk) == rownames(percent.ribo_bulk))
percent.ribo_bulk$subject <- meta_bulk$Subject
percent.ribo_bulk$subject <- factor(percent.ribo_bulk$subject)

min(percent.ribo_bulk$percent.ribo_bulk)
mean(percent.ribo_bulk$percent.ribo_bulk)
max(percent.ribo_bulk$percent.ribo_bulk)

plot_ribo_bulk <- ggplot(percent.ribo_bulk, aes(x=subject, y=percent.ribo_bulk)) +
  geom_violin(aes(fill=subject), trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  ylab("Relative expression (% total)") +
  xlab("Subject") +
  ggtitle("Relative expression ribosomal genes") +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(plot_ribo_bulk, filename = "11 Expression metrics/Fibers at rest//Ribo/VlnPlot ribo percentage bulk.pdf", width = 8, height = 5, scale=1.5)

################################################################################################################################################
#################################################      ASSESS MYH GENES      ######################################################################
################################################################################################################################################


# Single fibers -----------------------------------------------------------
MYH.genes <- rownames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest)[grep("^MY[HL]",rownames(filtered_normalized_fibertype_seurat_wo_MSTRG_rest))]
percent.MYH <- colSums(counts[MYH.genes,])/Matrix::colSums(counts)*100
filtered_normalized_fibertype_seurat_wo_MSTRG_rest <- AddMetaData(filtered_normalized_fibertype_seurat_wo_MSTRG_rest, percent.MYH, col.name = "percent.MYH")

min(percent.MYH)
mean(percent.MYH)
max(percent.MYH)

pdf(file="11 Expression metrics/Fibers at rest/MYH/VlnPlot MYH percentage single fibers.pdf")
VlnPlot(filtered_normalized_fibertype_seurat_wo_MSTRG_rest, features = "percent.MYH", pt.size = 0.1) + NoLegend()
dev.off()

# Bulk muscle -----------------------------------------------------------
MYH.genes_bulk <- rownames(tpm.mat)[grep("^MY[HL]",rownames(tpm.mat))]
percent.MYH_bulk <- colSums(tpm.mat[MYH.genes_bulk,])/Matrix::colSums(tpm.mat)*100
percent.MYH_bulk <- as.data.frame(percent.MYH_bulk)

all(rownames(meta_bulk) %in% rownames(percent.MYH_bulk))
all(rownames(meta_bulk) == rownames(percent.MYH_bulk))
percent.MYH_bulk$subject <- meta_bulk$Subject
percent.MYH_bulk$subject <- factor(percent.MYH_bulk$subject)

min(percent.MYH_bulk$percent.MYH_bulk)
mean(percent.MYH_bulk$percent.MYH_bulk)
max(percent.MYH_bulk$percent.MYH_bulk)

plot_MYH_bulk <- ggplot(percent.MYH_bulk, aes(x=subject, y=percent.MYH_bulk)) +
  geom_violin(aes(fill=subject), trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  ylab("Relative expression (% total)") +
  xlab("Subject") +
  ggtitle("Relative expression MYH genes") +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(plot_MYH_bulk, filename = "11 Expression metrics/Fibers at rest/MYH/VlnPlot MYH percentage bulk.pdf", width = 8, height = 5, scale=1.5)

################################################################################################################################################
#################################################      ASSESS HIGHEST EXPRESSED GENES      ######################################################################
################################################################################################################################################

# Single fibers -----------------------------------------------------------
counts_topexpressed <- counts
counts_topexpressed@x <- (counts_topexpressed@x/rep.int(colSums(counts_topexpressed), diff(counts_topexpressed@p))) * 100
most_expressed <- order(Matrix::rowSums(counts_topexpressed), decreasing = T)[20:1]

pdf(file="11 Expression metrics/Fibers at rest/Relative expression/Highest expressed genes.pdf")
boxplot(as.matrix(t(counts_topexpressed[most_expressed, ])),
        cex = 0.1,
        las = 1,
        xlab = "Percentage of total UMIs per fiber",
        main = "Top 20 expressed genes",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

# Bulk muscle -----------------------------------------------------------
counts_bulk_topexpressed <- tpm.mat
counts_bulk_topexpressed <- counts_bulk_topexpressed/colSums(counts_bulk_topexpressed) * 100
most_expressed_bulk <- order(Matrix::rowSums(counts_bulk_topexpressed), decreasing = T)[20:1]

pdf(file="11 Expression metrics/Fibers at rest/Relative expression/Highest expressed genes bulk.pdf")
boxplot(as.matrix(t(counts_bulk_topexpressed[most_expressed_bulk, ])),
        cex = 0.1,
        las = 1,
        xlab = "Percentage of total counts per fiber (TPM)",
        main = "Top 20 expressed genes (TPM)",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

################################################################################################################################################
#################################################     SATELLITE CELL MARKERS      ##############################################################
################################################################################################################################################


# Single fibers -----------------------------------------------------------
# PAX7
PAX7 <- "PAX7"
percent.PAX7 <- counts_nonfiltered[PAX7,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.PAX7, col.name = "percent.PAX7")

min(percent.PAX7)
mean(percent.PAX7)
max(percent.PAX7)

pdf(file="11 Expression metrics/Fibers at rest/Satellite cells/VlnPlot PAX7 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.PAX7", pt.size = 0.1) + NoLegend()
dev.off()

# NCAM1
NCAM1 <- "NCAM1"
percent.NCAM1 <- counts_nonfiltered[NCAM1,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.NCAM1, col.name = "percent.NCAM1")

min(percent.NCAM1)
mean(percent.NCAM1)
max(percent.NCAM1)

pdf(file="11 Expression metrics/Fibers at rest/Satellite cells/VlnPlot NCAM1 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.NCAM1", pt.size = 0.1) + NoLegend()
dev.off()

################################################################################################################################################
#################################################     MACROPHAGE MARKERS      ##############################################################
################################################################################################################################################

# MRC1 -----------------------------------------
MRC1 <- "MRC1"
percent.MRC1 <- counts_nonfiltered[MRC1,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.MRC1, col.name = "percent.MRC1")

min(percent.MRC1)
mean(percent.MRC1)
max(percent.MRC1)

pdf(file="11 Expression metrics/Fibers at rest/Macrophages/VlnPlot MRC1 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.MRC1", pt.size = 0.1) + NoLegend()
dev.off()

# C1QA -----------------------------------------
C1QA <- "C1QA"
percent.C1QA <- counts_nonfiltered[C1QA,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.C1QA, col.name = "percent.C1QA")

min(percent.C1QA)
mean(percent.C1QA)
max(percent.C1QA)

pdf(file="11 Expression metrics/Fibers at rest/Macrophages/VlnPlot C1QA percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.C1QA", pt.size = 0.1) + NoLegend()
dev.off()

################################################################################################################################################
#################################################     ENDOTHELIAL MARKERS      ##############################################################
################################################################################################################################################

# CDH5 -----------------------------------------
CDH5 <- "CDH5"
percent.CDH5 <- counts_nonfiltered[CDH5,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.CDH5, col.name = "percent.CDH5")

min(percent.CDH5)
mean(percent.CDH5)
max(percent.CDH5)

pdf(file="11 Expression metrics/Fibers at rest/Endothelial cells/VlnPlot CDH5 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.CDH5", pt.size = 0.1) + NoLegend()
dev.off()

# PECAM1 -----------------------------------------
PECAM1 <- "PECAM1"
percent.PECAM1 <- counts_nonfiltered[PECAM1,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.PECAM1, col.name = "percent.PECAM1")

min(percent.PECAM1)
mean(percent.PECAM1)
max(percent.PECAM1)

pdf(file="11 Expression metrics/Fibers at rest/Endothelial cells/VlnPlot PECAM1 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.PECAM1", pt.size = 0.1) + NoLegend()
dev.off()

################################################################################################################################################
#################################################     SMOOTH MUSCLE CELL MARKERS      ##############################################################
################################################################################################################################################

# ACTA2 -----------------------------------------
ACTA2 <- "ACTA2"
percent.ACTA2 <- counts_nonfiltered[ACTA2,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.ACTA2, col.name = "percent.ACTA2")

min(percent.ACTA2)
mean(percent.ACTA2)
max(percent.ACTA2)

pdf(file="11 Expression metrics/Fibers at rest/Smooth muscle cells/VlnPlot ACTA2 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.ACTA2", pt.size = 0.1) + NoLegend()
dev.off()

# MYH11 -----------------------------------------
MYH11 <- "MYH11"
percent.MYH11 <- counts_nonfiltered[MYH11,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.MYH11, col.name = "percent.MYH11")

min(percent.MYH11)
mean(percent.MYH11)
max(percent.MYH11)

pdf(file="11 Expression metrics/Fibers at rest/Smooth muscle cells/VlnPlot MYH11 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.MYH11", pt.size = 0.1) + NoLegend()
dev.off()

################################################################################################################################################
#################################################     FAP CELL MARKERS      ##############################################################
################################################################################################################################################

# PDGFRA -----------------------------------------
PDGFRA <- "PDGFRA"
percent.PDGFRA <- counts_nonfiltered[PDGFRA,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.PDGFRA, col.name = "percent.PDGFRA")

min(percent.PDGFRA)
mean(percent.PDGFRA)
max(percent.PDGFRA)

pdf(file="11 Expression metrics/Fibers at rest/FAP cells/VlnPlot PDGFRA percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.PDGFRA", pt.size = 0.1) + NoLegend()
dev.off()

# DCN -----------------------------------------
DCN <- "DCN"
percent.DCN <- counts_nonfiltered[DCN,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.DCN, col.name = "percent.DCN")

min(percent.DCN)
mean(percent.DCN)
max(percent.DCN)

pdf(file="11 Expression metrics/Fibers at rest/FAP cells/VlnPlot DCN percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.DCN", pt.size = 0.1) + NoLegend()
dev.off()

# CD34 -----------------------------------------
CD34 <- "CD34"
percent.CD34 <- counts_nonfiltered[CD34,]/Matrix::colSums(counts_nonfiltered)*100
seurat_nonfiltered <- AddMetaData(seurat_nonfiltered, percent.CD34, col.name = "percent.CD34")

min(percent.CD34)
mean(percent.CD34)
max(percent.CD34)

pdf(file="11 Expression metrics/Fibers at rest/FAP cells/VlnPlot CD34 percentage.pdf")
VlnPlot(seurat_nonfiltered, features = "percent.CD34", pt.size = 0.1) + NoLegend()
dev.off()

################################################################################################################################################
#################################################      DYNAMIC RANGE PLOT      #################################################################
################################################################################################################################################

# All genes ---------------------------------------------------------------
counts_nonfiltered_plot <- counts_nonfiltered

# Calculate relative expression for each gene for each fiber
counts_nonfiltered_plot@x <- (counts_nonfiltered_plot@x/rep.int(colSums(counts_nonfiltered_plot), diff(counts_nonfiltered_plot@p))) * 100

# Take average across fibers
mean_expression_all <- rowMeans(counts_nonfiltered_plot)
mean_expression_all <- as.data.frame(mean_expression_all)

# Rename mean expression column
mean_expression_all <- mean_expression_all %>% dplyr::rename(mean_expression = mean_expression_all)

# Add column with gene names
mean_expression_all$gene <- rownames(mean_expression_all)

# Create rank order
mean_expression_all <- mean_expression_all %>%
  dplyr::arrange(desc(mean_expression))
mean_expression_all$order <- seq_len(nrow(mean_expression_all))

# Filtered genes ---------------------------------------------------------------
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


# Create plot with all genes -------------------------------------------------------------

nonmuscle <- mean_expression_all %>%
  dplyr::filter(gene == "PAX7" | gene == "NCAM1" | gene ==  "MRC1" | gene ==  "C1QA" | gene ==  "CDH5" | gene ==  "PECAM1" | gene == "ACTA2" | gene ==  "MYH11" | gene ==  "PDGFRA" | gene ==  "DCN")
nonmuscle$type <- c("FAP", "SMC", "EC", "EC", "SC", "SMC", "MAC", "FAP", "SC", "MAC")


plot_expression_all <-  ggplot() +

  # Add all genes
  geom_point(data = mean_expression_all %>% filter(mean_expression != 0), aes(x=order, y=log(mean_expression,10)), colour="#B7DFB3", size=0.25) +

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

  # Add labels top 5 most expressed
  geom_label_repel(data = mean_expression_filtered %>% dplyr::filter(gene == "ACTA1" | gene == "MYH2" | gene == "MYH7" | gene == "TNNT1" | gene == "TNNT3"),
                   mapping = aes(order, log(mean_expression,10), label = gene),
                   size = 1.5, max.overlaps = Inf, label.padding=0.1, min.segment.length=0.1, segment.size=0.2, max.time = 10, fill = "white", force=80) +

  # Add labels non muscle cell markers
  geom_label_repel(nonmuscle,
                   mapping = aes(order+3, log(mean_expression,10), label = gene, fill=type),
                   size = 1.5, max.overlaps = Inf, label.padding=0.1, min.segment.length=0.1, segment.size=0.2,  force=80) +
  scale_fill_manual(values = c("#E8DFF5", "#FCE1E4", "#FCF4DD", "#9CADCE", "#DAEAF6")) +

  # Add dots non muscle cell markers

  # Add Satellite cells
  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PAX7"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "NCAM1"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

  # Add Macrophages
  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MRC1"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +


  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "C1QA"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +


  # Add Endothelial cells
  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "CDH5"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PECAM1"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

  # Add Smooth muscle cells
  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "ACTA2"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MYH11"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +


  # Add FAPs
  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PDGFRA"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

  geom_point(data = mean_expression_all %>% dplyr::filter(gene == "DCN"),
             aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +


  # Change design
  ylab("% total counts, 10log") +
  xlab("Gene rank") +
  theme_classic() +
  theme(
    text = element_text(face="bold", colour="black", size=6),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
  ) +
  scale_y_continuous(limits = c(-6.5,3), expand = c(0, 0))

ggsave(plot_expression_all, filename = "11 Expression metrics/Fibers at rest/Relative expression/Dynamic range Tx.png", width = 60, height =40, units="mm")

# Create plot with only one gene per non-muscle cell type -------------------------------------------------------------

# Satellite cell: PAX7
# Macrophage: MRC1
# Smooth muscle cell: MYH11
# Endothelial cell: CDH5
# FAP: PDGFRA
# NMJ: Chrne, Etv5, Colq, Musk and Lrp4

nonmuscle_simplified <- mean_expression_all %>%
    dplyr::filter(gene == "PAX7" | gene ==  "MRC1" | gene ==  "CDH5" | gene ==  "ACTA2" | gene ==  "PDGFRA" | gene ==  "CHRNE")
nonmuscle_simplified$type <- c("SMC", "EC", "FAP", "SC", "MAC", "NMJ")


plot_expression_all_simplified <-  ggplot() +

    # Add all genes
    geom_point(data = mean_expression_all %>% filter(mean_expression != 0), aes(x=order, y=log(mean_expression,10)), colour="#B7DFB3", size=0.25) +

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

    # Add labels non muscle cell markers
    geom_label_repel(nonmuscle_simplified,
                     mapping = aes(order, log(mean_expression,10), label = gene, fill=type),
                     size = 1.5, max.overlaps = Inf, label.padding=0.1, min.segment.length=0.1, segment.size=0.2,  force=20) +
    scale_fill_manual(values = c("#E8DFF5", "#FCE1E4", "#FCF4DD", "#C1E1C1", "#9CADCE", "#DAEAF6")) +

    # Add dots
    # Add Muscle cell markers
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "ACTA1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MYH2"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MYH7"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "TNNT1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "TNNT3"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Add Satellite cells
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PAX7"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Add Macrophages
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MRC1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Add Endothelial cells
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "CDH5"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Add Smooth muscle cells
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "ACTA2"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +


    # Add FAPs
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "PDGFRA"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Add NMJ
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "CHRNE"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +

    # Change design
    ylab("% total counts, 10log") +
    xlab("Gene rank") +
    theme_classic() +
    theme(
        text = element_text(face="bold", colour="black", size=6),
        strip.text = element_text(colour = "white"),
        strip.background = element_rect(fill="black"),
        legend.position = "none",
    ) +
    scale_y_continuous(limits = c(-6.5,3), expand = c(0, 0))

ggsave(plot_expression_all_simplified, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/Dynamic range Tx v2.png", width = 60, height = 40, units="mm")

# Create plot without non-muscle cell markers  -------------------------------------------------------------

plot_expression_all_without_nonmuscle <-  ggplot() +

    # Add all genes
    geom_point(data = mean_expression_all %>% filter(mean_expression != 0), aes(x=order, y=log(mean_expression,10)), colour="#B7DFB3", size=0.25) +

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
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "ACTA1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MYH2"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "MYH7"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "TNNT1"),
               aes(x=order, y=log(mean_expression,10)), colour="black", size=0.25) +
    geom_point(data = mean_expression_all %>% dplyr::filter(gene == "TNNT3"),
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

ggsave(plot_expression_all_without_nonmuscle, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/Dynamic_range_Tx_without_nonmuscle.png", width = 60, height = 60, units="mm")


      ###############################################################################################################################################
############################################      PLOT RELATIVE EXPRESSION FILTERED      ######################################################
###############################################################################################################################################

plot_expression_filter <-  ggplot() +

  # Add only filtered genes
  geom_point(data = mean_expression_filtered, aes(x=order, y=log(mean_expression,10)), colour="#134057") +

  # Top 10 most expressed genes
  geom_point(data = mean_expression_filtered %>% dplyr::slice(1:10), aes(x=order, y=log(mean_expression,10)), colour="#BC4749") +

  # Add labels top 10 most expressed
  geom_label_repel(data = mean_expression_filtered %>% dplyr::slice(1:10),
                   mapping = aes(order, log(mean_expression,10), label = gene),
                   size = 3, max.overlaps = Inf) +

  # Change design
  ylab("Relative expression (% total, 10log)") +
  xlab("Rank") +
  ggtitle("Relative expression filtered genes") +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=15, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(plot_expression_filter, filename = "11 Expression metrics/Fibers at rest/Relative expression/Relative expression filtered genes.pdf", width = 8, height = 5, scale=1.5)


################################################################################################################################################
#################################################      PLOT SINGLE FIBER VS BULK TRANSCRIPTOMICS      ##########################################
################################################################################################################################################

# All genes ---------------------------------------------------------------
counts_nonfiltered_plot <- counts_nonfiltered

# Calculate relative expression for each gene for each fiber
counts_nonfiltered_plot@x <- (counts_nonfiltered_plot@x/rep.int(colSums(counts_nonfiltered_plot), diff(counts_nonfiltered_plot@p))) * 100

# Take average across fibers
mean_expression_all <- rowMeans(counts_nonfiltered_plot)
mean_expression_all <- as.data.frame(mean_expression_all)

# Rename mean expression column
mean_expression_all <- mean_expression_all %>% dplyr::rename(mean_expression = mean_expression_all)

# Add column with gene names
mean_expression_all$gene <- rownames(mean_expression_all)

# Create rank order
mean_expression_all <- mean_expression_all %>%
  dplyr::arrange(desc(mean_expression))
mean_expression_all$order <- seq_len(nrow(mean_expression_all))

