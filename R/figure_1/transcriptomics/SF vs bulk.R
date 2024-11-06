################################################################################################################################################
################################################       PREPARATION      ####@###################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(scater)
library(rtracklayer)
library(GenomicFeatures)
library(ggrepel)
library(VennDiagram)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(eulerr)

BiocManager::install("GenomicFeatures")


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("7 Targeted fiber type markers/Fibers at rest/Inflection/filtered_normalized_fibertype_seurat_wo_MSTRG_rest.RData")

# Load bulk data
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Bulk muscle RNAseq/Own data")

counts_bulk <- read.table("2 Data preparation/counts.txt", header = TRUE, sep = "\t", dec = ".")
meta_bulk <- read.table("2 Data preparation/meta.txt", header = TRUE, sep = "\t", dec = ".")

txdb <- makeTxDbFromGFF("2 Data preparation/TotalRNA_polyA_transcriptome_filtered.gtf",format="gtf")

setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

################################################################################################################################################
#################################################      EXTRACT FILTERED COUNTS      ############################################################
################################################################################################################################################

singlefiber <- GetAssayData(object = filtered_normalized_fibertype_seurat_wo_MSTRG_rest, assay = "RNA", slot = "counts")

################################################################################################################################################
#################################################      CREATE DF WITH SF AND BULK DATA      ####################################################
################################################################################################################################################

# Single fiber
singlefiber_raw <- singlefiber

mean_counts_singlefiber <- rowMeans(singlefiber_raw)
mean_counts_singlefiber <- as.data.frame(mean_counts_singlefiber)
mean_counts_singlefiber$gene <- rownames(mean_counts_singlefiber)


mean_counts_singlefiber <- mean_counts_singlefiber %>%
  dplyr::arrange(desc(mean_counts_singlefiber))
mean_counts_singlefiber$order_T <- seq_len(nrow(mean_counts_singlefiber))

# Bulk

# Remove MSTRG artefact genes from counts ------------------------------------------------------------
counts_bulk <- counts_bulk[!grepl('MSTRG',rownames(counts_bulk)),]

# Only use samples at rest
meta_bulk$time <- gsub(".*_", "", meta_bulk$Group)
pre_samples_bulk <- meta_bulk %>% dplyr::filter(time == "Pre") %>% rownames()
counts_bulk_pre <- counts_bulk[,pre_samples_bulk]

# Filter genes  ------------------------------------------------------------

# Calculate % of samples that express each gene
genes.percent.expression <- rowMeans(counts_bulk_pre>9 )*100  # Minimum of 10 counts to be assumed expressed with non-UMI RNAseq

# Keep only genes that are expressed in 30% of samples
genes.filter <- names(genes.percent.expression[genes.percent.expression>30])   # 13820 genes retained

# Filter counts
counts_bulk_pre_filter <- counts_bulk_pre[genes.filter,]

# TPM-normalize to take into account gene length  ------------------------------------------------------------

# Get exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")

# For each gene: reduce all exons to a set of non-overlapping exons, calculate lengths (widhts) and sum these
# Still only estimate, since no information of expression of which isoform per gene (or % of different isoforms)
exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))

# Filter for only expressed genes
exonic.gene.sizes.filter <- exonic.gene.sizes %>% dplyr::filter(rownames(exonic.gene.sizes) %in% rownames(counts_bulk_pre_filter))

# Get exonic.gene.sizes in same order as counts
all(rownames(exonic.gene.sizes.filter) %in% rownames(counts_bulk_pre_filter))
all(rownames(exonic.gene.sizes.filter) == rownames(counts_bulk_pre_filter))

# Calculate TPM matrix
TPM <- counts_bulk_pre_filter / exonic.gene.sizes.filter$`sum(width(GenomicRanges::reduce(exons.list.per.gene)))`
tpm.mat <- t( t(TPM) * 1e6 / colSums(TPM) )

mean_counts_bulk <- rowMeans(tpm.mat, na.rm=TRUE)
mean_counts_bulk <- as.data.frame(mean_counts_bulk)
mean_counts_bulk <- na.omit(mean_counts_bulk)
mean_counts_bulk$gene <- rownames(mean_counts_bulk)

mean_counts_bulk <- mean_counts_bulk %>%
  dplyr::arrange(desc(mean_counts_bulk))
mean_counts_bulk$order_bulk <- seq_len(nrow(mean_counts_bulk))

# Filter singlefiber and bulk to only keep genes detected in both datasets
mean_counts_singlefiber_corr <- subset(mean_counts_singlefiber, gene %in% mean_counts_bulk$gene)
mean_counts_bulk_corr <- subset(mean_counts_bulk, gene %in% mean_counts_singlefiber$gene)

# Combine dataframes
combined_raw <- merge(mean_counts_singlefiber_corr, mean_counts_bulk_corr, by.x="gene", by.y="gene")


################################################################################################################################################
#############################################     Venn diagram detected genes singlefiber vs bulk      #########################################
################################################################################################################################################

Venn_detected_genes <- venn.diagram(
  # General
  filename=NULL,
  disable.logging=T,
  x = list(
    mean_counts_singlefiber %>% dplyr::select(gene) %>% unlist(use.names=F),
    rownames(counts_bulk_pre_filter) %>% unlist(use.names=F)
  ),
  category.names = c("Single fiber" , "Bulk"),
  main.fontface = "bold",
  main.fontfamily = "sans",
  main.cex = 0.5,

  # Circles
  lwd = 2,
  col=c("#4F7A5D", "#B7DFB3"),
  fill = c(alpha("#4F7A5D",0.3), alpha('#B7DFB3',0.3)),

  # Numbers
  cex = 0.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.distance = c(0.05, 0.02),

  # Names
  cat.cex = 0.75,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25,25),
  cat.dist = c(-0.09, -0.05),
  cat.col = c("#4F7A5D", "#B7DFB3")
)

ggsave(Venn_detected_genes, filename = "11 Expression metrics/Fibers at rest/Single fiber vs bulk/Detected sf vs bulk.png", width = 40, height = 40, units="mm")



################################################################################################################################################
###############################################     CREATE PLOT WITH RAW COUNTS CORRELATION      ###############################################
################################################################################################################################################

# Calculate Spearman correlation
cor.test(combined_raw$mean_counts_singlefiber, combined_raw$mean_counts_bulk,  method = "spearman") # p = 0.54 (non-log)
cor.test(combined_raw$mean_counts_singlefiber, combined_raw$mean_counts_bulk,  method = "pearson") # p = 0.91 (non-log)
cor.test(log(combined_raw$mean_counts_singlefiber,10), log(combined_raw$mean_counts_bulk,10),  method = "pearson") # p = 0.63 (log)

# Create plot

SFvsBulk_plot <- combined_raw %>%
  ggplot() +
  geom_point(aes(x=log(mean_counts_singlefiber,10), y=log(mean_counts_bulk,10)), size=0.1) +
  geom_smooth(aes(x=log(mean_counts_singlefiber,10), y=log(mean_counts_bulk,10)), method='lm', formula= y~x, se=F, colour = "#28666E", size=0.75) +
  xlab("Single fiber (avg counts, log10)") +
  ylab("Bukl (avg intensity, log10)") +
  annotate("text", x=3, y=0, label= "r = 0.63", colour="black", fontface=2, size=2) +
  theme_classic() +
  theme(
    text = element_text(face="bold", size=6, colour="black"),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none"
  )

ggsave(SFvsBulk_plot, filename = "~/single_fiber_heterogeneity/doc/figures/figure_1/Correlation sf vs bulk.png", width = 60, height = 60, units="mm")


################################################################################################################################################
#####################################################     NON-CODING RNA SINGLE FIBER       ####################################################
################################################################################################################################################

# Load in gene annotation file --------------------------------------------
annotations_sf <- read.csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

genes_sf <- rownames(mean_counts_singlefiber)

genes_sf <- data.frame(
  gene = genes_sf
)

# Add gene annotations
genes_sf <- left_join(genes_sf, annotations_sf, by = c("gene" = "GENENAME"))




# Visual RNA biotype distribution -----------------------------------------------

# Determine how many genes per type
genes_sf_type <- genes_sf %>%
  group_by(GENEBIOTYPE) %>%
  summarise(n = n())

# Filter
genes_sf_type <- genes_sf_type %>% dplyr::filter(GENEBIOTYPE == "antisense" | GENEBIOTYPE == "lincRNA" | GENEBIOTYPE == "protein_coding")

# Compute percentages
genes_sf_type$fraction = genes_sf_type$n / sum(genes_sf_type$n)

# Compute the cumulative percentages (top of each rectangle)
genes_sf_type$ymax = cumsum(genes_sf_type$fraction)

# Compute the bottom of each rectangle
genes_sf_type$ymin = c(0, head(genes_sf_type$ymax, n=-1))

# Compute label position
genes_sf_type$labelPosition.y <- (genes_sf_type$ymax + genes_sf_type$ymin) / 2

genes_sf_type$labelPosition.x <- c(1,2.2,2)

# Compute a good label
genes_sf_type$label <- paste0(genes_sf_type$GENEBIOTYPE, "\n n: ", genes_sf_type$n)

# Make the plot
plot_biotype <- ggplot(genes_sf_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=GENEBIOTYPE)) +
  geom_rect() +
  geom_label(aes(x= labelPosition.x, y=labelPosition.y, label=label), size=2, label.padding= unit(0.1, "lines")) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-0.1, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  theme(
    text = element_text(face="bold", size=8, colour="black"),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(colour = "white"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(plot_biotype, file="11 Expression metrics/Fibers at rest/Single fiber vs bulk/Gene biotype sf.png", width = 60, height = 60, units="mm")


################################################################################################################################################
#####################################################     NON-CODING RNA BULK       ############################################################
################################################################################################################################################

# Load in gene annotation file --------------------------------------------

genes_bulk <- rownames(counts_bulk_pre_filter)

genes_bulk <- data.frame(
  gene = genes_bulk
)

# Add gene annotations
genes_bulk <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                  keys = genes_bulk$gene,
                                  columns = c("GENENAME","GENEID","GENEBIOTYPE"),
                                  keytype = "GENENAME")

genes_bulk <- genes_bulk[!duplicated(genes_bulk$GENENAME),]

# Check if correct amount of genes
all(counts_bulk_pre_filter$gene == genes_bulk$GENENAME)

# Visual RNA biotype distribution -----------------------------------------------

# Determine how many genes per type
genes_bulk_type <- genes_bulk %>%
  group_by(GENEBIOTYPE) %>%
  summarise(n = n())

# Filter
genes_bulk_type <- genes_bulk_type %>% dplyr::filter(GENEBIOTYPE == "antisense" | GENEBIOTYPE == "lincRNA" | GENEBIOTYPE == "protein_coding")

# Compute percentages
genes_bulk_type$fraction = genes_bulk_type$n / sum(genes_bulk_type$n)

# Compute the cumulative percentages (top of each rectangle)
genes_bulk_type$ymax = cumsum(genes_bulk_type$fraction)

# Compute the bottom of each rectangle
genes_bulk_type$ymin = c(0, head(genes_bulk_type$ymax, n=-1))

# Compute label position
genes_bulk_type$labelPosition.y <- (genes_bulk_type$ymax + genes_bulk_type$ymin) / 2

genes_bulk_type$labelPosition.x <- c(1,2.2,2)

# Compute a good label
genes_bulk_type$label <- paste0(genes_bulk_type$GENEBIOTYPE, "\n n: ", genes_bulk_type$n)

# Make the plot
plot_biotype_bulk <- ggplot(genes_bulk_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=GENEBIOTYPE)) +
  geom_rect() +
  geom_label(aes(x= labelPosition.x, y=labelPosition.y, label=label), size=2, label.padding= unit(0.1, "lines")) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-0.1, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  theme(
    text = element_text(face="bold", size=8, colour="black"),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(plot_biotype_bulk, file="11 Expression metrics/Fibers at rest/Single fiber vs bulk/Gene biotype bulk.png", width = 60, height = 60, units="mm")


################################################################################################################################################
#####################################################     NON-CODING RNA SINGLE FIBER VS BULK       ############################################
################################################################################################################################################

# Combine both data frames
combined <- data.frame(
  Experiment = rep(c("bulk", "sf"), each=3),
  Biotype = rep(c("Antisense", "LincRNA", "Protein coding"),2),
  Fraction = c(genes_bulk_type$fraction*100,genes_sf_type$fraction*100)
)

# Get label coordinates
combined_label <- combined %>% dplyr::filter(Experiment == "bulk")
combined_label$labelPosition.y <- c(1, 92, 10)
combined_label$labelPosition.x <- c(4.3,4,4.3)

# Create plot
Donut <- ggplot(combined, aes(x = Experiment, y = Fraction, fill = Biotype)) +
  geom_col(aes(fill=Biotype)) +
  geom_label(data = combined_label, aes(x= labelPosition.x, y=labelPosition.y, label=Biotype), size=2, label.padding= unit(0.1, "lines")) +
  annotate("text",x = 2, y = 50,label = "Single \nfiber",colour = "black",fontface = 2,size=1.5) +
  annotate("text",x = 3, y = 50,label = "Bulk",colour = "black",fontface = 2,size=1.5) +
  scale_x_discrete(limits = c(" ", "bulk","sf")) +
  scale_fill_brewer(palette=2) +
  coord_polar("y") +
  theme_void() +
  theme(
    text = element_text(face="bold", size=2, colour="black"),
    legend.position="none"
  )

ggsave(Donut, file="11 Expression metrics/Fibers at rest/Single fiber vs bulk/Biotype sf vs bulk.png", width = 40, height = 60, units="mm")



