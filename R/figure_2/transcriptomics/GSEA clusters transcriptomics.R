################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

# Packages ----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(viridis)
library(RCurl)
library(cowplot)
library(ggpubr)
library(VennDiagram)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(org.Hs.eg.db)


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

# Assign identity of final clusters ---------------------------------------------
Idents(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"

# Load in gene annotation file --------------------------------------------
annotations <- read.csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

# Set to SCT assay for visualization ----------------------------
DefaultAssay(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "SCT"

# Load file with all markers, via FindAllMarkers, for each cluster
markers_all <- read.csv("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/Markers/markers_all.csv")

# All genes in Seurat object as background for ORA analysis:
vargenes_SCT <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@assays$SCT@var.features

# Average expression per gene per cluster, used for heatmaps to compare clusters---------------------------------------------------------------------
seurat_averaged <- AverageExpression(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, 
                                     assays = "SCT",
                                     return.seurat = T)

################################################################################################################################################
##############################################      CHECK CORRECT CLUSTERING     ###############################################################
################################################################################################################################################

# Check TSNE to confirm correct cluster assignment
DimPlot(object = filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, 
        reduction = "tsne", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)



################################################################################################################################################
#####################################################      FILTER MARKER GENES     #############################################################
################################################################################################################################################


# Filter genes only with at least absolute value of 1 L2FC------------------------------------------------------------------
markers_all_filtered <- markers_all %>% dplyr::filter(avg_log2FC >= 1 | avg_log2FC <= -1)

################################################################################################################################################
#############################################      GET TOP 50 MARKERS FOR EACH CLUSTER     #####################################################
################################################################################################################################################

# Split per cluster
markers_Fast1 <- markers_all_filtered %>% dplyr::filter(cluster == "Fast1") 
markers_Fast2 <- markers_all_filtered %>% dplyr::filter(cluster == "Fast2") 
markers_Slow1 <- markers_all_filtered %>% dplyr::filter(cluster == "Slow1")
markers_Slow2 <- markers_all_filtered %>% dplyr::filter(cluster == "Slow2")


################################################################################################################################################
##############################################     VENN PLOT OVERLAP CLUSTERS - L2FC    ########################################################
################################################################################################################################################

Venn_clusters <- venn.diagram(
  # General
  filename=NULL,
  disable.logging=T,
  x = list(
    markers_Fast1 %>% as.data.frame() %>% dplyr::select(gene) %>% unlist(use.names=F), 
    markers_Fast2 %>% as.data.frame() %>% dplyr::select(gene) %>% unlist(use.names=F),
    markers_Slow1 %>% as.data.frame() %>% dplyr::select(gene) %>% unlist(use.names=F),
    markers_Slow2 %>% as.data.frame() %>% dplyr::select(gene) %>% unlist(use.names=F)
  ),
  category.names = c("Fast1" , "Fast2", "Slow1", "Slow2"),
  main = "Markers per cluster (L2FC > 1)",
  main.fontface = "bold",
  main.fontfamily = "sans",
  main.cex = 2,
  
  
  # Circles
  lwd = 2,
  col=c("#134057", "#8CB3E8", "#BC4749", "#417B5A"),
  fill = c(alpha("#134057",0.3), alpha('#8CB3E8',0.3), alpha("#BC4749",0.3), alpha("#417B5A",0.3)),
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Names
  cat.cex = 2,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25,25, 0, 0),
  cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.col = c("#134057", "#8CB3E8", "#BC4749", "#417B5A")
)

grid.newpage()
grid.draw(Venn_clusters)

pdf(file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Venn/Venn clusters.pdf")
grid.draw(Venn_clusters)
dev.off()

################################################################################################################################################
#######################################################     SLOW1 VS SLOW2    ##################################################################
################################################################################################################################################

# Compare markers Slow1 vs Slow2 ----------------------------

# Find markers
markers_Slow1vsSlow2 <- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    ident.1 = "Slow1",
                                    ident.2 = "Slow2")    

markers_Slow1vsSlow2 <- markers_Slow1vsSlow2 %>% 
  rownames_to_column(var = "gene")

# Add gene annotations 
markers_Slow1vsSlow2 <- left_join(markers_Slow1vsSlow2, annotations, by = c("gene" = "GENENAME"))


# Make column with pct diff
markers_Slow1vsSlow2$pct.diff <- markers_Slow1vsSlow2$pct.1 - markers_Slow1vsSlow2$pct.2

# Save markers
write_csv(markers_Slow1vsSlow2, file = "8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Markers/Slow1 vs Slow2/Slow1 vs Slow2.csv")

# Filter genes only with at least absolute value of 1 L2FC------------------------------------------------------------------
markers_Slow1vsSlow2_filtered_all <- markers_Slow1vsSlow2 %>% dplyr::filter(avg_log2FC >= 1 | avg_log2FC <= -1) 

markers_Slow1vsSlow2_filtered_top50 <- markers_Slow1vsSlow2_filtered_all %>% top_n(n = 50, wt = abs(avg_log2FC)) %>% arrange(avg_log2FC)

# Subset Slow1 and Slow2 clusters ---------------------------------------------------------------------
Seurat_averaged_Slow1_Slow2 <- subset(x = seurat_averaged, idents = c("Slow1", "Slow2"))


# Create heatmap ----------------------------------------------------------
pdf(file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Markers/Slow1 vs Slow2/heatmap top50.pdf", width=15, height=9)
Seurat::DoHeatmap(Seurat_averaged_Slow1_Slow2, features = markers_Slow1vsSlow2_filtered_top50$gene, raster=F) + 
  scale_fill_gradientn(colors = c("#134057", "white", "#BC4749"))
dev.off()


# Perform GSEA with clusterProfiler for Slow 2 --------------------------------------------------------------------

# Marker genes Slow 2:
markers_Slow2 <- markers_Slow1vsSlow2 %>% dplyr::filter(avg_log2FC <= -1) 

# Perform ORA with background only top 3000 genes SCT -------------------------------------------------------------
ORA_Slow2 <- enrichGO(gene = markers_Slow2$gene,
                        universe      = vargenes_SCT,
                        keyType       = "SYMBOL",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)


# Simplify GO terms redundancy --------------------------------------------
ORA_Slow2_simplify <- simplify(ORA_Slow2, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
ORA_Slow2_simplify_res <- as.data.frame(ORA_Slow2_simplify)

# Rich factor: ratio of PC genes that are annotated in term to all genes annotated in this term --------------------------------------------
ORA_Slow2_simplify_res <- mutate(ORA_Slow2_simplify_res, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
ORA_Slow2_simplify_res <- mutate(ORA_Slow2_simplify_res, foldEnrich = 
                                     (as.numeric(sub("/\\d+", "", ORA_Slow2_simplify_res$GeneRatio)) / as.numeric(sub(".*/", "", ORA_Slow2_simplify_res$GeneRatio))) /
                                     (as.numeric(sub("/\\d+", "", ORA_Slow2_simplify_res$BgRatio)) / as.numeric(sub(".*/", "", ORA_Slow2_simplify_res$BgRatio)))
)


# Slice first 15 GO terms based on Count-------------------------------------------------
ORA_Slow2_simplify_res <- ORA_Slow2_simplify_res %>% arrange(desc(Count)) %>% dplyr::slice(1:15)

# Visualize --------------------------------------------
plot_Slow2vsSlow1 <- ggplot(ORA_Slow2_simplify_res,aes(foldEnrich,
                                                        fct_reorder(Description, foldEnrich))) + 
  geom_segment(aes(xend=0, yend = Description)) + 
  geom_point(aes(color=p.adjust, size = richFactor)) + 
  scale_color_gradientn(
    colours=c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10", 
    guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) + 
  theme_classic() +
  xlab("Fold Enrichment") +
  ggtitle("Slow2 vs Slow1") +
  theme(
    axis.title.y= element_blank(),
    plot.title = element_text(hjust = 0.5, face="bold")
  ) 

# Save results ------------------------------------------------------------
ggsave(plot_Slow2vsSlow1, file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/ORA/Slow2 vs Slow1/Slow2 vs Slow1.pdf", width=10, height=4)
write_csv(ORA_Slow2_simplify_res, file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/ORA/Slow2 vs Slow1/Slow2 vs Slow1.csv")


# Visual for biotype of RNA -----------------------------------------------

# Determine how many genes per type
markers_Slow1vsSlow2_type <- markers_Slow1vsSlow2_filtered_all %>% 
  group_by(GENEBIOTYPE) %>% 
  summarise(n = n())

# Compute percentages
markers_Slow1vsSlow2_type$fraction = markers_Slow1vsSlow2_type$n / sum(markers_Slow1vsSlow2_type$n)

# Compute the cumulative percentages (top of each rectangle)
markers_Slow1vsSlow2_type$ymax = cumsum(markers_Slow1vsSlow2_type$fraction)

# Compute the bottom of each rectangle
markers_Slow1vsSlow2_type$ymin = c(0, head(markers_Slow1vsSlow2_type$ymax, n=-1))

# Compute label position
markers_Slow1vsSlow2_type$labelPosition <- (markers_Slow1vsSlow2_type$ymax + markers_Slow1vsSlow2_type$ymin) / 2

# Compute a good label
markers_Slow1vsSlow2_type$label <- paste0(markers_Slow1vsSlow2_type$GENEBIOTYPE, "\n n: ", markers_Slow1vsSlow2_type$n)

# Make the plot
plot_slow2_type <- ggplot(markers_Slow1vsSlow2_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=GENEBIOTYPE)) +
  geom_rect() +
  geom_label( x=2, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-0.1, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("RNA Biotype: Slow2 markers (L2FC > 1)") + 
theme(
  text = element_text(face="bold", size=25, colour="black"),
  axis.text = element_blank(),
  axis.title.x = element_blank(),
  strip.text = element_text(colour = "white"),
  strip.background = element_rect(fill="black"),
  legend.position = "none",
  plot.title = element_text(hjust = 0.5)
)

# Save results ------------------------------------------------------------
ggsave(plot_slow2_type, file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Biotype/Slow2 markers biotype.pdf")

################################################################################################################################################
#######################################################     Fast1 VS Fast2    ##################################################################
################################################################################################################################################

# Compare markers Fast1 vs Fast2 ----------------------------

# Find markers
markers_Fast1vsFast2 <- FindMarkers(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest,
                                    ident.1 = "Fast1",
                                    ident.2 = "Fast2")    

markers_Fast1vsFast2 <- markers_Fast1vsFast2 %>% 
  rownames_to_column(var = "gene")

# Add gene annotations 
markers_Fast1vsFast2 <- left_join(markers_Fast1vsFast2, annotations, by = c("gene" = "GENENAME"))


# Make column with pct diff
markers_Fast1vsFast2$pct.diff <- markers_Fast1vsFast2$pct.1 - markers_Fast1vsFast2$pct.2

# Save markers
write_csv(markers_Fast1vsFast2, file = "8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Markers/Fast1 vs Fast2/Fast1 vs Fast2.csv")

# Filter genes only with at least absolute value of 1 L2FC------------------------------------------------------------------
markers_Fast1vsFast2_filtered_all <- markers_Fast1vsFast2 %>% dplyr::filter(avg_log2FC >= 1 | avg_log2FC <= -1) 

markers_Fast1vsFast2_filtered_top50 <- markers_Fast1vsFast2_filtered_all %>% top_n(n = 50, wt = abs(avg_log2FC)) %>% arrange(avg_log2FC)

# Subset Fast1 and Fast2 clusters ---------------------------------------------------------------------
Seurat_averaged_Fast1_Fast2 <- subset(x = seurat_averaged, idents = c("Fast1", "Fast2"))


# Create heatmap ----------------------------------------------------------
pdf(file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Markers/Fast1 vs Fast2/heatmap top50.pdf", width=15, height=9)
Seurat::DoHeatmap(Seurat_averaged_Fast1_Fast2, features = markers_Fast1vsFast2_filtered_top50$gene, raster=F) + 
  scale_fill_gradientn(colors = c("#134057", "white", "#BC4749"))
dev.off()


# Perform GSEA with clusterProfiler for Fast 2 --------------------------------------------------------------------

# Marker genes Fast2 (only top 200 genes):
markers_Fast2 <- markers_Fast1vsFast2 %>% dplyr::filter(avg_log2FC <= -1) %>% top_n(n = 200, wt = abs(avg_log2FC))

# Perform ORA with background only top 3000 genes SCT -------------------------------------------------------------
ORA_Fast2 <- enrichGO(gene = markers_Fast2$gene,
                      universe      = vargenes_SCT,
                      keyType       = "SYMBOL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)


# Simplify GO terms redundancy --------------------------------------------
ORA_Fast2_simplify <- simplify(ORA_Fast2, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
ORA_Fast2_simplify_res <- as.data.frame(ORA_Fast2_simplify)

# Rich factor: ratio of PC genes that are annotated in term to all genes annotated in this term --------------------------------------------
ORA_Fast2_simplify_res <- mutate(ORA_Fast2_simplify_res, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
ORA_Fast2_simplify_res <- mutate(ORA_Fast2_simplify_res, foldEnrich = 
                                   (as.numeric(sub("/\\d+", "", ORA_Fast2_simplify_res$GeneRatio)) / as.numeric(sub(".*/", "", ORA_Fast2_simplify_res$GeneRatio))) /
                                   (as.numeric(sub("/\\d+", "", ORA_Fast2_simplify_res$BgRatio)) / as.numeric(sub(".*/", "", ORA_Fast2_simplify_res$BgRatio)))
)


# Slice first 15 GO terms based on Count-------------------------------------------------
ORA_Fast2_simplify_res <- ORA_Fast2_simplify_res %>% arrange(p.adjust) %>% dplyr::slice(1:15)

# Visualize --------------------------------------------
plot_Fast2vsFast1 <- ggplot(ORA_Fast2_simplify_res,aes(Count,
                                                       fct_reorder(Description, Count))) + 
  geom_segment(aes(xend=0, yend = Description)) + 
  geom_point(aes(color=p.adjust, size = foldEnrich)) + 
  scale_color_gradientn(
    colours=c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10", 
    guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) + 
  theme_classic() +
  xlab("GeneRatio") +
  ggtitle("Fast2 vs Fast1") +
  theme(
    axis.title.y= element_blank(),
    plot.title = element_text(hjust = 0.5, face="bold")
  ) 


# Save results ------------------------------------------------------------
ggsave(plot_Fast2vsFast1, file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/ORA/Fast2 vs Fast1/Fast2 vs Fast1.pdf", width=10, height=5)
write_csv(ORA_Fast2_simplify_res, file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/ORA/Fast2 vs Fast1/Fast2 vs Fast1.csv")

# Visual for biotype of RNA -----------------------------------------------

# Determine how many genes per type
markers_Fast1vsFast2_type <- markers_Fast1vsFast2_filtered_all %>% 
  group_by(GENEBIOTYPE) %>% 
  summarise(n = n()) %>% 
  dplyr::filter(GENEBIOTYPE == "antisense" | GENEBIOTYPE == "lincRNA"| GENEBIOTYPE == "protein_coding")

# Compute percentages
markers_Fast1vsFast2_type$fraction = markers_Fast1vsFast2_type$n / sum(markers_Fast1vsFast2_type$n)

# Compute the cumulative percentages (top of each rectangle)
markers_Fast1vsFast2_type$ymax = cumsum(markers_Fast1vsFast2_type$fraction)

# Compute the bottom of each rectangle
markers_Fast1vsFast2_type$ymin = c(0, head(markers_Fast1vsFast2_type$ymax, n=-1))

# Compute label position
markers_Fast1vsFast2_type$labelPosition <- (markers_Fast1vsFast2_type$ymax + markers_Fast1vsFast2_type$ymin) / 2

# Compute a good label
markers_Fast1vsFast2_type$label <- paste0(markers_Fast1vsFast2_type$GENEBIOTYPE, "\n n: ", markers_Fast1vsFast2_type$n)

# Make the plot
plot_Fast2_type <- ggplot(markers_Fast1vsFast2_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=GENEBIOTYPE)) +
  geom_rect() +
  geom_label( x=2, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-0.1, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("RNA Biotype: Fast2 markers") + 
  theme(
    text = element_text(face="bold", size=25, colour="black"),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill="black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Save results ------------------------------------------------------------
ggsave(plot_Fast2_type, file="8 Fiber heterogeneity (only rested samples)/4 GSEA Clusters/Biotype/Fast2 markers biotype.pdf")



