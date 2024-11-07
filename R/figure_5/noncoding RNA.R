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
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(eulerr)
library(biomaRt)

BiocManager::install("biomaRt")


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load annotations
genes <- read_csv("12 Annotation table creation/Fibers at rest/Annotation_rest.csv")

################################################################################################################################################
#################################################      EXTRACT LIST WITH NON-CODING RNA      ###################################################
################################################################################################################################################

# Filter only non-coding
nc_genes <- genes %>% dplyr::filter(GENEBIOTYPE != "protein_coding")
write_csv(nc_genes, file = "~/single_fiber_heterogeneity/data/noncoding_RNA_list_transcriptomics/single_fiber_noncoding_genes_transcriptomics.csv")

################################################################################################################################################
#######################################################      GET FULL SEQUENCES      ###########################################################
################################################################################################################################################

# Acquire ensmebl biomaRt object
ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")

# Get sequences for gene 1 to 30
nc_genes_1to30 <- nc_genes %>% dplyr::slice(1:30) %>% dplyr::filter(GENEID != "ENSG00000228651")
files_1to30 <- nc_genes_1to30$GENEID

sequence_list_noncoding_1to30 <- list()

for(i in files_1to30){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_1to30[[symbol]] <- tmp
}

# Get sequences for gene 31 to 60
nc_genes_31to60 <- nc_genes %>% dplyr::slice(31:60) %>% dplyr::filter(GENEID != "ENSG00000255864")  %>% dplyr::filter(GENEID != "ENSG00000237152")
files_31to60 <- nc_genes_31to60$GENEID

sequence_list_noncoding_31to60 <- list()

for(i in files_31to60){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_31to60[[symbol]] <- tmp
}

# Get sequences for gene 61 to 70
nc_genes_61to70 <- nc_genes %>% dplyr::slice(61:70) %>% dplyr::filter(GENEID != "ENSG00000257151")
files_61to70 <- nc_genes_61to70$GENEID

sequence_list_noncoding_61to70 <- list()

for(i in files_61to70){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_61to70[[symbol]] <- tmp
}

# Get sequences for gene 71 to 90
nc_genes_71to90 <- nc_genes %>% dplyr::slice(71:90)
files_71to90 <- nc_genes_71to90$GENEID

sequence_list_noncoding_71to90 <- list()

for(i in files_71to90){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_71to90[[symbol]] <- tmp
}

# Get sequences for gene 91 to 110
nc_genes_91to110 <- nc_genes %>% dplyr::slice(91:110)  %>% dplyr::filter(GENEID != "ENSG00000267194")
files_91to110 <- nc_genes_91to110$GENEID

sequence_list_noncoding_91to110 <- list()

for(i in files_91to110){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_91to110[[symbol]] <- tmp
}

# Get sequences for gene 111 to 130
nc_genes_111to130 <- nc_genes %>% dplyr::slice(111:130)
files_111to130 <- nc_genes_111to130$GENEID

sequence_list_noncoding_111to130 <- list()

for(i in files_111to130){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_111to130[[symbol]] <- tmp
}

# Get sequences for gene 131 to 150
nc_genes_131to150 <- nc_genes %>% dplyr::slice(131:150) %>% dplyr::filter(GENEID != "ENSG00000163009") %>% dplyr::filter(GENEID != "ENSG00000203386")
files_131to150 <- nc_genes_131to150$GENEID

sequence_list_noncoding_131to150 <- list()

for(i in files_131to150){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_131to150[[symbol]] <- tmp
}

# Get sequences for gene 151 to 160
nc_genes_151to160 <- nc_genes %>% dplyr::slice(151:160)
files_151to160 <- nc_genes_151to160$GENEID

sequence_list_noncoding_151to160 <- list()

for(i in files_151to160){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_151to160[[symbol]] <- tmp
}

# Get sequences for gene 161 to 170
nc_genes_161to170 <- nc_genes %>% dplyr::slice(161:170)
files_161to170 <- nc_genes_161to170$GENEID

sequence_list_noncoding_161to170 <- list()

for(i in files_161to170){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_161to170[[symbol]] <- tmp
}

# Get sequences for gene 171 to 180
nc_genes_171to180 <- nc_genes %>% dplyr::slice(171:180)
files_171to180 <- nc_genes_171to180$GENEID

sequence_list_noncoding_171to180 <- list()

for(i in files_171to180){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_171to180[[symbol]] <- tmp
}

# Get sequences for gene 181 to 190
nc_genes_181to190 <- nc_genes %>% dplyr::slice(181:190)
files_181to190 <- nc_genes_181to190$GENEID

sequence_list_noncoding_181to190 <- list()

for(i in files_181to190){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_181to190[[symbol]] <- tmp
}

# Get sequences for gene 191 to 200
nc_genes_191to200 <- nc_genes %>% dplyr::slice(191:200)
files_191to200 <- nc_genes_191to200$GENEID

sequence_list_noncoding_191to200 <- list()

for(i in files_191to200){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_191to200[[symbol]] <- tmp
}

# Get sequences for gene 201 to 210
nc_genes_201to210 <- nc_genes %>% dplyr::slice(201:210)
files_201to210 <- nc_genes_201to210$GENEID

sequence_list_noncoding_201to210 <- list()

for(i in files_201to210){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_201to210[[symbol]] <- tmp
}

# Get sequences for gene 211 to 215
nc_genes_211to215 <- nc_genes %>% dplyr::slice(211:215)
files_211to215 <- nc_genes_211to215$GENEID

sequence_list_noncoding_211to215 <- list()

for(i in files_211to215){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_211to215[[symbol]] <- tmp
}

# Get sequences for gene 216 to 220
nc_genes_216to220 <- nc_genes %>% dplyr::slice(216:220) %>% dplyr::filter(GENEID != "ENSG00000251603")
files_216to220 <- nc_genes_216to220$GENEID

sequence_list_noncoding_216to220 <- list()

for(i in files_216to220){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_216to220[[symbol]] <- tmp
}

# Get sequences for gene 221 to 230
nc_genes_221to230 <- nc_genes %>% dplyr::slice(221:230)
files_221to230 <- nc_genes_221to230$GENEID

sequence_list_noncoding_221to230 <- list()

for(i in files_221to230){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_221to230[[symbol]] <- tmp
}

# Get sequences for gene 231 to 240
nc_genes_231to240 <- nc_genes %>% dplyr::slice(231:240)
files_231to240 <- nc_genes_231to240$GENEID

sequence_list_noncoding_231to240 <- list()

for(i in files_231to240){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_231to240[[symbol]] <- tmp
}

# Get sequences for gene 241 to 250
nc_genes_241to250 <- nc_genes %>% dplyr::slice(241:250)
files_241to250 <- nc_genes_241to250$GENEID

sequence_list_noncoding_241to250 <- list()

for(i in files_241to250){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_241to250[[symbol]] <- tmp
}

# Get sequences for gene 251 to 260
nc_genes_251to260 <- nc_genes %>% dplyr::slice(251:260)
files_251to260 <- nc_genes_251to260$GENEID

sequence_list_noncoding_251to260 <- list()

for(i in files_251to260){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_251to260[[symbol]] <- tmp
}

# Get sequences for gene 261 to 270
nc_genes_261to270 <- nc_genes %>% dplyr::slice(261:270) %>% dplyr::filter(GENEID != "ENSG00000223414")
files_261to270 <- nc_genes_261to270$GENEID

sequence_list_noncoding_261to270 <- list()

for(i in files_261to270){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_261to270[[symbol]] <- tmp
}

# Get sequences for gene 271 to 280
nc_genes_271to280 <- nc_genes %>% dplyr::slice(271:280) %>% dplyr::filter(GENEID != "ENSG00000235475")
files_271to280 <- nc_genes_271to280$GENEID

sequence_list_noncoding_271to280 <- list()

for(i in files_271to280){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_271to280[[symbol]] <- tmp
}

# Get sequences for gene 281 to 290
nc_genes_281to290 <- nc_genes %>% dplyr::slice(281:290) %>% dplyr::filter(GENEID != "ENSG00000226380")
files_281to290 <- nc_genes_281to290$GENEID

sequence_list_noncoding_281to290 <- list()

for(i in files_281to290){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_281to290[[symbol]] <- tmp
}

# Get sequences for gene 291 to 300
nc_genes_291to300 <- nc_genes %>% dplyr::slice(291:300) %>% dplyr::filter(GENEID != "ENSG00000281657")
files_291to300 <- nc_genes_291to300$GENEID

sequence_list_noncoding_291to300 <- list()

for(i in files_291to300){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_291to300[[symbol]] <- tmp
}

# Get sequences for gene 301 to 310
nc_genes_301to310 <- nc_genes %>% dplyr::slice(301:310)
files_301to310 <- nc_genes_301to310$GENEID

sequence_list_noncoding_301to310 <- list()

for(i in files_301to310){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_301to310[[symbol]] <- tmp
}

# Get sequences for gene 311 to 317
nc_genes_311to317 <- nc_genes %>% dplyr::slice(311:317)
files_311to317 <- nc_genes_311to317$GENEID

sequence_list_noncoding_311to317 <- list()

for(i in files_311to317){
    tmp <- biomaRt::getSequence(id = i, type="ensembl_gene_id", mart = ensembl, seqType = "transcript_exon_intron")
    tmp_symbol <- nc_genes %>% dplyr::filter(GENEID == i) %>% dplyr::select(GENENAME)
    symbol <- as.character(as.vector(tmp_symbol[1,]))
    name <- i
    tmp$GENENAME <- symbol
    sequence_list_noncoding_311to317[[symbol]] <- tmp
}

################################################################################################################################################
#######################################################     COMBINE LISTS      #################################################################
################################################################################################################################################

sequence_list_noncoding_all <- c(sequence_list_noncoding_1to30,
  sequence_list_noncoding_31to60,
  sequence_list_noncoding_61to70,
  sequence_list_noncoding_71to90,
  sequence_list_noncoding_91to110,
  sequence_list_noncoding_111to130,
  sequence_list_noncoding_131to150,
  sequence_list_noncoding_151to160,
  sequence_list_noncoding_161to170,
  sequence_list_noncoding_171to180,
  sequence_list_noncoding_181to190,
  sequence_list_noncoding_191to200,
  sequence_list_noncoding_201to210,
  sequence_list_noncoding_211to215,
  sequence_list_noncoding_216to220,
  sequence_list_noncoding_221to230,
  sequence_list_noncoding_231to240,
  sequence_list_noncoding_241to250,
  sequence_list_noncoding_251to260,
  sequence_list_noncoding_261to270,
  sequence_list_noncoding_271to280,
  sequence_list_noncoding_281to290,
  sequence_list_noncoding_291to300,
  sequence_list_noncoding_301to310,
  sequence_list_noncoding_311to317
)

# Export combined list as .Rdata file
save(sequence_list_noncoding_all, file="11 Expression metrics/Fibers at rest/Non-coding RNA/Sequences_noncoding.RData")



