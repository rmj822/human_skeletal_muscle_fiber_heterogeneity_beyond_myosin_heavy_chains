
################################################################################################################################################
################################################       PREPARATION      ########################################################################
################################################################################################################################################

library(Seurat)
library(devtools)
library(SeuratWrappers)
library(SeuratDisk)
library(monocle3)


# Set working directory to own folder -------------------------------------
setwd("~/OneDrive - UGent/PhD/Projects/2018 CHH pathway/8 Single fiber transcriptomics/Single fiber RNAseq")

# Load filtered Seurat object ---------------------------------------------
load("8 Fiber heterogeneity (only rested samples)/2 Reclustering/Reclustering/filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest.Rdata")

DefaultAssay(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "SCT"
Idents(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest) <- "final_cluster"


SaveH5Seurat(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest, filename = "~/Python-RNA-velocity/seuratobj.h5Seurat")
Convert("seuratobj.h5Seurat", dest = "~/Python-RNA-velocity/h5ad")


################################################################################################################################################
#################################################       PSEUDOTIME WITH MONOCLE3     ###########################################################
################################################################################################################################################

# Converting a Seurat object into monocle:

data_monocle_transcriptomics <-
    SeuratWrappers::as.cell_data_set(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest)

# Visualization of the object metadata:

head(SummarizedExperiment::colData(data_monocle_transcriptomics))

# Double checking that the gene names are there:

rownames(monocle3::fData(data_monocle_transcriptomics))[1:10]

# Double checking that the counts are there:

head(counts(data_monocle_transcriptomics))

# Double checking that the coLData has row names:

row.names(data_monocle_transcriptomics@int_colData)
row.names(data_monocle_transcriptomics@int_colData) <- rownames(filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@meta.data)
row.names(data_monocle_transcriptomics@int_colData)

row.names(data_monocle_transcriptomics@colData)
row.names(data_monocle_transcriptomics@metadata)
row.names(data_monocle_transcriptomics@int_metadata)

# Recreating object clusters:

recreate.partitions_transcriptomics <- c(rep(1, length(data_monocle_transcriptomics@colData@rownames)))

recreate.partitions_transcriptomics <-
    SummarizedExperiment::colData(data_monocle_transcriptomics) |>
    as.data.frame() |>
    dplyr::pull(final_cluster)

names(recreate.partitions_transcriptomics) <- data_monocle_transcriptomics@colData@rownames
recreate.partitions_transcriptomics <- as.factor(recreate.partitions_transcriptomics)
recreate.partitions_transcriptomics

data_monocle_transcriptomics@clusters@listData[["UMAP"]][["partitions"]] <-
    recreate.partitions_transcriptomics

# Assign cluster information

list.cluster_transcriptomics <- filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@active.ident

data_monocle_transcriptomics@clusters@listData[["UMAP"]][["clusters"]] <-
    list.cluster_transcriptomics

# Assign UMAP coordinates
data_monocle_transcriptomics@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-
    filtered_normalized_fibertype_reclustered_seurat_wo_MSTRG_rest@reductions$umap@cell.embeddings

# Showing cluster of cells before trajectory: check if UMAP is correct:

monocle3::plot_cells(
    data_monocle_transcriptomics,
    color_cells_by = "final_cluster",
    label_groups_by_cluster = F,
    group_label_size = 10,
    cell_size = 0.75
) +
    ggplot2::scale_color_viridis_d(option = "plasma") +
    ggplot2::theme(legend.position = "right") +
    ggplot2::theme_minimal()


# Calculate clusters, needed for learn_graph althoud not used
data_monocle_transcriptomics <- monocle3::cluster_cells(data_monocle_transcriptomics)

# Calculate pseudotime - learn trajectory:

data_monocle_transcriptomics <-
    monocle3::learn_graph(data_monocle_transcriptomics, use_partition = F)

# Plot pseudotime

monocle3::plot_cells(
    data_monocle_transcriptomics,
    color_cells_by = "final_cluster",
    label_groups_by_cluster = F,
    group_label_size = 7,
    cell_size = 1,
    label_cell_groups = FALSE,
    trajectory_graph_color = "black",
    alpha = 0.5,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE
)


# Identify root principal points
get_earliest_principal_node <- function(cds, time_bin="130-170"){
    cell_ids <- which(colData(cds)[, "final_cluster"] == time_bin)

    closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]

    root_pr_nodes
}

data_monocle_transcriptomics <- order_cells(data_monocle_transcriptomics, root_pr_nodes=get_earliest_principal_node(data_monocle_transcriptomics, time_bin = "Fast2"))

plot_cells(data_monocle_transcriptomics,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


data_monocle_transcriptomics <- order_cells(data_monocle_transcriptomics, reduction_method = "UMAP", root_cells = colnames(data_monocle_transcriptomics[, clusters(data_monocle_transcriptomics) == 2]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)


colnames(data_monocle_transcriptomics)







# ggplot2::ggsave(
#     here::here("doc/figures/figure_4/UMAP_pseudotime_proteomics.png"),
#     height = 60,
#     width = 90,
#     units = "mm"
# )

seurat_clusters <- seurat_proteome@meta.data |>
    as.data.frame() |>
    dplyr::select(fiberID, RNA_snn_res.0.4) |>
    dplyr::rename("seurat_clusters" = RNA_snn_res.0.4) |>
    tibble::remove_rownames()

# readr::write_csv(seurat_clusters,
#                  here::here("data/proteomics clustering/seurat_clusters_6PC_res04.csv"))

