rm(list=ls())
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(clustree)
library(cowplot)
library(stringr)
library(monocle3)
library(viridis)
main_dir = "Stromal_KITL/"
data_dir = '01.batch/'
output_dir <- file.path(main_dir, '02.monocle3/')

data <- GetAssayData(mouse_data, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
cds <- new_cell_data_set(data,
                        cell_metadata = mouse_data@meta.data,
                        gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
int.embed <- Embeddings(mouse_data, reduction = "harmony")
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$PCA),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds, resolution=0.005,reduction_method = 'UMAP')

plot_cells(cds, show_trajectory_graph = FALSE) + 
        ggtitle("label by clusterID")
plot_cells(cds, color_cells_by = "condition", show_trajectory_graph = FALSE) + 
        ggtitle("label by condition")

cds <- learn_graph(cds)

plot_cells(cds,
        color_cells_by = "condition",
        label_groups_by_cluster=T,
        label_leaves=T,
        label_branch_points=T)

cds <- order_cells(cds) 

p8 <- plot_cells(cds,
        genes=c("Kitl"),           
        label_cell_groups=F,
        show_trajectory_graph=F, 
        cell_size=1, trajectory_graph_color = "black", 
        label_branch_points = F, 
        label_roots = F, label_leaves = F)+
        scale_color_viridis(option="inferno")



