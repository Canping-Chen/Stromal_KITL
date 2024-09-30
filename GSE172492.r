rm(list=ls())
gc()
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
set.seed(12)

data <- Read10X_h5('GSM5354107_RS014.h5')
sce1 <- CreateSeuratObject(counts = data, project = 'Normal')
Normal <- sce1 %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", verbose = F) %>% 
    ScaleData(verbose = F) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters()  


data <- Read10X_h5('GSM5354108_RS015.h5')
sce2 = CreateSeuratObject(counts = data, project = 'Injury')
Injury <- sce2 %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", verbose = F) %>% 
    ScaleData(verbose = F) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters()  


Samples <- c('Normal', 'Injury')
merged_object <- merge(Normal, Injury, 
                        add.cell.ids = Samples, 
                        project = "Merged_Project")

sce <- merged_object %>%
            NormalizeData(verbose = F) %>%
            FindVariableFeatures(selection.method = "vst", verbose = F) %>% 
            ScaleData(verbose = F) %>%
            RunPCA() %>%
            FindNeighbors(dims = 1:30) %>%
            RunUMAP(dims = 1:30) %>%
            FindClusters() 

sce$orig.ident <- factor( sce$orig.ident, levels = c( 'Normal', 'Injury'))


FeaturePlot(object = sce, 
            features = c("Lrat", "Col1a1", "Pecam1", "Vwf", "Kdr",  "Cdh5"))

