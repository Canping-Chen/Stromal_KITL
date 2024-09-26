rm(list=ls())
gc()
require(Seurat)
require(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(tidydr)
library(Cairo)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(future)
library(future.apply)
set.seed(12)

main_dir = "Stromal_KITL/"
data_dir = '00.data/'
output_dir <- file.path(main_dir, '01.batch/')

if (!dir.exists(output_dir)){
    dir.create(output_dir,recursive = TRUE)
} else {
    print("Dir already exists!")
}        

folders <- list.dirs(file.path(main_dir, data_dir), recursive = FALSE)
samples <- c('Tumor', 'Normal')
sceList = lapply(folders,function(pro){ 
    sample_name <- basename(pro)
    ind_10X = Read10X(pro)         
    ind_sce = CreateSeuratObject(counts =  ind_10X ,
                    project <- sample_name)                        
        return(ind_sce)
    })

int_genes <- Reduce(intersect, list(rownames(sceList[[1]]) , rownames(sceList[[2]])))
sceList[[1]] <- sceList[[1]][int_genes,]
sceList[[1]]$condition <- samples[1]
sceList[[2]] <- sceList[[2]][int_genes,]
sceList[[2]]$condition <- samples[2]    

do.call(rbind,lapply(sceList, dim))    
names(sceList) = samples

sce.all <- merge(sceList[[1]], y= sceList[ -1 ] ,
        add.cell.ids=samples) 
table(sce.all$orig.ident)

scRemoveBatch = function(seurat.data, 
                        batchID, 
                        n.pcs = 50){
    library(Seurat)
    library(harmony)
    seurat_int <- seurat.data %>% 
        SCTransform( verbose = FALSE) %>%
        RunPCA(npcs = n.pcs, verbose = F) %>%
        RunHarmony(batchID, plot_convergence = T) %>%             
        RunUMAP(reduction = "harmony", dims = 1:20, verbose = F)
    return(seurat_int)
}

seurat.harmony <- scRemoveBatch(seurat.data = sce.all, batchID = "orig.ident")
seurat.harmony$condition <- sub("_.*", "", colnames(seurat.harmony))

