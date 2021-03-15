# Load packages
library("Seurat")
library("cowplot")
library("dplyr")
library("ggplot2")
library("cowplot")
library("gridExtra")


# Load datasets
d1 <- Read10X("/Users/vinva957/Desktop/NBIS/Projects/project_5507/data/day1_data/filtered_feature_bc_matrix/")
d1_gene <- CreateSeuratObject(counts = d1$`Gene Expression`, project = "D1_Gene")
d1_adt <- CreateSeuratObject(counts = d1$`Antibody Capture`, project = "D1_ADT")

d24 <- Read10X("/Users/vinva957/Desktop/NBIS/Projects/project_5507/data/day24_data/")
d24_gene <- CreateSeuratObject(counts = d24$`Gene Expression`, project = "D24_Gene")
d24_adt <- CreateSeuratObject(counts = d24$`Antibody Capture`, project = "D24_ADT")

