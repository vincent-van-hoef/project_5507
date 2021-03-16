# Load packages
library("Seurat")
library("cowplot")
library("dplyr")
library("ggplot2")
library("cowplot")
library("gridExtra")
library("patchwork")

# Set working directory to script location (when sourcing...)
#proj_dir <- dirname(sys.frame(1)$ofile)
# or interactively
proj_dir <- "/Users/vinva957/Desktop/NBIS/Projects/project_5507"
setwd(proj_dir)

# Load several custom functions
# source("helpers.R")

# Set up a report folder, remove existing one first
res_dir <- paste0(proj_dir, "/Report/")
unlink(res_dir, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE)

#############
# Load Data #
#############

# Load datasets after the mapping pipeline for day1 and day24 separately and split into the different modalities
d1 <- Read10X("/Users/vinva957/Desktop/NBIS/Projects/project_5507/data/day1_data/filtered_feature_bc_matrix/")
d1_gene <- CreateSeuratObject(counts = d1$`Gene Expression`, project = "D1_Gene")
d1_adt <- CreateSeuratObject(counts = d1$`Antibody Capture`, project = "D1_ADT")

d24 <- Read10X("/Users/vinva957/Desktop/NBIS/Projects/project_5507/data/day24_data/")
d24_gene <- CreateSeuratObject(counts = d24$`Gene Expression`, project = "D24_Gene")
d24_adt <- CreateSeuratObject(counts = d24$`Antibody Capture`, project = "D24_ADT")

# Preprocess (qc, normalize, merge days, ...) each modality spearately, then integrate.

#######
# RNA #
#######

# Create a QC folder in the Report dir
qc_rna_dir <- paste0(res_dir, "QC/RNA/")
dir.create(qc_rna_dir, showWarnings = FALSE, recursive = TRUE)

# Calculate the mitochondrial gene percentage
d1_gene[["percent.mt"]] <- PercentageFeatureSet(d1_gene, pattern = "^mt-")
d24_gene[["percent.mt"]] <- PercentageFeatureSet(d24_gene, pattern = "^mt-")

# Visualize the descriptive statistics
p1 <- VlnPlot(d1_gene, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, combine = FALSE)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + NoLegend()
}

p2 <- VlnPlot(d24_gene, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, combine = FALSE)
for(i in 1:length(p2)) {
  p2[[i]] <- p2[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + NoLegend()
}

pdf(paste0(qc_rna_dir, "Stats_day1_24.pdf"))
cowplot::plot_grid(plotlist = c(p1, p2), nrow = 2, labels = c("A", "", "", "B", "",""))
dev.off()

# Filter to have between 200 and 4000 genes, less than 20000 molecules anda mito percentage less than 20%
d1_gene <- subset(d1_gene, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 20)
d24_gene <- subset(d24_gene, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000  & percent.mt < 20)

# Normalize the data using the logNormalize method
d1_gene <- NormalizeData(d1_gene, normalization.method = "LogNormalize")
d1_gene <- FindVariableFeatures(d1_gene, selection.method = "vst", nfeatures = 2000)

d24_gene <- NormalizeData(d24_gene, normalization.method = "LogNormalize")
d24_gene <- FindVariableFeatures(d24_gene, selection.method = "vst", nfeatures = 2000)

# Merge the two days naively to show need for integration 
pca.combo <- merge(d1_gene, y = d24_gene, add.cell.ids = c("d1", "d24"))
# Run the standard workflow for visualization and clustering
pca.combo <- FindVariableFeatures(pca.combo)
pca.combo <- ScaleData(pca.combo, verbose = FALSE)
pca.combo <- RunPCA(pca.combo, npcs = 30, verbose = FALSE, reduction.name = "comboPCA")
pca.combo <- RunUMAP(pca.combo, reduction = "comboPCA", dims = 1:20)

# More sophisticated integration
rna.anchors <- FindIntegrationAnchors(object.list = list(d1_gene, d24_gene), dims = 1:20)
gene.combined <- IntegrateData(anchorset = rna.anchors, dims = 1:20, new.assay.name = "integrated.RNA")

DefaultAssay(gene.combined) <- "integrated.RNA"

# Run the standard workflow for visualization and clustering
gene.combined <- ScaleData(gene.combined, verbose = FALSE)
gene.combined <- RunPCA(gene.combined, npcs = 30, verbose = FALSE, reduction.name = "pca.RNA")
# t-SNE and Clustering
gene.combined <- RunUMAP(gene.combined, reduction = "pca.RNA", dims = 1:20, reduction.name = "umap.RNA")

# Visualization
p1 <- DimPlot(pca.combo, reduction = "umap", group.by = "orig.ident") + ggtitle("") + NoLegend() +FontSize(
  x.text = 6,
  y.text = 6,
  x.title = 6,
  y.title = 6,
  main = NULL
)
p2 <- DimPlot(gene.combined, reduction = "umap.RNA", group.by = "orig.ident") + ggtitle("") + theme(legend.position = "bottom")  +FontSize(
  x.text = 6,
  y.text = 6,
  x.title = 6,
  y.title = 6,
  main = NULL
)
pdf(paste0(qc_rna_dir, "Merging_vs_Integration_RNA.pdf"))
plot_grid(p1, p2, ncol = 1, labels = c("A", "B"))
dev.off()


#######
# ADT #
#######

# Create a QC folder in the Report dir
qc_adt_dir <- paste0(res_dir, "QC/ADT/")
dir.create(qc_adt_dir, showWarnings = FALSE, recursive = TRUE)

# Normalize the ADT data according to the CLR method
d1_adt <- NormalizeData(d1_adt, normalization.method = "CLR", margin = 2)
d24_adt <- NormalizeData(d24_adt, normalization.method = "CLR", margin = 2)

# Run rest of workflow using all ADT as variable features
features <- rownames(d1_adt)
# Day 1
d1_adt <- ScaleData(d1_adt, features = features)
d1_adt <- RunPCA(d1_adt, features = features, verbose = FALSE, npcs = 10, nfeatures.print = 5, approx = FALSE)
# Day 24
d24_adt <- ScaleData(d24_adt, features = features)
d24_adt <- RunPCA(d24_adt, features = features, verbose = FALSE, npcs = 10, nfeatures.print = 5, approx = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(d1_adt, d24_adt), 
                                  anchor.features = features, 
                                  scale = FALSE,
                                  normalization.method = "LogNormalize",
                                  reduction = "rpca",
                                  dims = 1:9)
adt.combined <- IntegrateData(anchorset = anchors, new.assay.name = "integrated.ADT", dims = 1:9)
adt.combined <- ScaleData(adt.combined, verbose = FALSE, features = features)
adt.combined <- RunPCA(adt.combined, verbose = FALSE, reduction.name = "pca.ADT")
adt.combined <- RunUMAP(adt.combined, dims = 1:9, reduction = "pca.ADT", reduction.name = "umap.ADT")

adt.combined <- adt.combined[,colnames(gene.combined)]

pdf(paste0(qc_adt_dir, "Integrated_ADT.pdf"))
DimPlot(adt.combined)
dev.off()

############
# Integrate#
############

# Create a QC folder in the Report dir
wnn_dir <- paste0(res_dir, "Clustering/")
dir.create(wnn_dir, showWarnings = FALSE, recursive = TRUE)

all.combined <- gene.combined

all.combined[["integrated.ADT"]] <- adt.combined[["integrated.ADT"]]
all.combined[["pca.ADT"]] <- adt.combined[["pca.ADT"]]
all.combined[["umap.ADT"]] <- adt.combined[["umap.ADT"]]

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
all.combined <- FindMultiModalNeighbors(
  all.combined, 
  reduction.list = list("pca.RNA", "pca.ADT"), 
  dims.list = list(1:30, 1:9), 
  modality.weight.name = "RNA.weight")
all.combined <- RunUMAP(all.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

resolutions <- c(0.05,0.15,0.25,0.5,1,2)

for(resolution in resolutions){
  
  # Create a folder for each resolution
  resol_dir <- paste0(wnn_dir, "Resolution_", gsub("\\.", "_", resolution), "/")
  dir.create(resol_dir, showWarnings = FALSE, recursive = TRUE)
    
  all.combined <- FindClusters(all.combined, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)

  pdf(paste0(resol_dir, "Clustering.pdf"))
  p1 <- DimPlot(all.combined, reduction = 'wnn.umap', group.by = paste0("wsnn_res.", resolution)) + 
    plot_annotation(title = paste0("Resolution_", gsub("\\.","_", resolution)))
  print(p1)
  dev.off()
  
  marker_features <-  c("Itgb7","Cd34","Il10rb","Fcgr1","Cd44","P2rx7","Fcgr3","Kit","Itgam", "Itgae")  
  names(marker_features) <- features
  plotList <- list()
  for (adt in names(marker_features)){
    DefaultAssay(all.combined) <- "integrated.ADT"
    p1 <- FeaturePlot(all.combined, reduction = "wnn.umap", features = adt, label=TRUE)
    p2 <- VlnPlot(all.combined, features = adt)
    DefaultAssay(all.combined) <- "integrated.RNA"
    p3 <- FeaturePlot(all.combined, reduction = "wnn.umap", features = marker_features[adt], label=TRUE)
    p4 <- VlnPlot(all.combined, features = marker_features[adt])
    plotList[[adt]] <- list(p1,p2, p3, p4)
  }
  pdf("Draft_Clustering.pdf", onefile=TRUE)
  for(i in 1:length(marker_features)){
    grid.arrange(plot_grid(plotlist = plotList[[i]], nrow=2))
  } 
}




# Basophils: c-kit- CD11b+ integrin B7 lo/negative
# Protein
marker_features <-  c("Itgb7","Cd34","Il10rb","Fcgr1","Cd44","P2rx7","Fcgr3","Kit","Itgam", "Itgae")  
names(marker_features) <- features
plotList <- list()
for (adt in names(marker_features)){
  DefaultAssay(all.combined) <- "integrated.ADT"
  p1 <- FeaturePlot(all.combined, reduction = "wnn.umap", features = adt, label=TRUE)
  p2 <- VlnPlot(all.combined, features = adt)
  DefaultAssay(all.combined) <- "integrated.RNA"
  p3 <- FeaturePlot(all.combined, reduction = "wnn.umap", features = marker_features[adt], label=TRUE)
  p4 <- VlnPlot(all.combined, features = marker_features[adt])
  plotList[[adt]] <- list(p1,p2, p3, p4)
}
pdf("Draft_Clustering.pdf", onefile=TRUE)
for(i in 1:length(marker_features)){
  grid.arrange(plot_grid(plotlist = plotList[[i]], nrow=2))
} 
dev.off()


marker_features <-  c("Tpsb2","Cpa3","Mcpt4","Mcpt1","Srgn","Hpgds","Mcpt8","Il3ra")  
plotList <- list()
for (marker in marker_features){
  DefaultAssay(all.combined) <- "integrated.RNA"
  p1 <- FeaturePlot(all.combined, reduction = "wnn.umap", features = marker, label=TRUE)
  p2 <- VlnPlot(all.combined, features = marker)
  plotList[[marker]] <- list(p1,p2)
}
pdf("Draft_Clustering_v2.pdf", onefile = TRUE)
for(i in c(1,3,5,7)){
  grid.arrange(
    plot_grid(plotlist = plotList[[i]], nrow=1),
    plot_grid(plotlist = plotList[[i+1]], nrow=1))
} 
dev.off()




