
library(Seurat)
library(zellkonverter)
library(dplyr)
library(ggplot2)
library(scater)
library(slingshot)

path <- "/Users/liuy39/Research/H026/scRNAseq/"
meta <- read.table(paste0(path,"metadata/scRNAseq_meta_summary.txt"), header=T)
dir_list <- list.dirs(paste0(path,"cellranger_out"), recursive = FALSE)
dir_list <- paste0(dir_list,"/filtered_feature_bc_matrix.h5")
all(file.exists(dir_list))

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

h5ad  <- "colon_macrophages.h5ad"
sce   <- readH5AD(h5ad)
meta_macrophage <- colData(sce)
meta_macrophage <- data.frame(meta_macrophage)

sobj_list     <- list()
sample_list   <- list()

for (h5_file in dir_list) {
  sample1     <- split_path(h5_file)[2]
  count_data  <- Read10X_h5(filename = h5_file)
  seurat_obj  <- CreateSeuratObject(counts = count_data, project=sample1, min.cells = 3, min.features = 200)
  metadata    <- seurat_obj@meta.data
  metadata$Sex     <- meta[ meta$Animal_ID == sample1, "Sex"]
  metadata$Genoty  <- meta[ meta$Animal_ID == sample1, "Genoty"]
  metadata$DOB  <- meta[ meta$Animal_ID == sample1, "DOB"]
  metadata$batch  <- meta[ meta$Animal_ID == sample1, "batch"]
  metadata$tissue  <- meta[ meta$Animal_ID == sample1, "tissue"]
  metadata$macrophage <- c("noMacro", "macrophage")[1+as.numeric(rownames(metadata) %in% rownames(meta_macrophage))]  
  seurat_obj@meta.data <- metadata
  print(seurat_obj)
  sample_list <- append(sample_list,sample1 )
  sobj_list <- append(sobj_list, seurat_obj)
}

merged_seurat <- merge(x = sobj_list[[1]], y = sobj_list[2:length(sobj_list)], add.cell.ids = as.vector(sample_list))
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")

colon <- subset(merged_seurat, subset = (tissue == "colon"))
macrophage <- subset(merged_seurat, subset = (macrophage == "macrophage") & (tissue == "colon")) 
obj <- macrophage

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# without integration
#obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
#obj <- FindClusters(obj, resolution = 0.5, cluster.name = "unintegrated_clusters")
#obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
#DimPlot(obj, reduction = "umap.unintegrated", group.by = "seurat_clusters")
#DimPlot(obj, reduction = "umap.unintegrated", group.by = c("orig.ident","seurat_clusters")) 

# integration
obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.2)
obj <- RunUMAP(obj, dims = 1:30, reduction = "integrated.cca")
# Visualization
#DimPlot(obj, reduction = "umap", label=T, group.by = c("orig.ident", "seurat_clusters"))
DimPlot(obj, reduction = "umap", label=T, group.by =  "seurat_clusters")   + NoLegend()
VlnPlot(obj, features = c("Ccr2","H2-Ab1","Ly6c1","Ly6c2"), ncol=2 ) 
FeaturePlot(obj, features = c("Ccr2","H2-Ab1","Ly6c1","Ly6c2"), ncol=2 ) 

saveRDS(obj, file="colon_macrophage.rds") 

#slingshot for clusters 0 - 4
obj_slingshot1 <- subset(obj, subset = (seurat_clusters %in% c(seq(0,4) ))) 
obj_slingshot1$seurat_clusters2 <- paste0("umap_cluster_", obj_slingshot1$seurat_clusters)
obj.sce1 <- as.SingleCellExperiment(obj_slingshot1)

#sce.sling2 <- slingshot(obj.sce1, clusterLabels='seurat_clusters2', reducedDim='PCA')  
sce.sling2 <- slingshot(obj.sce1, clusterLabels='seurat_clusters2', start.clus="umap_cluster_2", reducedDim='PCA')
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)

# list lineages
SlingshotDataSet(sce.sling2 )
# class: SlingshotDataSet 
# Samples Dimensions
#lineages: 3 
#Lineage1: umap_cluster_2  umap_cluster_4  umap_cluster_3  
#Lineage2: umap_cluster_2  umap_cluster_4  umap_cluster_1  
#Lineage3: umap_cluster_2  umap_cluster_0  

obj.sce1 <- runUMAP(obj.sce1 , dimred="PCA")
reducedDim(sce.sling2, "UMAP") <- reducedDim(obj.sce1 , "UMAP")
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
# gg <- plotUMAP(sce.sling2, colour_by=I(sce.sling2$seurat_clusters)) 
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)

for (path in embedded) {
     embedded <- data.frame(path$s[path$ord,])
     gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
 }
gg



# change clusters name/label
new.cluster.ids <- c(paste0("cluster_", seq(0,11)))
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
DimPlot(obj, reduction = "umap", label=T ) + NoLegend()


new.cluster.ids <- c(paste0("umap_cluster_", seq(0,11)))
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
levels(obj)
DimPlot(obj, reduction = "umap", label=T ) + NoLegend()


