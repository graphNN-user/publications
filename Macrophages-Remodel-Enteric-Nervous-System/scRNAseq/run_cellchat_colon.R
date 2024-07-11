
library(Seurat)
library(zellkonverter)
library(dplyr)
library(ggplot2)
library(scater)
library(slingshot)
library(CellChat)
library(patchwork)

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

h5ad  <- "colon_neuron.h5ad"
sce   <- readH5AD(h5ad)
meta_neuron <- colData(sce)
meta_neuron <- data.frame(meta_neuron)

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
  metadata$neuron     <- c("noNeuron", "neuron")[1+as.numeric(rownames(metadata) %in% rownames(meta_neuron))]
  seurat_obj@meta.data <- metadata
  print(seurat_obj)
  sample_list <- append(sample_list,sample1 )
  sobj_list <- append(sobj_list, seurat_obj)
}

merged_seurat <- merge(x = sobj_list[[1]], y = sobj_list[2:length(sobj_list)], add.cell.ids = as.vector(sample_list))
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")

colon <- subset(merged_seurat, subset = (macrophage == "macrophage" | neuron == "neuron" ) & (tissue == "colon"))

obj <- colon
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# integration
obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.2)
obj <- RunUMAP(obj, dims = 1:30, reduction = "integrated.cca")

macrophage <- readRDS("colon_macrophage.rds") 
neuron <- readRDS("colon_neuron.rds") 
meta_macro <- macrophage@meta.data
meta_macro <- dplyr::rename(meta_macro, cell_type=macrophage)
meta_macro$seurat_clusters0 <- paste0("macrophage_cluster_", meta_macro$seurat_clusters) 
meta_neuron <- neuron@meta.data
meta_neuron <- dplyr::rename(meta_neuron, cell_type=neuron)
#meta_neuron$seurat_clusters <- paste0("neuron_cluster_", meta_neuron$seurat_clusters)
meta_neuron$seurat_clusters0 <- "neuron_cluster"


meta  <- rbind(meta_macro, meta_neuron) 
meta1 <- meta[, "seurat_clusters0",drop=F] 
meta0 <- obj@meta.data
meta1 <- meta1[rownames(meta0), , drop=F]
#meta$samples <- meta$orig.ident

meta2 <- merge(meta0, meta1, by=0)
rownames(meta2) <- meta2$Row.names
meta2 <- meta2[,-1]

obj@meta.data <- meta2
#selected <- c("macrophage_0","macrophage_1","macrophage_2","macrophage_3","macrophage_4","macrophage_5","neuron_0","neuron_1","neuron_2","neuron_3","neuron_4")

#selected <- c("macrophage_cluster_0", "macrophage_cluster_1",  "macrophage_cluster_2","macrophage_cluster_3","macrophage_cluster_4","neuron_cluster")
#selected <- c("macrophage_cluster_0","macrophage_cluster_2","neuron_cluster_0","neuron_cluster_1","neuron_cluster_2","neuron_cluster_3","neuron_cluster_4")
selected <- c("macrophage_cluster_0", "macrophage_cluster_2","neuron_cluster")
obj1 <- subset(obj, subset = (seurat_clusters0 %in% selected ))

#data.input <- obj1[["RNA"]]$data
saveRDS(obj1, file="macrophage_02_clusters_neuron_colon_seuratobj.rds")

#obj1 <- readRDS("macrophage_cluster02_neuron_colon_seuratobj.rds") 
cellchat <- createCellChat(object = obj1, group.by = "seurat_clusters0", assay = "RNA")
CellChatDB <- CellChatDB.mouse
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 2447 

#ptm = Sys.time()
#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))
#> [1] 13.20763
cellchat <- computeCommunProb(cellchat, type = "triMean")  # "triMean", "truncatedMean", "thresholdedMean", "median"
cellchat <- filterCommunication(cellchat, min.cells = 100)
df.net <- subsetCommunication(cellchat)
head(df.net)
#df.interest <- df.net[ (grepl("macrophage", df.net$source) & grepl("neuron", df.net$target)) | (grepl("macrophage", df.net$target) & grepl("neuron", df.net$source)), ] 

#macrophage <- c("macrophage_0","macrophage_1","macrophage_2","macrophage_3","macrophage_4","macrophage_5") 
macrophage <- c("macrophage_0","macrophage_2")
neuron <- c("neuron_0","neuron_1","neuron_2","neuron_3","neuron_4")
df.net1 <- subsetCommunication(cellchat, sources.use = macrophage, targets.use = neuron)
#df.net2 <- subsetCommunication(cellchat, sources.use = neuron, targets.use =macrophage  )

#calculates the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
#visualize the aggregated cell-cell communication network.
groupSize <- as.numeric(table(cellchat@idents))
pdf("cell_cell_communication.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

mat <- cellchat@net$weight

#par(mfrow = c(2,4), xpd=TRUE)
#for (i in 1:nrow(mat)) {
for (i in 1:2) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max =
  max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf("significant_interactions.pdf")
# show all the significant interactions (L-R pairs)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(3:7), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(3:7), remove.isolate = FALSE) 
dev.off()

netVisual_chord_gene(cellchat, sources.use = c(1:2), targets.use = 3,small.gap=0.1, lab.cex = 0.5,legend.pos.y = 5)
netVisual_chord_gene(cellchat, sources.use = 3, targets.use = c(1:2),small.gap=0.1, lab.cex = 0.5,legend.pos.y = 5)

library(dplyr)

df.net_neuron_macro0 <- df.net[(df.net$source == "neuron_cluster" & df.net$target == "macrophage_cluster_0")|(df.net$target == "neuron_cluster" & df.net$source == "macrophage_cluster_0"),]  
df.net_neuron_macro0b <- arrange(df.net_neuron_macro0, desc(prob))
df.net_neuron_macro0c <- head(df.net_neuron_macro0b, 25)
netVisual_chord_gene(cellchat, net=df.net_neuron_macro0c)
pdf("neuron_macro0.pdf")
netVisual_chord_gene(cellchat, net=df.net_neuron_macro0c, small.gap=0.1, lab.cex = 0.5,legend.pos.y = 5) 
dev.off()


df.net_neuron_macro2 <- df.net[(df.net$source == "neuron_cluster" & df.net$target == "macrophage_cluster_2")|(df.net$target == "neuron_cluster" & df.net$source == "macrophage_cluster_2"),]
df.net_neuron_macro2b <- arrange(df.net_neuron_macro2, desc(prob))
df.net_neuron_macro2c <- head(df.net_neuron_macro2b, 25)
netVisual_chord_gene(cellchat, net=df.net_neuron_macro2c)
pdf("neuron_macro2.pdf")
netVisual_chord_gene(cellchat, net=df.net_neuron_macro2c, small.gap=0.1, lab.cex = 0.5,legend.pos.y = 5)
dev.off()




saveRDS(cellchat, file="macrophage_neuron_colon_cellchatobj.rds")




