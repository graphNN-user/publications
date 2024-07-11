
library(Seurat)
library(zellkonverter)
library(dplyr)
library(ggplot2)
library(scater)
library(slingshot)
library(scrattch.vis)
options(stringsAsFactors = F)

colon <- readRDS("colon_macrophage.rds")

obj <- colon
DimPlot(obj, reduction = "umap", label=T, group.by =  "seurat_clusters")   + NoLegend()
#markers <- read.table("marker_genes2.txt", header=F)
#obj1 <- subset(obj, subset= (seurat_clusters %in% c(seq(0,5))))
#VlnPlot(obj1, features=markers[1:3], ncol=3)  

sce <- as.SingleCellExperiment(obj)
sce$seurat_clusters2 <- factor( paste0("umap_cluster_", sce$seurat_clusters))

#colors <-  c("#ebcb2e", "#9ec22f", "#a9961b", "#cc3a1b", "#cc8778" , "#d14c8d", "#4cabdc", "#5ab793")
#colors <- c("#007500","#0094C2","#00B7FF","#00FFFF","#0640FF","#0DFF00","#104F00","#21F4CC","#282828","#505050")
#colors <- c("#579D53","#5F79C4","#6B03FF","#74A9CF","#74D16F","#787878","#7F5564","#840C56","#8800FF","#8C6BB1") 
#colors <- c("#8D72D1","#8E170E","#9400E2","#A0A0A0","#A1990F","#A6E2A3","#ADBF00","#B28D9A","#BB36BF","#C77963")
#      #C8C8C8 #C994C7 #CD2626 #CD6090 #CEFF00 #D9F0A3 #E54E03 #E7298A #EE2C2C #F08080 
#      #F800FF #F99EB4 #FDDCFF #FF3030 #FF6A6A #FF8C00 #FFA500 #FFD700 #FFFF00 

colors <-  c("#ebcb2e","#9ec22f","#a9961b","#cc3a1b","#cc8778","#d14c8d","#4cabdc","#5ab793","#bb36bf","#c77963") 

colorsident <- cbind(colors = colors, id = levels(sce$seurat_clusters2))

anno.df <- data.frame(sample_name = colnames(sce), primary_type_id = sce$seurat_clusters2, primary_type_label = sce$seurat_clusters2, primary_type_color = colorsident[match(sce$seurat_clusters2, colorsident[,2]), 1]) 

exprs <- as.data.frame(t(as.matrix(logcounts(sce))))
exprs <- cbind(rownames(exprs), exprs)
colnames(exprs)[1] <- "sample_name"

#genes <- gene_list$Gene[grep("Complement & coagulation cascade", gene_list$Annotation)]
#genes <- c("Ccr2","Ly6c1","Ly6c2")
#genes <- read.csv("markergene_3.csv", header=T)
genes <- read.table("addiitonalgene.txt", header=T, sep="\t") 
genes <- genes$markergene
genes <- genes$addiitonal_genes
exprs$H2_Ab1 <- exprs$"H2-Ab1" 
exprs$H2_Eb1 <- exprs$"H2-Eb1" 

pdf("sample_bar_colon.pdf")
sample_bar_plot(exprs, anno.df, genes = genes, grouping = "primary_type", bg_color ="#f7f7f7")
dev.off()


