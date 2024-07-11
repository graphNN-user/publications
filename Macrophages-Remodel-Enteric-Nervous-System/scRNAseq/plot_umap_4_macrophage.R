
library(Seurat)
library(zellkonverter)
library(dplyr)
library(ggplot2)
library(scater)

markers <- read.delim("marker_genes4.txt", sep="\t", header=F)
markers <- markers$V1

macrophages <- readRDS("colon_macrophage.rds") 
obj <- macrophages

#subsets <- c(paste0("subset_", c(4,0,2,3,"0+4")))
subsets <- c(paste0("macro", seq(1,5))) 
#new.cluster.ids <- c(subsets, "dendritic cells", "unknown","dendritic cells","dendritic cells","dendritic cells" ) 
new.cluster.ids <- c(subsets, "dc1", "unknown","dc2","dc3","dc4" )

#png("macrophage_clusters_1.png", width = 480, height = 480, units = "px", pointsize = 12, bg = "white",  res = 600)
#png("macrophage_clusters_1.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4 )
#png("macrophage_clusters_1.png", width = 5, height = 5, units = "in", res = 1200, pointsize = 1 )

png("macrophage_umap.png", width = 6, height = 6, units = "in", res = 1800, pointsize = 2 )

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
levels(obj)
DimPlot(obj, reduction = "umap", label=F )  + NoLegend()
dev.off()


for ( i in seq(1,12) ) {
  tiff(paste0("marker_",i,".tiff"), width = 7, height = 6, units = "in", res = 800, pointsize = 1, compression = "lzw" )
  temp <-  FeaturePlot(obj, features=markers[[i]], label=F, label.size = 2, alpha=0.3)
  print(temp)
  dev.off()
}

png("macrophage_clusters_2.png", width = 5, height = 5, units = "in", res = 1200, pointsize = 1 )
# names(new.cluster.ids) <- levels(obj)
# obj <- RenameIdents(obj, new.cluster.ids)
# levels(obj)
DimPlot(obj, reduction = "umap", label=T ) + NoLegend()
dev.off()

png("markergenes_1.png", width = 7, height = 7, units = "in", res = 1200, pointsize = 1 )
FeaturePlot(obj, features=markers[1:4], ncol=2, label=F, label.size = 2, alpha=0.3)
dev.off()
png("markergenes_2.png", width = 7, height = 7, units = "in", res = 1200, pointsize = 1 )
FeaturePlot(obj, features=markers[5:8], ncol=2, label=F, label.size = 2, alpha=0.3)
dev.off()
png("markergenes_3.png", width = 7, height = 7, units = "in", res = 1200, pointsize = 1 )
FeaturePlot(obj, features=markers[9:12], ncol=2, label=F, label.size = 2, alpha=0.3)
dev.off()
png("markergenes_4.png", width = 7, height = 7, units = "in", res = 1200, pointsize = 1 )
FeaturePlot(obj, features=markers[13:16], ncol=2, label=F, label.size = 2, alpha=0.3)
dev.off()



pdf("macrophage_clusters.pdf")
DimPlot(obj, reduction = "umap", label=T, group.by = c("seurat_clusters"), raster.dpi = c(512, 512)) 
dev.off()

pdf("macrophage_subsets.pdf")
DimPlot(obj, reduction = "umap", label=T, group.by = c("ident"), raster.dpi = c(512, 1024)) + NoLegend() 
dev.off()

pdf("marker_genes.pdf")
FeaturePlot(obj, features=markers[1:4], ncol=2, label=T, label.size = 2, alpha=0.3)
FeaturePlot(obj, features=markers[5:8], ncol=2, label=T, label.size = 2, alpha=0.3)
FeaturePlot(obj, features=markers[9:12], ncol=2, label=T, label.size = 2, alpha=0.3)
FeaturePlot(obj, features=markers[13:16], ncol=2, label=T, label.size = 2, alpha=0.3)
dev.off()


macrophage1 <- subset(macrophages, subset= (seurat_clusters %in% c(0,1,2,3,4) ) ) 
obj <- macrophage1
new.cluster.ids <- c(paste0("water_subset_", c(4,0,2,3,"0+4")))

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
levels(obj)
DimPlot(obj, reduction = "umap", label=T ) + NoLegend()



