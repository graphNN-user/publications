
library(Seurat)
library(zellkonverter)
library(dplyr)
library(ggplot2)
library(scater)
library(slingshot)

macrophages <- readRDS("colon_macrophage.rds") 
macrophage1 <- subset(macrophages, subset= (seurat_clusters %in% c(0,1,2,3,4) ) ) 
obj <- macrophage1
#new.cluster.ids <- c(paste0("subset_", c(4,0,2,3,"0+4")))
new.cluster.ids <- c(paste0("macro", c(1,2,3,4,5)))
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
 
obj.sce1 <- as.SingleCellExperiment(obj)
sce.sling2 <- slingshot(obj.sce1, clusterLabels='ident', start.clus="macro3", reducedDim='PCA')
pseudo.paths <- slingPseudotime(sce.sling2)
SlingshotDataSet(sce.sling2 )


obj.sce1 <- runUMAP(obj.sce1 , dimred="PCA")
reducedDim(sce.sling2, "UMAP") <- reducedDim(obj.sce1 , "UMAP")
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
# gg <- plotUMAP(sce.sling2, colour_by=I(sce.sling2$seurat_clusters))
# gg <- plotUMAP(sce.sling2, colour_by=I(sce.sling2$ident))
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)

for (path in embedded) {
     embedded <- data.frame(path$s[path$ord,])
     gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
 }
#pdf("slingshot_1.pdf")
png("slingshot_time2.png", width = 6, height = 6, units = "in", res = 1800, pointsize = 2 )
gg
dev.off()

gg <- plotUMAP(sce.sling2, colour_by=I(sce.sling2$ident))
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)

for (path in embedded) {
     embedded <- data.frame(path$s[path$ord,])
     gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
 }
png("slingshot_cluster2.png", width = 6, height = 6, units = "in", res = 1800, pointsize = 2 )
gg
dev.off()





