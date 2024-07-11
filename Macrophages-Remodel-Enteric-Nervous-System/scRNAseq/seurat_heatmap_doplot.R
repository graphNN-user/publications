
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
#names(new.cluster.ids) <- levels(obj)
new.cluster.ids <- c(paste0("macro", seq(1,5) )) 
names(new.cluster.ids) <- levels(obj)

obj <- RenameIdents(obj, new.cluster.ids)

path  <- "/Users/liuy39/Research/H026/gene_list_from_Milena/"
files <- c("a_Grow_diff.csv","b_Scavanger_Rec.csv","c_Anti_infl_Axonal_guidance_Demyelination_Adhesion.csv","d_Microglia.csv","e_Neurotransmitters_Vasular_growth.csv")
files <- paste0(path, files)

my_features <- list()
for (f in files) {
  print(f)
  genes <- read.csv(f)
  my_features <- append(my_features, genes[[1]])
}
my_features <- unlist(my_features)
my_features <- rev(my_features)
my_features <- c(my_features, "Ccr2")

levels(obj) <- c("macro2", "macro3", "macro4", "macro1", "macro5")
# Clust2/Mac3 => Clust3/Mac4 => Clust0/Mac1 => Clust4/Mac5 => Clust1/Mac2
levels(obj) <- c("macro3", "macro4", "macro1", "macro5", "macro2")

pdf("heatmap_singleCell.pdf")
DoHeatmap(obj, features=my_features)
dev.off()

library(ggplot2)

png("markergene_dotplot4.png", width = 5, height = 9, units = "in", res = 1600, pointsize = 1 )
#DotPlot(object = obj, features = my_features, cols = c("blue","grey", "red")) + coord_flip()  + theme(axis.text.x = element_text(angle = 45, hjust=1 ))
DotPlot(object = obj, features = my_features ) + coord_flip()  + theme(axis.text.x = element_text(angle = 45, hjust=1 ))   + scale_colour_gradient2(low = "blue", mid = "grey", high = "red" )
dev.off()

 
