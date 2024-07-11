
library(Seurat)
#library("scCustomize")
library(dplyr)
library(ggplot2)

DoubleMAD <- function(x){
   x         <- x[!is.na(x)]
   m         <- median(x)
   abs.dev   <- abs(x - m)
   left.mad  <- median(abs.dev[x<=m])
   right.mad <- median(abs.dev[x>=m])
   return(c(left.mad, right.mad))
}

DoubleMADsFromMedian <- function(x){
   two.sided.mad <- DoubleMAD(x)
   m <- median(x, na.rm=TRUE)
   x.mad <- rep(two.sided.mad[1], length(x))
   x.mad[x > m] <- two.sided.mad[2]
   mad.distance <- abs(x - m) / x.mad
   mad.distance[x==m] <- 0
   return(mad.distance)
}

sobj_list <- readRDS("seurat_obj_list.rds")
sobj_list2 <- list()
sample_list   <- list()
path <- "/Users/liuy39/Research/H026/scRNAseq/seurat_QC/doublet_consensus"

for (i in 1:length(sobj_list )) {

  s_obj <- sobj_list[[i]]
  metadata <- s_obj@meta.data
  orig.ident <- as.character(metadata[1,"orig.ident"])
  db_file = paste0(path,"/metadata_doublet_consensus_",orig.ident,".csv")
  db_meta   <- read.csv(db_file, row.names=1)
  db_meta2  <- db_meta[,"consensus", drop=F]
  metadata2 <- merge(metadata, db_meta2, by=0)
  rownames(metadata2) <- metadata2$Row.names
  metadata2 <- metadata2[,-1]
  s_obj@meta.data <- metadata2
  sobj_list2 <- append(sobj_list2, s_obj)
  sample_list <- append(sample_list,orig.ident)
}

merged_seurat <- merge(x = sobj_list2[[1]], y = sobj_list2[2:length(sobj_list2)], add.cell.ids = as.vector(sample_list))

metadata <- merged_seurat@meta.data 
nGene    <- metadata$nFeature_RNA
outliers <- nGene[DoubleMADsFromMedian(nGene ) > 2] 
outliers[order(outliers)][[1]]
# nFeature_RNA < 4487
# mitoratio < .20
merged_seurat
#84211

filtered_seurat <- subset(merged_seurat, subset= (mitoRatio < 20))                                    # 76892
filtered_seurat <- subset(filtered_seurat, subset= (nFeature_RNA < 4487 ) & (nFeature_RNA > 300))     # 68512
filtered_seurat <- subset(filtered_seurat, subset = (consensus == "singlet_consesus") )               # 64797

saveRDS(filtered_seurat, file="merged_seurat_doublet_removed_mito20_nGene4487.rds" )


