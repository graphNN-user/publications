
library(Seurat)
library(dplyr)
library(ggplot2)

path <- "/Users/liuy39/Research/H026/scRNAseq/"
meta <- read.table(paste0(path,"metadata/scRNAseq_meta_summary.txt"), header=T)
dir_list <- list.dirs(paste0(path,"cellranger_out"), recursive = FALSE)
dir_list <- paste0(dir_list,"/filtered_feature_bc_matrix.h5")
all(file.exists(dir_list))

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

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
  seurat_obj@meta.data <- metadata
  print(seurat_obj)
  sample_list <- append(sample_list,sample1 )
  sobj_list <- append(sobj_list, seurat_obj)
}

merged_seurat <- merge(x = sobj_list[[1]], y = sobj_list[2:length(sobj_list)], add.cell.ids = as.vector(sample_list))
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
#merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

VlnPlot(merged_seurat, features = c("nCount_RNA", "nFeature_RNA", "mitoRatio" )) 

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
metadata <- rename(metadata, sample = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA)
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

ggplot(metadata, aes(x=sample, fill=sample )) + geom_bar() + ggtitle("Number of Cells")

ggplot(metadata, aes(color=sample, x=nUMI )) + geom_density(alpha = 0.2) + scale_x_log10() +  ylab("Cell density") + ggtitle("density of nUMI")
ggplot(metadata, aes(color=sample, x=nGene )) + geom_density(alpha = 0.2) + scale_x_log10() +  ylab("Cell density") + ggtitle("density of nGene")
ggplot(metadata, aes(x=sample, y=log10(nGene), fill=sample )) + geom_boxplot() + ggtitle("boxplot of nGenes")
ggplot(metadata, aes(x=nUMI, y=nGene, color=mitoRatio)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + 
 stat_smooth(method=lm) + scale_x_log10() +  scale_y_log10() + facet_wrap(~sample)

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

nGene <- metadata$nGene
outliers <- nGene[DoubleMADsFromMedian(nGene ) > 3]
outliers[order(outliers)]
