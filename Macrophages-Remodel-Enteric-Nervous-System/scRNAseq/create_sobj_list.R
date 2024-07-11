
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
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")
#  filtered_seurat <- subset(x =seurat_obj, subset= (nCount_RNA >=500 ) & (nFeature_RNA >= 250) & (nFeature_RNA < 4000) & (mitoRatio < 20))
#  print(seurat_obj)
  sample_list <- append(sample_list,sample1 )
  sobj_list <- append(sobj_list, seurat_obj) 
}

saveRDS(sobj_list, "seurat_obj_list.rds")


