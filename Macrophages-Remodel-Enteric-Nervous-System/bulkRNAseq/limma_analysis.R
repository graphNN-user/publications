library(dplyr)
library(edgeR)
library(ggplot2)

data_dir <- "/Users/liuy39/Research/H026/rnaseqc"
dirs     <- list.dirs(data_dir)
dirs     <- dirs[grepl("H26", dirs)]
samples  <- basename(dirs)
samples  <- samples[grepl("^H26", samples)]
files    <- paste0(dirs, "/", samples, "R1_sorted.bam.gene_reads.gct")
all(file.exists(files))

mylist <- list()

for (i in c(1:length(files ))) {
  exprM <- read.table(files[[i]],sep="\t",header=T,row.names=1,skip=2)
  exprM <- exprM[,-1, drop=F]
  colnames(exprM) <- samples[[i]]
  mylist[[i]] = exprM
}
df <- do.call("cbind",mylist) #combine all vectors into a matrix

info <- read.table("coldata_rin_exon2.txt", header=T, sep="\t")
info1<- dplyr::select(info, names, condition,cell_status2, RINe, exonmapping)

info2 <- info1[ (info1$condition == "DSS" & info1$cell_status2 == "C_pos_L_pos_M_neg") | (info1$condition == "Tap" & info1$cell_status2 == "C_neg_L_neg_M_neg"), ]
info2$combined <- paste0(info2$condition,"_", info2$cell_status2 )

rownames(info2) <- info2$names
df2   <- df[,info2$names]

d0 <- DGEList(df2)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d    <- d0[-drop,]

condition  <- factor(info2$condition, levels=c("Tap","DSS"))
cellstatus <- as.factor(info2$cell_status2)
exonmapping<- info2$exonmapping
cellstatus <- factor(cellstatus,levels=c("C_pos_L_pos_M_neg","C_pos_L_pos_M_pos","C_pos_L_neg_M_pos","C_neg_L_neg_M_pos"))
combined <- factor(info2$combined, levels=c("Tap_C_neg_L_neg_M_neg", "DSS_C_pos_L_pos_M_neg"))

design <- model.matrix(~ 0 + combined + exonmapping)
y <- voom(d, design, plot=T)

fit <- lmFit(y, design)
head(coef(fit))

# boxplot
boxplot_fun <- function(targetgene) {
  expr      <- d$counts
#  gene1     <- "ENSMUSG00000064368.1"
  gene1     <- targetgene
  gene1_df  <- expr[rownames(expr) == gene1,,drop=F]
  gene1_dft <- t(gene1_df)
  colnames(gene1_dft) <- "genename"
  gene1_dft2 <- merge(gene1_dft, info2, by=0)
  ggplot(gene1_dft2, aes(x=cell_status2,y=genename,fill=combined)) + geom_boxplot() + ggtitle(targetgene) + theme(axis.text.x = element_text(angle = 45))
}

contr <- makeContrasts(combinedTap_C_neg_L_neg_M_neg  - combinedDSS_C_pos_L_pos_M_neg, levels=colnames(design))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n= Inf)
head(top.table, 10)

p1 <- boxplot_fun("ENSMUSG00000040152.9")
p1 <- boxplot_fun("ENSMUSG00000006360.12")
genemap <- read.table("H26-02_S9_R1_sorted.bam.gene_reads.gct", skip=2, header=T,sep="\t", row.names=1)
genemap <- genemap[,-2,drop=F]
top.table2 <- merge(genemap, top.table, by=0)
top.table3 <- arrange(top.table2, adj.P.Val)
write.csv(top.table3, file="limma_water0_vs_Subset1_C_pos_L_pos_M_neg.csv",row.names=F)

top.table3_up <- top.table3[ top.table3$adj.P.Val < 0.05 & top.table3$logFC > 1.5,]
write.csv(top.table3_up, file="limma_water0_vs_Subset1_C_pos_L_pos_M_neg_up.csv",row.names=F)

top.table3_down <- top.table3[ top.table3$adj.P.Val < 0.05 & top.table3$logFC < -1.5,]
write.csv(top.table3_down, file="limma_water0_vs_Subset1_C_pos_L_pos_M_neg_down.csv",row.names=F)

