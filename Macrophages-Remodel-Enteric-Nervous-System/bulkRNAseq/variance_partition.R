
library(variancePartition)
library(dplyr)
library(DESeq2)
library(tximeta)
suppressPackageStartupMessages(library(SummarizedExperiment))
library(org.Mm.eg.db)


info <- read.table("coldata_rin_exon.txt", header=T, sep="\t")
info$files <- as.character(info$files)
all(file.exists(info$files))
se <- tximeta(info, dropInfReps=TRUE)
assayNames(se)
rowRanges(se)

gse <- summarizeToGene(se)
rowRanges(gse)
x <- GRanges("chr1", IRanges(10e6, 11e6))
gse[gse %over% x, ]
gse_dss <- gse[, gse$condition == "DSS"]
gse <- addIds(gse, column="SYMBOL")
columns(org.Mm.eg.db)
assayNames(gse)
cs <- colSums(assay(gse, "counts"))
hist(cs/1e6, col="grey", border="white",main="", xlab="column sums (per million)")
dds <- DESeqDataSet(gse, design=~condition + cell_status)
keep <- rowSums(counts(dds) >= 5) >= 3
table(keep)

dds <- dds[keep,]
dds <- DESeq(dds)

counts_normalized <- counts(dds, normalized=T)
class(counts_normalized)
colnames(counts_normalized)
counts_normalized2 <- counts_normalized[, info$names]

info$cell_status2  <- as.factor(info$cell_status2)
info$condition <- as.factor(info$condition)

info1 <- dplyr::select(info, names, condition,cell_status2, RINe, exonmapping)
rownames(info1) <- info1$names
form <- ~ RINe + exonmapping  + (1|condition) + (1|cell_status2)
varPart <- fitExtractVarPartModel(counts_normalized2, form, info1)
vp <- sortCols( varPart )
plotVarPart( vp )

pdf("varpart.pdf")
plotVarPart( vp )
dev.off()
