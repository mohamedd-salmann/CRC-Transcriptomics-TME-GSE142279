# ==============================================================================
# Script 2: DESeq2 Object Construction, QC & DEG Analysis
# ==============================================================================
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)
library(dplyr)

# 1. Load preprocessed data

counts <- readRDS("counts_processed.rds")
meta <- readRDS("meta_processed.rds")

# 2. Build DESeqDataSet and filter low counts

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ patient_id + condition
)

dim(dds)

keep <- rowSums(counts(dds) >= 10) >= 3
dds<- dds[keep, ]

dim(dds)

# 3. Fit DESeq2 model
dds <- DESeq(dds)
#########################################################

# 4. QC: PCA & Sample Distance Heatmap
vsd <- vst(dds, blind = F)
plotPCA(vsd, intgroup= "condition")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap::pheatmap(sampleDistMatrix, annotation_col = meta["condition"])
##########################################################

# 5. Differential Expression

res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

res <- res[order(res$padj), ]

res <- as.data.frame(res)

head(res)
#########################################################


# 6. Annotate with Gene Symbols
annots <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(res),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

res <- merge(res, annots, by.x = 0, by.y = "ENTREZID")
rownames(res) <- res[, 1]
res <- res[, -1]
res <- res[order(res$padj), ]
head(res)
######################################################

# 7. Subset significant genes

sig_genes <- subset(res, padj < 0.05)
Genes <- dplyr::select(sig_genes, "SYMBOL") %>% na.omit()
up_genes <- subset(sig_genes, log2FoldChange > 0)
down_genes <- subset(sig_genes, log2FoldChange < 0)

# Save DE results and fitted dds object
saveRDS(dds, file = "dds.rds")
saveRDS(res, file = "DEGs.rds")
write.csv(res, file = "DEGs.csv")
write.csv(Genes, file = "Sig_Genes.csv")
#####################################################

# 8. Visualizations: Volcano & Top 50 Heatmap

EnhancedVolcano(
  res, 
  lab= res$SYMBOL, 
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 10e-15,
  FCcutoff = 1,
  pointSize = 3,
  labSize = 2,
  colAlpha = 1
)

#Heatmap 

top50 <- head(sig_genes[order(sig_genes$padj), ], 50)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

top50_mat <- vsd_mat[rownames(top50), ]

gene_symbols <- top50$SYMBOL
rownames(top50_mat) <- make.unique(gene_symbols)

annotation_col <- meta[, c("patient_id", "condition")]
annotation_col <- annotation_col[colnames(top50_mat), , drop = FALSE]

pheatmap(
  top50_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8
)
