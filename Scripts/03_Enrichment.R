# ==============================================================================
# Script 3: Functional Enrichment Analysis (GO, KEGG, GSEA)
# ==============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(cowplot)

# 1. Load DEG results
res <- readRDS("DEGs.rds")

sig_genes <- subset(res, padj < 0.05)
up_genes <- subset(sig_genes, log2FoldChange > 0)
down_genes <- subset(sig_genes, log2FoldChange < 0)

sig_genes$ENTREZID <- rownames(sig_genes)
up_genes$ENTREZID <- rownames(up_genes)
down_genes$ENTREZID <- rownames(down_genes)


# Remove NA Entrez IDs
up_entrez <- na.omit(up_genes$ENTREZID)
down_entrez <- na.omit(down_genes$ENTREZID)
#######################################################

# 2. GO Over-representation Analysis
ego_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_up)

ego_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_down)

up_plot <- dotplot(ego_up, showCategory = 10, font.size = 7.5, title = "GO BP - Upregulated")
down_plot <- dotplot(ego_down, showCategory = 10, font.size = 7.5, title = "GO BP - Downregulated")
plot_grid(up_plot, down_plot, ncol = 2)
#############################################################

# 3. KEGG Pathway Analysis

kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
head(kegg_up)

kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

head(kegg_down)

kegg_up_plot <- dotplot(kegg_up, showCategory = 10, title = "KEGG - Upregulated")
kegg_down_plot <- dotplot(kegg_down, showCategory = 10, title = "KEGG - Downregulated")

plot_grid(kegg_up_plot, kegg_down_plot, nrow = 2)
###########################################################



# 4. Gene Set Enrichment Analysis (GSEA)

geneList <- res$stat
names(geneList) <- rownames(res)
head(geneList)
sum(is.na(geneList))
sum(is.na(names(geneList)))

#If there are NAs, remove it
#geneList <- geneList[!is.na(names(geneList))]

geneList <- sort(geneList, decreasing = TRUE)


gsea_kegg <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

dotplot(gsea_kegg, showCategory = 15, font.siz= 10, title = "Enriched KEGG Pathways")

gseaplot2(
  gsea_kegg,
  geneSetID = 1,
  title = "KEGG: Ribosome (Downregulated in Tumor)"
)

# Inspect top GSEA results
head(gsea_kegg@result[, c("ID", "Description", "NES", "p.adjust")])
