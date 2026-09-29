# ==============================================================================
# Script 4: Cell Type Deconvolution (xCell2)
# ==============================================================================
library(xCell2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(rstatix)

# 1. Load preprocessed counts and metadata
counts <- readRDS("counts_processed.rds")
meta <- readRDS("meta_processed.rds")

# 2. Map Entrez IDs to Gene Symbols for xCell2
xcell_expr <- counts
rownames(xcell_expr) <- mapIds(
  org.Hs.eg.db,
  keys = rownames(counts), 
  column = "SYMBOL", 
  keytype = "ENTREZID", 
  multiVals = "first"
)
head(xcell_expr)

#Remove NAs
xcell_expr <- xcell_expr[!is.na(rownames(xcell_expr)), ]
xcell_expr <- xcell_expr[!duplicated(rownames(xcell_expr)), ]

# 3. Run xCell2 deconvolution
data("BlueprintEncode.xCell2Ref")
xcell_res <- xCell2Analysis(mix = xcell_expr, xcell2object = BlueprintEncode.xCell2Ref)
################################################################################

# 4. Format and merge with metadata
xcell_df <- as.data.frame(t(xcell_res))
xcell_df$sample_id <- rownames(xcell_df)

merged_cell <- merge(xcell_df, meta, by.x = "sample_id", by.y = 0)

# Identify cell-type columns
cell_types <- setdiff(
  colnames(merged_cell), 
  c("sample_id", "sample_name", "patient_id", "condition", "tissue", "sizeFactor")
)

# Reshape into long format
cell_long <- merged_cell %>%
  dplyr::select(patient_id, condition, all_of(cell_types)) %>%
  pivot_longer(
    cols = all_of(cell_types), 
    names_to = "cell_type", 
    values_to = "abundance"
  )

################################################################################

# 5. Summary statistics and paired Wilcoxon tests
cell_summary <- cell_long %>%
  group_by(cell_type, condition) %>%
  summarize(
    mean = mean(abundance), 
    median = median(abundance), 
    .groups = "drop"
  )

cell_stats <- cell_long %>%
  group_by(cell_type) %>%
  wilcox_test(abundance ~ condition, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  arrange(p.adj)

head(cell_stats)

################################################################################

# 6. Visualization: Top 10 differentially abundant cell types
top_cells <- cell_stats %>%
  slice_head(n = 10) %>%
  pull(cell_type)

ggplot(cell_long %>% filter(cell_type %in% top_cells), aes(condition, abundance)) + 
  geom_boxplot() + 
  geom_point(aes(group = patient_id), position = position_jitter(width = 0.05)) + 
  facet_wrap(~cell_type, scales = "free_y") + 
  theme_bw() + 
  labs(
    title = "Normal vs Tumor Cell-Type Abundance", 
    x = NULL,
    y = "Estimated Abundance"
  )