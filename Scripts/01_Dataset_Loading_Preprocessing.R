# ==============================================================================
# Script 1: Data Loading & Preprocessing
# ==============================================================================

# 1. Load Data
counts <- as.matrix(read.delim("Data/GSE142279_raw_counts_GRCh38.p13_NCBI.tsv", header = T, row.names = 1))
meta <- read.csv("Data/Metadata.csv", header = T, row.names = 1)

dim(counts)
dim(meta)

# 2. Check sample alignment

all(rownames(meta) == colnames(counts))

# 3. Format factor variables

meta$condition <- factor(meta$condition, levels = c("Tumor", "Normal"))

meta$patient_id <- factor(meta$patient_id)

class(meta$patient_id)
class(meta$condition)

# 4. Save processed artifacts for downstream scripts
saveRDS(counts, file = "counts_processed.rds")
saveRDS(meta, file = "meta_processed.rds")
