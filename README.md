# Transcriptomic Characterization of Colorectal Cancer

**Identification of dysregulated biological programs and tumor-microenvironment signatures from paired tumor/normal RNA-seq (GSE142279)**

![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3?logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-%E2%89%A5%203.18-1F8F4E)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![GEO](https://img.shields.io/badge/GEO-GSE142279-orange)

## Summary

This repository contains a modular, reproducible R pipeline that analyzes **paired tumor vs. normal RNA-seq** data from GEO dataset **GSE142279** (GRCh38, NCBI-generated raw counts) to:

1. Assess data quality and sample structure (PCA, sample clustering).
2. Identify differentially expressed genes (DEGs) with DESeq2 using a paired design.
3. Characterize dysregulated biological programs via functional enrichment (GO BP, KEGG, GSEA).
4. Estimate immune/stromal composition (xCell2 deconvolution) and test tumor-vs-normal differences with paired Wilcoxon signed-rank tests.

## Biological Objective

To define the transcriptional programs and tumor-microenvironment signatures that distinguish colorectal tumors from matched normal mucosa, accounting for patient-specific effects.

## Repository Structure

```
.
├── Data/
│   ├── GSE142279_raw_counts_GRCh38.p13_NCBI.tsv   # raw counts (download from GEO)
│   ├── Metadata.csv                               # sample annotation
│   └── *.rds                                      # intermediate objects (generated, git-ignored)
├── Scripts/
│   ├── 01_Dataset_Loading_Preprocessing.R         # loading, factor formatting, QC
│   ├── 02_DEGs.R                                  # DESeq2, shrinkage, annotation
│   ├── 03_Enrichment.R                            # GO BP, KEGG, GSEA
│   └── 04_Cell_Deconvolution.R                    # xCell2 + paired Wilcoxon
├── Figures/                                       # PCA, volcano, heatmaps, dotplots, boxplots
├── README.md
└── .gitignore
```

## Prerequisites

- R >= 4.3, Bioconductor >= 3.18

| Purpose | Packages |
|---|---|
| Differential expression | `DESeq2`, `apeglm` |
| Annotation | `org.Hs.eg.db`, `AnnotationDbi` |
| Enrichment | `clusterProfiler`, `enrichplot` |
| Deconvolution | `xCell2` |
| Visualization | `ggplot2`, `pheatmap`, `EnhancedVolcano`, `ggpubr` |
| Wrangling | `tidyverse` |

```r
install.packages(c("tidyverse", "BiocManager", "remotes", "pheatmap", "ggpubr"))
BiocManager::install(c("DESeq2", "apeglm", "org.Hs.eg.db", "AnnotationDbi",
                       "clusterProfiler", "enrichplot", "EnhancedVolcano"))
remotes::install_github("dviraran/xCell2")
```

Tip: use `renv::init()` and commit `renv.lock` to pin package versions.

## Running the Pipeline

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
```

Download the count matrix from GEO into `Data/`, then run the scripts **in order** from the repo root:

```r
source("Scripts/01_Dataset_Loading_Preprocessing.R")  # load counts + metadata, factors, QC
source("Scripts/02_DEGs.R")                           # DESeq2, shrinkage, annotation
source("Scripts/03_Enrichment.R")                     # GO BP, KEGG, GSEA
source("Scripts/04_Cell_Deconvolution.R")             # xCell2 + paired Wilcoxon
```

Or from the terminal:

```bash
for s in Scripts/0{1,2,3,4}_*.R; do Rscript "$s" || break; done
```

| Step | Inputs | Outputs |
|---|---|---|
| 01 | raw counts, `Metadata.csv` | cleaned `.rds` objects, PCA/QC plots |
| 02 | `.rds` from 01 | DEG tables, volcano plot, heatmap |
| 03 | DEG results | GO/KEGG dotplots, GSEA plots and tables |
| 04 | normalized expression, metadata | xCell2 scores, tumor-vs-normal boxplots |

## Outputs

- **QC:** PCA and sample clustering showing tumor/normal separation.
- **DEGs:** annotated table (shrunken log2FC, adjusted p-values, symbols), volcano plot, top-DEG heatmap.
- **Enrichment:** GO BP and KEGG dotplots; GSEA enrichment plots.
- **Microenvironment:** xCell2 cell-type scores and boxplots with paired Wilcoxon p-values.

![Volcano plot](Figures/3.Volcano%20Plot.png)
![Cell-type abundance](Figures/9.Cell-Type%20Abundance.png)

## Data Availability

Data come from NCBI GEO accession [GSE142279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142279). The raw count matrix (GRCh38.p13, NCBI-generated) is included in `Data/`, along with `Metadata.csv`. Intermediate `.rds` files are not included (see `.gitignore`) and are regenerated when you run the pipeline.

## Citation

If you use this code, please cite the dataset and the tools used:

- **Dataset:** GEO accession GSE142279 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142279)
- **DESeq2:** Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550 (2014).
- **clusterProfiler:** Wu T, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2(3), 100141 (2021).
- **xCell2:** see https://github.com/dviraran/xCell2 for the recommended citation.

## License

Released under the [MIT License](LICENSE).
