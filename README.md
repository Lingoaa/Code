# Full Analysis Pipelines for Colorectal Cancer (CRC) Evaluation
This repository contains self-written customized R analysis scripts integrating published public R packages for single-cell RNA-seq processing and CellChat-based cell-cell communication analysis described in the manuscript.

## 1. System Requirements
- **Operating System:** Windows 10/11, macOS, or Linux.
- **Software Dependencies:** R (version 4.2.0 or higher) and RStudio.
- **Required public R Packages:** Seurat, CellChat, tidyverse, patchwork, ggplot2, pheatmap, ggpubr, reshape2, BiocManager, preprocessCore.
- No non-standard dedicated hardware is required; systematic cross-version testing of the pipeline was not implemented.

## 2. Installation guide
1. Preinstall R (≥4.2.0) and RStudio locally.
2. Run installation commands embedded at the start of scripts to automatically install dependent public packages.
3. Full dependency installation costs around 5–10 min on a regular desktop PC.

## 3. Instructions for Custom Data Use
1. Clone or download all .R scripts from the repository.
2. Modify `data_dir` file path (around line 29) inside scripts to your local expression matrix folder.
3. Execute codes sequentially.

## 4. Reproduction Instructions
All figures and quantitative results in this paper can be reproduced by running provided scripts with raw/processed sequencing data listed in Data Availability statement of manuscript.
