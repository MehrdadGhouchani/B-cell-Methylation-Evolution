# B-cell-Methylation-Evolution

# R pipeline for DNA methylation analysis of B cell neoplasms (426 samples).

[![R](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

| # | Script | Status | Analysis |
|---|--------|--------|----------|
| 01 | 01_qc_normalization.R | ⏳ | Quality control |
| 02 | 02_heatmap_variability.R | ✅ | Variable CpGs |
| 03 | 03_pca_tsne.R | ⏳ | Dimensionality reduction |
| 04 | 04_epiCMIT.R | ⏳ | Epigenetic clocks |
| 05 | 05_wilcoxon_cll.R | ⏳ | CLL subtype tests |
| 06 | 06_correlation.R | ⏳ | Methylation correlations |
| 07 | 07_differential.R | ⏳ | DM analysis |
| 08 | 08_DMRcate.R | ✅ | Genome annotation |
| 09 | 09_kaplan_meier.R | ⏳ | Survival curves |

# Project
Decoding evolutionary dynamics of B cell neoplasms using DNA methylation profiling.

- 426 B-cell lymphoma samples
- Illumina MethylationEPIC v2

---
Mehrdad Ghouchani
