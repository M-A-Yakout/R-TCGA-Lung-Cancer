# R-TCGA-Lung-Cancer: Multi-Omics Differential Expression Analysis Pipeline

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.15+-green.svg)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub](https://img.shields.io/badge/GitHub-M--A--Yakout-181717?logo=github)](https://github.com/M-A-Yakout)

> A comprehensive R-based pipeline for differential gene expression analysis integrating TCGA and GEO datasets with pathway enrichment and visualization

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Datasets](#datasets)
- [Results Interpretation](#results-interpretation)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## 🔬 Overview

This pipeline performs comprehensive differential expression analysis across multiple datasets, integrating TCGA cancer genomics data with GEO viral infection studies. It implements state-of-the-art bioinformatics methods including DESeq2 for RNA-seq analysis, limma for microarray data, and clusterProfiler for functional enrichment.

### Key Analyses

1. **TCGA-LUAD Analysis** - Tumor vs Normal tissue comparison in lung adenocarcinoma
2. **GEO GSE188427** - RSV infection transcriptomic profiling
3. **GEO GSE197364** - HRSV infection in HEp-2 cells
4. **Pathway Enrichment** - Gene Ontology biological processes
5. **Visualization** - Volcano plots, heatmaps, and enrichment plots

## ✨ Features

- 🧬 **Multi-dataset integration**: TCGA and GEO platforms
- 📊 **Robust statistical methods**: DESeq2 with apeglm shrinkage for RNA-seq, limma for microarray
- 🎨 **Publication-ready visualizations**: Volcano plots, heatmaps, GO dotplots
- 🔍 **Pathway enrichment**: Gene Ontology biological processes analysis
- 📁 **Automated output management**: Organized results directory structure
- 🛡️ **Error handling**: Comprehensive try-catch blocks for robust execution
- 📝 **Detailed logging**: Step-by-step progress messages

## 🔧 Prerequisites

### System Requirements

- **R version**: ≥ 4.0.0
- **RAM**: ≥ 8GB recommended (16GB for large TCGA datasets)
- **Storage**: ≥ 10GB free space for data downloads
- **OS**: Linux, macOS, or Windows with Rtools

### Required R Packages

#### Bioconductor Packages
```r
TCGAbiolinks, GEOquery, DESeq2, apeglm, limma, 
SummarizedExperiment, clusterProfiler, org.Hs.eg.db,
pheatmap, ComplexHeatmap
```

#### CRAN Packages
```r
dplyr, ggplot2, BiocManager, remotes
```

## 📥 Installation

### Quick Start

```bash
# Clone the repository
git clone https://github.com/M-A-Yakout/R-TCGA-Lung-Cancer.git
cd R-TCGA-Lung-Cancer

# Run the setup script
Rscript step_0.R
```

### Manual Installation

```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "TCGAbiolinks", "GEOquery", "DESeq2", "apeglm", 
  "limma", "SummarizedExperiment", "clusterProfiler", 
  "org.Hs.eg.db", "pheatmap", "ComplexHeatmap"
))

# Install CRAN packages
install.packages(c("dplyr", "ggplot2"))
```

## 🚀 Usage

### Basic Usage

```bash
# Run the complete pipeline
Rscript step_0.R  # Setup and package installation
Rscript step_1.R  # TCGA-LUAD analysis
Rscript step_2.R  # GEO GSE188427 analysis
Rscript step_3.R  # GEO GSE197364 analysis
Rscript step_4.R  # Pathway enrichment
Rscript step_5.R  # Heatmap generation
Rscript step_results.R  # Summary report
```

### R Session Usage

```r
# Source all scripts sequentially
source("step_0.R")
source("step_1.R")
source("step_2.R")
source("step_3.R")
source("step_4.R")
source("step_5.R")
source("step_results.R")
```

### Running Specific Analyses

```r
# Run only TCGA analysis
source("step_0.R")
source("step_1.R")

# Run only GEO datasets
source("step_0.R")
source("step_2.R")
source("step_3.R")
```

## 📊 Pipeline Steps

### Step 0: Environment Setup
- Installs required packages
- Loads necessary libraries
- Creates results directory

### Step 1: TCGA-LUAD Analysis
- Downloads TCGA lung adenocarcinoma RNA-seq data (STAR counts)
- Performs differential expression analysis (Tumor vs Normal)
- Applies apeglm log fold change shrinkage
- Generates volcano plot
- Filters significant DEGs (padj < 0.05, |log2FC| ≥ 1)

### Step 2: GEO GSE188427 Analysis
- Retrieves microarray data from GEO
- Processes RSV infection samples
- Performs limma differential expression analysis
- Creates volcano plot visualization

### Step 3: GEO GSE197364 Analysis
- Downloads HRSV infection dataset
- Compares HRSV-infected vs Mock control
- Generates differential expression results
- Produces volcano plot

### Step 4: Pathway Enrichment
- Performs Gene Ontology enrichment analysis
- Identifies significant biological processes
- Creates dotplot and barplot visualizations

### Step 5: Heatmap Generation
- Visualizes top 50 DEGs from TCGA analysis
- Clusters samples and genes hierarchically
- Annotates sample groups

## 📁 Output Files

All results are saved in the `Results/` directory:

### CSV Files
- `TCGA_LUAD_sigDEGs.csv` - Significant DEGs from TCGA analysis
- `GSE188427_sigDEGs.csv` - Significant DEGs from GSE188427
- `GSE197364_sigDEGs.csv` - Significant DEGs from GSE197364
- `TCGA_GO_enrichment.csv` - Gene Ontology enrichment results

### PDF Visualizations
- `TCGA_LUAD_volcano.pdf` - Volcano plot for TCGA data
- `GSE188427_volcano.pdf` - Volcano plot for GSE188427
- `GSE197364_volcano.pdf` - Volcano plot for GSE197364
- `TCGA_GO_dotplot.pdf` - GO enrichment dotplot
- `TCGA_GO_barplot.pdf` - GO enrichment barplot
- `TCGA_LUAD_heatmap_top50.pdf` - Heatmap of top 50 DEGs

### RDS Objects
- `TCGA_LUAD_results.rds` - DESeq2 results object
- `TCGA_LUAD_dds.rds` - DESeqDataSet object

## 🗂️ Datasets

### TCGA-LUAD
- **Source**: The Cancer Genome Atlas
- **Project**: TCGA-LUAD (Lung Adenocarcinoma)
- **Data Type**: RNA-Seq (STAR - Counts)
- **Samples**: Primary Tumor and Solid Tissue Normal
- **Platform**: Illumina HiSeq

### GSE188427
- **Source**: Gene Expression Omnibus
- **Study**: RSV infection transcriptomic analysis
- **Platform**: Microarray
- **Comparison**: RSV-infected vs Healthy samples

### GSE197364
- **Source**: Gene Expression Omnibus
- **Study**: HRSV infection in HEp-2 cells
- **Platform**: Microarray
- **Comparison**: HRSV-infected vs Mock control

## 📖 Results Interpretation

### Differential Expression Criteria
- **Adjusted p-value**: < 0.05 (FDR correction)
- **Log fold change**: |log2FC| ≥ 1 (2-fold change)
- **Normalization**: DESeq2 size factors (TCGA), quantile normalization (GEO)

### Volcano Plot Interpretation
- **Red points**: Significant DEGs
- **Grey points**: Non-significant genes
- **Blue dashed lines**: Significance thresholds (padj = 0.05, |log2FC| = 1)

### Heatmap Features
- **Color scale**: Blue (low) → White (median) → Red (high)
- **Clustering**: Hierarchical clustering of both genes and samples
- **Scaling**: Row-wise z-score normalization

## 🔍 Troubleshooting

### Common Issues

**Issue**: Package installation fails
```r
# Solution: Update Bioconductor
BiocManager::install(version = "3.15", ask = FALSE, update = TRUE)
```

**Issue**: TCGA download timeout
```r
# Solution: Reduce chunk size
GDCdownload(query, method = "api", files.per.chunk = 5)
```

**Issue**: Memory error with large datasets
```r
# Solution: Increase R memory limit (Windows)
memory.limit(size = 16000)
```

**Issue**: GEO data retrieval fails
```r
# Solution: Use alternative GEO mirror
options(GEOquery.inmemory.gpl = FALSE)
```

## 📚 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{r_tcga_lung_cancer,
  author = {Mohamed Mostafa},
  title = {R-TCGA-Lung-Cancer: Multi-Omics Differential Expression Analysis Pipeline},
  year = {2024},
  url = {https://github.com/M-A-Yakout/R-TCGA-Lung-Cancer},
  version = {1.0.0}
}
```

### Key Publications

- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014). Genome Biology.
- **limma**: Ritchie, M.E., et al. (2015). Nucleic Acids Research.
- **clusterProfiler**: Yu, G., et al. (2012). OMICS.
- **TCGAbiolinks**: Colaprico, A., et al. (2016). Nucleic Acids Research.

## 🤝 Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Guidelines
- Follow R coding style conventions
- Add appropriate error handling
- Update documentation for new features
- Include test cases where applicable

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact

**Project Maintainer**: Mohamed Mostafa  
**Email**: mohamed.ashraf.y.s.m@gmail.com  
**GitHub**: [@M-A-Yakout](https://github.com/M-A-Yakout)  
**Issues**: [GitHub Issues](https://github.com/M-A-Yakout/R-TCGA-Lung-Cancer/issues)

---

## 🙏 Acknowledgments

- TCGA Research Network for providing cancer genomics data
- GEO database contributors and submitters
- Bioconductor community for excellent tools and documentation
- R Core Team for the R statistical environment

---

**Last Updated**: December 15, 2024  
**Pipeline Version**: 1.0.0  
**Compatibility**: R ≥ 4.0.0, Bioconductor ≥ 3.15

<p align="center">Made with ❤️ for the bioinformatics community</p>
