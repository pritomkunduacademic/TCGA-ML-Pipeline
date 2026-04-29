# 🧬 TCGA-BRCA Transcriptomics: An End-to-End Pipeline for Differential Expression & ML Feature Engineering

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-Materials-lightgrey?style=for-the-badge)
![Genomics](https://img.shields.io/badge/Domain-Bioinformatics-success?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-ML--Ready-blue?style=for-the-badge)

## 🎯 Overview
This repository contains a professional-grade R workflow designed for the analysis of **Breast Invasive Carcinoma (TCGA-BRCA)** data. The pipeline automates the journey from raw GDC (Genomic Data Commons) retrieval to high-resolution publication figures and datasets optimized for **Machine Learning (ML)** classification. 

The workflow integrates advanced **Genomic Annotation** and **Biotype Stratification**, ensuring that statistical outputs are mapped to biological functions and filtered for high-confidence feature engineering.

## 🚀 Key Features
* **📡 Automated Data Sourcing:** Direct API integration with GDC using `TCGAbiolinks`.
* **⚖️ Balanced Cohort Selection:** Intelligent sampling logic to ensure equal representation of **Primary Tumor** and **Solid Tissue Normal** samples.
* **🔬 Statistical Rigor:** Robust Differential Expression Analysis (DEA) powered by `DESeq2`.
* **🏷️ Intelligent Annotation:** Automated mapping of Ensembl IDs to HGNC symbols using `EnsDb.Hsapiens.v86`.
* **🧬 Biotype Filtering:** Specialized exports for **Protein-Coding** vs. **Non-Coding** genes to refine biological focus.
* **🤖 ML-Ready Exports:** Generates feature-engineered VST (Variance Stabilizing Transformation) matrices with integrated labels.
* **📊 Publication Graphics:** Automated 600 DPI plots including PCA, Volcano, and Heatmaps.

## 🛠️ Tech Stack & Dependencies
### **Environment**
* **R / RStudio** (Version 4.0 or higher recommended)

### **Bioconductor Packages**
* `TCGAbiolinks`: GDC Data Interface
* `DESeq2`: Differential Expression Analysis
* `SummarizedExperiment`: Genomic Container Management
* `biomaRt` / `org.Hs.eg.db` / `EnsDb.Hsapiens.v86`: Gene Annotation & Mapping

### **CRAN Packages**
* `tidyverse` (dplyr, ggplot2): Data Wrangling & Visualization
* `pheatmap`: Complex Heatmaps
* `openxlsx`: Excel Reporting

## 📋 Workflow Summary
1.  **Query:** Filter for RNA-Seq `STAR - Counts` and `open` access metadata.
2.  **Sampling:** Randomly select balanced groups (e.g., 10 Tumor vs 10 Normal).
3.  **Analysis:** Filter low-count genes and execute the `DESeq()` pipeline.
4.  **Annotation:** Strip Ensembl version tags and map IDs to Gene Symbols.
5.  **Stratification:** Categorize genes by Biotype (Protein-Coding) and Expression Status (Over-expressed/Down-regulated).
6.  **Transformation:** Convert raw counts to VST for downstream ML modeling.

## 🧬 Academic & Research Insights
### **Expression Categorization Logic**
The pipeline applies a dual-threshold filter to ensure both statistical significance and biological magnitude:
* **Significance:** $padj < 0.05$ (False Discovery Rate control).
* **Magnitude:** $|log2FoldChange| > 1$ (Minimum 2-fold change).

### **Feature Engineering for ML**
By isolating `protein_coding` genes, the workflow reduces noise and dimensionality for predictive modeling. This ensures that ML classifiers focus on the proteome-contributing subset most relevant to pharmaceutical research and biomarker discovery.

## 📂 Project Structure
```text
├── TCGABiolink.R              # Core Analysis & Normalization
├── Gene_Mapping_Filtering.R   # Annotation & Biotype Stratification
├── PCA_plot.png               # Sample Clustering (600 DPI)
├── Volcano_plot.png           # DEG Distribution
├── 02 Output/                 # Processed Data
│   ├── DESeq2_results.xlsx    # Full DEG Statistics
│   ├── Annotated_Full_Results.csv # Results with Gene Symbols
│   ├── 1_Protein_Coding_Genes.csv # ML-ready Protein Coding Features
│   ├── 2_Protein_Non_Coding_Genes.csv # Regulatory RNA subset
│   └── 3_Over_expressed_Genes.csv # Potential Oncogenic Drivers
