# 🧬 TCGA-BRCA Transcriptomics: An End-to-End Pipeline for Differential Expression & ML Feature Engineering

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-Materials-lightgrey?style=for-the-badge)
![Genomics](https://img.shields.io/badge/Domain-Bioinformatics-success?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-ML--Ready-blue?style=for-the-badge)

## 🎯 Overview
This repository contains a professional-grade R workflow designed for the analysis of **Breast Invasive Carcinoma (TCGA-BRCA)** data. The pipeline automates the journey from raw GDC (Genomic Data Commons) retrieval to high-resolution publication figures and datasets optimized for **Machine Learning (ML)** classification.

## 🚀 Key Features
* **📡 Automated Data Sourcing:** Direct API integration with GDC using `TCGAbiolinks`.
* **⚖️ Balanced Cohort Selection:** Intelligent sampling logic to ensure equal representation of **Primary Tumor** and **Solid Tissue Normal** samples.
* **🔬 Statistical Rigor:** Robust Differential Expression Analysis (DEA) powered by `DESeq2`.
* **🤖 ML-Ready Exports:** Generates feature-engineered VST (Variance Stabilizing Transformation) matrices with integrated labels.
* **📊 Publication Graphics:** Automated 600 DPI plots (PCA, Volcano, Heatmap, MA, and Dispersion).

## 🛠️ Tech Stack & Dependencies
### **Environment**
* **R / RStudio** (Version 4.0 or higher recommended)

### **Bioconductor Packages**
* `TCGAbiolinks`: GDC Data Interface
* `DESeq2`: Differential Expression Analysis
* `SummarizedExperiment`: Genomic Container Management
* `maftools`: Mutation & Genomic Analysis

### **CRAN Packages**
* `tidyverse`: Data Wrangling
* `pheatmap`: Complex Heatmaps
* `ggrepel`: Advanced Labeling
* `openxlsx`: Excel Reporting

## 📋 Workflow Summary
1.  **Query:** Filter for RNA-Seq `STAR - Counts` and `open` access metadata.
2.  **Sampling:** Randomly select balanced groups (e.g., 10 Tumor vs 10 Normal).
3.  **Analysis:** Filter low-count genes and execute the `DESeq()` pipeline.
4.  **Transformation:** Convert raw counts to VST for downstream ML modeling.
5.  **Reporting:** Export statistically significant genes to `.xlsx`.

## 📂 Project Structure
```text
├── TCGABiolink.R             # Core Analysis Script
├── PCA_plot.png              # Sample Clustering (600 DPI)
├── Volcano_plot.png          # DEG Distribution
├── Heatmap_top20.png         # Top 20 Expressed Genes
├── 02 Output/
│   ├── DESeq2_results.xlsx   # Full DEG Statistics
│   └── VST_ML_dataset.xlsx   # ML Feature Matrix (Samples x Genes)
