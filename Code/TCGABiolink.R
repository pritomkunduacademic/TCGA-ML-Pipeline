# 01 Install Needed Packages
# 1.1 Biocmanager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 1.2 Maftool
BiocManager::install("maftools")

# 1.3 Summarized Experiment
BiocManager::install("SummarizedExperiment")

# 1.4 pheatmap
install.packages("pheatmap")
install.packages("RColorBrewer")

# 2 Loading needed Library
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(SummarizedExperiment)
library(pheatmap)
library(DESeq2)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)

# 3 Getting The list of Projects
df <- getGDCprojects()

# 3.1 Geting The Project Summary

 getProjectSummary("TCGA-BRCA")

# 4 Building Up a GDC Query
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling")
# 4.1 Getting the query output
out_query <- getResults(query)

# 5 Define the total number of samples here

n_samples <- 20

if (n_samples %% 2 != 0) {
  stop("n_samples must be even for balanced classes")
}
n_each <- n_samples/2

# 6 Building GDC Query For RNA Seq Meta Data Retrial
query_meta <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  access = "open")

meta <- getResults(query_meta)

# 7 Extract Tumor and Normal Samples
tumor_samples <- unique(meta$cases[meta$sample_type == "Primary Tumor"])
normal_samples <- unique(meta$cases[meta$sample_type == "Solid Tissue Normal"])

# 8 Handle edge case (not enough normals)
if (length(normal_samples) < n_each) {
  warning("Not enough normal samples available. Adjusting sample size.")
  n_each <- length(normal_samples)
}
# 9 Random balanced sampling (reproducible)
set.seed(42)

selected_tumor <- sample(tumor_samples, n_each)
selected_normal <- sample(normal_samples, n_each)

selected_samples <- c(selected_tumor, selected_normal)

# 10 Final Query
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  access = "open",
                  barcode = selected_samples)

# 11 Download Data
GDCdownload(query, method = "api", files.per.chunk = 20)

# 12 Prepare Data (SummarizedExperiment)
data_se <- GDCprepare(query)

# 13 Inspect Data
data_se

# 14 # Count matrix
count_matrix <- assay(data_se)

# 14  Metadata
meta <- colData(data_se)

# 15  Check dimensions
dim(counts)
head(meta$sample_type)

# 16 Class Check
class(count_matrix)
dim(count_matrix)

# 17 Filter the data 
keep <- rowSums(count_matrix) >= 10
counts_filtered <- count_matrix[keep, ]

# 18 Prepare metadata
meta$condition <- factor(meta$sample_type)

# 19 Keep only Tumor vs Normal if needed
keep_samples <- meta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")

counts_filtered <- counts_filtered[, keep_samples]
meta <- meta[keep_samples, ]

# 20 Create DESeq2 Object
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = meta,
  design = ~ condition
)

# 21 Run Deseq2
dds <- DESeq(dds)

# 22 # Results
res <- results(dds)

# 23 View Data 
head(res[order(res$pvalue), ])

# 24 Transform counts (log scale)
vsd <- vst(dds, blind = FALSE)

ml_data <- assay(vsd)

# 25 Create labels
labels <- meta$condition

# ------------------------------------------------------------
# 26 Creating Publication Plots
# ------------------------------------------------------------
# 26.1 Global theme

theme_nature <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# 26.2 PCA Plot (VST-based)

png("PCA_plot.png", width = 10, height = 6, units = "in", res = 600)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("PCA of Samples") +
  theme_nature

dev.off()

# 26.3 Volcano Plot
library(ggrepel)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- na.omit(res_df)

res_df$group <- "NS"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# Top genes for labeling
top_genes <- head(res_df[order(res_df$padj), ], 10)

png("Volcano_plot.png", width = 10, height = 6, units = "in", res = 600)

ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = group), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "#0072B2", "NS" = "grey70", "Up" = "#D55E00")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot") +
  theme_nature

dev.off()

# 26.4 Top 20 Genes Heatmap

# Remove NA adjusted p-values
res_clean <- res[!is.na(res$padj), ]

#  Select top 20 DEGs safely
top20 <- rownames(res_clean)[order(res_clean$padj)][1:20]

#  Extract expression matrix
mat <- assay(vsd)[top20, ]

#  Remove any remaining NA rows (safety)
mat <- mat[complete.cases(mat), ]

#  Z-score scaling (very important for heatmap visualization)
mat <- t(scale(t(mat)))

#  Create annotation (sample labels)
annotation_col <- data.frame(condition = colData(dds)$condition)
rownames(annotation_col) <- colnames(mat)

#  Define color palette (publication standard)
heat_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

#  Generate and SAVE heatmap safely
pheatmap(mat,
         annotation_col = annotation_col,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         border_color = NA,
         filename = "Heatmap_top20.png",
         width = 10,
         height = 6)
# 26.5 MA Plot

png("MA_plot.png", width = 10, height = 6, units = "in", res = 600)

plotMA(res, ylim = c(-4, 4), main = "MA Plot")

dev.off()

# 26.6 Dispersion Plot

png("Dispersion_plot.png", width = 10, height = 6, units = "in", res = 600)

plotDispEsts(dds)

dev.off()

# ------------------------------------------------------------
# 27 File Export 
# ------------------------------------------------------------

# 27.1 Saving Deseq2 Results (Convert DESeq2 results to data frame)

res_df <- as.data.frame(res)

# Add gene column
res_df$gene <- rownames(res_df)

# Order by p-value (recommended for reporting)
res_df <- res_df[order(res_df$pvalue), ]

# Reorder columns (gene first)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

# Write to Excel
write.xlsx(res_df, file = "02 Output/02_DESeq2_results_BRCA_20.xlsx")

# 27.2 Saving the VST ML ready File
# Extract ML matrix (genes x samples → samples x genes)
ml_data <- t(assay(vsd))

#  Add labels
labels <- meta$condition

# Create final ML dataframe
ml_df <- as.data.frame(ml_data)
ml_df$label <- labels

# Add sample IDs as a column
ml_df$sample <- rownames(ml_df)

# Reorder columns (sample + label first)
ml_df <- ml_df[, c("sample", "label", setdiff(colnames(ml_df), c("sample", "label")))]


#  Save to Excel
write.xlsx(ml_df, file = "02 Output/BRCA_VST_ML_dataset.xlsx")