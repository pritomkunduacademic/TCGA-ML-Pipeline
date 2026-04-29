# 01 Installation of needed package
BiocManager::install("biomaRt")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnsDb.Hsapiens.v86")

# 02 Loading Needed Libraries
library(dplyr)
library(openxlsx)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

# 03 Reading The Data
df <- read.xlsx("02 Output/DESeq2_results_BRCA_20.xlsx") 

# 04 Pre-process Gene IDs
df$ensembl_id_clean <- gsub("\\..*", "", df$gene)

# 05 Retrieve the mapping
gene_annot <- genes(EnsDb.Hsapiens.v86, 
                    filter = GeneIdFilter(df$ensembl_id_clean),
                    return.type = "data.frame") %>%
  dplyr::select(gene_id, symbol, gene_biotype)

# 06 Add Symbols back to Dataframe
df <- df %>%
  dplyr::select(-any_of(c("gene_symbol", "symbol"))) %>% 
  dplyr::left_join(gene_annot, by = c("ensembl_id_clean" = "gene_id")) %>%
  dplyr::rename(gene_symbol = symbol)


# 07 Professional Formatting
final_results <- df %>%
  dplyr::select(gene, gene_symbol, gene_biotype, everything()) %>%
  dplyr::select(-ensembl_id_clean)

# 08 Replace NA symbols with "NA_Symbol" for clarity
final_results$gene_symbol[is.na(final_results$gene_symbol)] <- "NA_Symbol"

# 09 Define Expression Categories (Corrected logic)
df <- df %>%
  mutate(expression_status = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "Over-expressed",
    padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
    TRUE ~ "Not-Significant"
  ))

# 10 Filter: Protein Coding Genes
protein_coding_genes <- df %>%
  dplyr::filter(gene_biotype == "protein_coding")

# 11 Filter: Protein Non-Coding Genes
protein_non_coding_genes <- df %>%
  dplyr::filter(gene_biotype != "protein_coding" | is.na(gene_biotype))

# 12 Filter: Over-expressed Genes
over_expressed_genes <- df %>%
  dplyr::filter(expression_status == "Over-expressed")

# 13 Final Exports
write.csv(final_results, "02 Output/Annotated_Full_Results.csv", row.names = FALSE)
write.csv(protein_coding_genes, "02 Output/1_Protein_Coding_Genes.csv", row.names = FALSE)
write.csv(protein_non_coding_genes, "02 Output/2_Protein_Non_Coding_Genes.csv", row.names = FALSE)
write.csv(over_expressed_genes, "02 Output/3_Over_expressed_Genes.csv", row.names = FALSE)