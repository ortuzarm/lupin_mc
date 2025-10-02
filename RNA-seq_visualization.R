#RNA-seq visualization
#maiteortuzar@usal.es

#Libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(matrixStats)
library(RColorBrewer)
library(clusterProfiler)
library(AnnotationDbi)
library(fgsea)
library(msigdbr)

#Data
counts <- read.csv2("gene_count_1.csv", row.names = 1, check.names = FALSE)

condition <- factor(c(
  rep("Control_uninoculated", 3),
  rep("SynCom_1", 3),
  rep("SynCom_2", 3),
  rep("SynCom_3", 3),
  rep("SynCom_4", 3),
  rep("SynCom_5", 3),
  rep("SynCom_6", 3),
  rep("SynCom_7", 3)
))
colData <- data.frame(row.names = colnames(counts), condition = condition)

#DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

#Significant genes
res_sig <- res %>% as.data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
sig_genes <- rownames(res_sig)

write.csv(as.data.frame(res), "DESeq2_results.csv")
write.csv(res_sig, "DESeq2_significant_results.csv")

#rlog transformation and visualizations
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

#PCA
pca_df <- prcomp(t(rld_mat))$x[, 1:2] %>% as.data.frame() %>% mutate(Condition = condition)
ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA (rlog)")
ggsave("PCA_rlog.pdf", 6, 5)

#Sample-to-sample distances
sample_dists <- dist(t(rld_mat))
pheatmap(as.matrix(sample_dists),
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         main = "Sample-to-sample distances")
ggsave("Distance_heatmap.pdf", 6, 5)

#rlog means by condition
treatment_groups <- list(
  Control_uninoculated = grep("^Control_uninoculated", colnames(counts), value = TRUE),
  SynCom_1 = grep("^SynCom_1", colnames(counts), value = TRUE),
  SynCom_2 = grep("^SynCom_2", colnames(counts), value = TRUE),
  SynCom_3 = grep("^SynCom_3", colnames(counts), value = TRUE),
  SynCom_4 = grep("^SynCom_4", colnames(counts), value = TRUE),
  SynCom_5 = grep("^SynCom_5", colnames(counts), value = TRUE),
  SynCom_6 = grep("^SynCom_6", colnames(counts), value = TRUE),
  SynCom_7 = grep("^SynCom_7", colnames(counts), value = TRUE)
)
mean_data <- sapply(treatment_groups, function(samples) {
  rowMeans(rld_mat[sig_genes, samples])
})
mean_data <- as.data.frame(mean_data)

#Top 50 variable genes heatmap
mean_data <- mean_data[, setdiff(colnames(mean_data), "Control_uninoculated")]

#Select most variable genes
var_rank <- rowVars(as.matrix(mean_data))
top_n <- 50
top_idx <- order(var_rank, decreasing = TRUE)[seq_len(min(top_n, length(var_rank)))]
heat_mat <- mean_data[top_idx, ]
desired_order <- c("SynCom_1", "SynCom_2", "SynCom_3", 
                   "SynCom_4", "SynCom_5", "SynCom_6", "SynCom_7")
heat_mat <- heat_mat[, desired_order]

#Read file with genes and descriptions
gene_map <- read.csv2("code_genes_with_descriptions.csv")
mapped_names <- ifelse(
  rownames(heat_mat) %in% gene_map$gene_id,
  gene_map$gene_description[match(rownames(heat_mat), gene_map$gene_id)],
  rownames(heat_mat)
)
mapped_names[is.na(mapped_names)] <- rownames(heat_mat)[is.na(mapped_names)]
mapped_names <- make.unique(mapped_names)
rownames(heat_mat) <- mapped_names

#Column annotation
annot_col <- data.frame(Condition = colnames(heat_mat))
rownames(annot_col) <- colnames(heat_mat)
my_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)

#Heatmap
pheatmap(heat_mat, annotation_col = annot_col, scale = "row", color = my_pal,
         fontsize_row = 6, fontsize_col = 8,
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = paste("Top", top_n, "variable genes (excluding control)"))
ggsave("Heatmap_topVar_sigGenes_excluding_control.pdf", width = 6, height = 8)
scaled_mat <- t(scale(t(heat_mat)))
write.csv(scaled_mat, "Heatmap_topVar_sigGenes_excluding_control_scaled.csv")
write.csv(mapped_names, "mapped_names.csv", row.names = FALSE)
