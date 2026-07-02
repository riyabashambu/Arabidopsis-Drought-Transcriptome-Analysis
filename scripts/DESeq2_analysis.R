############################################################
# Arabidopsis Drought Transcriptome Analysis (DESeq2 Pipeline)
#
# Author: Riya
# Organism: Arabidopsis thaliana (TAIR10)
# Conditions: Control vs Drought
############################################################

############################################################
# STEP 0: Setup
############################################################

setwd("//wsl.localhost/Ubuntu-24.04/home/riya/Arabidopsis-Drought-Transcriptome-Analysis")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

############################################################
# STEP 1: Load Libraries
############################################################

library(DESeq2)

library(dplyr)
library(tibble)
library(readr)
library(tidyr)

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)

############################################################
# STEP 2: Load featureCounts data
############################################################

counts <- read.delim(
  "results/counts/gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

dim(counts)
colnames(counts)

############################################################
# STEP 3: Build count matrix
############################################################

count_matrix <- counts[, c(1, 7:12)]
rownames(count_matrix) <- count_matrix$Geneid
count_matrix <- count_matrix[, -1]

colnames(count_matrix) <- c(
  "Drought_1","Drought_2","Drought_3",
  "Control_1","Control_2","Control_3"
)

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

############################################################
# STEP 4: Metadata
############################################################

sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(
    c("Drought","Drought","Drought",
      "Control","Control","Control"),
    levels = c("Control","Drought")
  )
)

############################################################
# STEP 5: DESeq2 object
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

############################################################
# STEP 6: Filtering + DE analysis
############################################################

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds)

############################################################
# STEP 7: Variance stabilization + PCA
############################################################

vsd <- vst(dds, blind = FALSE)

pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()

ggsave("results/figures/PCA_plot.png")

############################################################
# STEP 8: Sample distance heatmap
############################################################

sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         main = "Sample Distance Heatmap")

############################################################
# STEP 9: Differential expression analysis
############################################################

res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$GeneID <- rownames(res_df)

res_df <- res_df %>% dplyr::select(GeneID, everything())
res_df <- res_df %>% filter(!is.na(padj))

sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

upregulated <- sig_genes %>% filter(log2FoldChange > 1)
downregulated <- sig_genes %>% filter(log2FoldChange < -1)

write.csv(res_df, "results/deseq2/all_genes.csv", row.names = FALSE)
write.csv(sig_genes, "results/deseq2/significant_genes.csv", row.names = FALSE)

############################################################
# STEP 10: Volcano plot
############################################################

EnhancedVolcano(
  res_df,
  lab = NA,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Drought vs Control"
)

############################################################
# STEP 11: GO Enrichment
############################################################

sig_genes_ids <- sig_genes$GeneID

gene_map <- bitr(
  sig_genes_ids,
  fromType = "TAIR",
  toType = "ENTREZID",
  OrgDb = org.At.tair.db
)

gene_map <- gene_map %>% filter(!is.na(ENTREZID))

go_bp <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.At.tair.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

dotplot(go_bp, showCategory = 20)

############################################################
# STEP 12: KEGG (robust version)
############################################################

genes <- unique(gene_map$ENTREZID)

kegg <- tryCatch({
  enrichKEGG(
    gene = genes,
    organism = "ath",
    pvalueCutoff = 0.1
  )
}, error = function(e) NULL)

if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  dotplot(kegg, showCategory = 15)
} else {
  cat("KEGG not significant for this dataset\n")
}


############################################################
# STEP 13: Advanced GO Enrichment Visualizations
############################################################

############################################################
# GO Barplot
############################################################

png(
  "results/figures/GO_BP_barplot.png",
  width = 2200,
  height = 1800,
  res = 300
)

barplot(
  go_bp,
  showCategory = 20,
  title = "GO Biological Process Enrichment"
)

dev.off()

png(
  "results/figures/GO_BP_dotplot.png",
  width = 2200,
  height = 1800,
  res = 300
)

dotplot(
  go_bp,
  showCategory = 20
)

dev.off()

############################################################
# GO Gene-Concept Network
############################################################

png(
  "results/figures/GO_BP_cnetplot.png",
  width = 2600,
  height = 2200,
  res = 300
)

cnetplot(
  go_bp,
  showCategory = 5
)

dev.off()

############################################################
# GO Enrichment Map
############################################################

pairwise_termsim(go_bp) -> go_bp_sim

png(
  "results/figures/GO_BP_emapplot.png",
  width = 2400,
  height = 2000,
  res = 300
)

emapplot(
  go_bp_sim,
  showCategory = 20
)

dev.off()

############################################################
# Step 14: Heatmap of Top Differentially Expressed Genes
############################################################

#-----------------------------------------------------------
# Select top 50 significant DEGs
#-----------------------------------------------------------

top50 <- head(sig_genes$GeneID, 50)

#-----------------------------------------------------------
# Extract normalized expression values
#-----------------------------------------------------------

vsd_mat <- assay(vsd)

top50_matrix <- vsd_mat[top50, ]

#-----------------------------------------------------------
# Scale expression by gene (row-wise Z-score)
#-----------------------------------------------------------

top50_matrix_scaled <- t(scale(t(top50_matrix)))

#-----------------------------------------------------------
# Sample annotation
#-----------------------------------------------------------

annotation_col <- data.frame(
  Condition = sample_info$condition
)

rownames(annotation_col) <- rownames(sample_info)

#-----------------------------------------------------------
# Heatmap colors
#-----------------------------------------------------------

heat_colors <- colorRampPalette(
  c("navy", "white", "firebrick3")
)(100)

#-----------------------------------------------------------
# Save heatmap
#-----------------------------------------------------------

############################################################
# Top 50 Differentially Expressed Genes Heatmap
############################################################

# Create figures directory
dir.create("results/figures",
           recursive = TRUE,
           showWarnings = FALSE)

# Open PNG device
png(
  filename = "results/figures/Top50_DEGs_heatmap.png",
  width = 2400,
  height = 3000,
  res = 300
)

# Draw heatmap
pheatmap(
  top50_matrix_scaled,
  color = heat_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 12,
  border_color = NA,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

############################################################
# END OF PIPELINE
############################################################