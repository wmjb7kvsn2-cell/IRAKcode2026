#!/usr/bin/env Rscript

# Differential expression analysis for THP-1 IRAK4 RNA-seq.
# Input: raw integer gene counts; NOT TPM.
# Recorded historical method: DESeq2 v1.34.0; DEG threshold Padj < 0.05
# and |log2FC| > 1.5. This is a reconstructed workflow, not recovered
# verbatim source code.

suppressPackageStartupMessages({
  library(DESeq2)
})

dir.create("results/DESeq2", recursive = TRUE, showWarnings = FALSE)

sample_order <- c(
  "WT_1", "WT_2", "WT_3",
  "D320Y_1", "D320Y_2", "D320Y_3",
  "Q29P_1", "Q29P_2", "Q29P_3"
)

count_file <- "data/counts_matrix.tsv"
if (!file.exists(count_file)) stop("Missing input: ", count_file)

counts_df <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts_df) != 10) {
  stop("Expected 1 gene-ID column + 9 sample columns in counts_matrix.tsv")
}
gene_ids <- counts_df[[1]]
count_mat <- as.matrix(counts_df[, sample_order, drop = FALSE])
rownames(count_mat) <- gene_ids
storage.mode(count_mat) <- "numeric"

if (anyNA(count_mat) || any(count_mat < 0) ||
    any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("DESeq2 input must contain non-negative raw integer counts.")
}
storage.mode(count_mat) <- "integer"
if (anyDuplicated(rownames(count_mat))) {
  stop("Duplicated gene identifiers detected; resolve explicitly before analysis.")
}

group <- factor(
  c(rep("WT", 3), rep("D320Y", 3), rep("Q29P", 3)),
  levels = c("WT", "D320Y", "Q29P")
)
coldata <- data.frame(row.names = sample_order, group = group)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ group
)

# Default DESeq2 independent filtering and Cook's distance handling are retained.
dds <- DESeq(dds)

write.csv(
  as.data.frame(counts(dds, normalized = TRUE)),
  "results/DESeq2/normalized_counts.csv",
  quote = TRUE
)

contrasts <- list(
  D320Y_vs_WT = c("group", "D320Y", "WT"),
  Q29P_vs_WT  = c("group", "Q29P", "WT"),
  D320Y_vs_Q29P = c("group", "D320Y", "Q29P")
)

for (nm in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[nm]])
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

  write.csv(res_df, file = paste0("results/DESeq2/", nm, "_all.csv"),
            row.names = FALSE, na = "")

  sig <- subset(res_df, !is.na(padj) & padj < 0.05 &
                  !is.na(log2FoldChange) & abs(log2FoldChange) > 1.5)
  write.csv(sig, file = paste0("results/DESeq2/", nm,
                               "_Padj0.05_absLFC1.5.csv"),
            row.names = FALSE, na = "")
}

capture.output(sessionInfo(), file = "results/DESeq2/sessionInfo.txt")
cat("DESeq2 analysis completed. DESeq2 version:",
    as.character(packageVersion("DESeq2")), "\n")
