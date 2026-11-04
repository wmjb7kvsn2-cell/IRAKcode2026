#!/usr/bin/env Rscript

# Targeted THP-1 RNA-seq heatmaps from RSEM TPM.
# Supplementary scale: log2(TPM + 1), no gene-wise scaling.
# Fig. 6D-F-style scale: gene-wise Z-score across all nine samples.
# Reconstructed workflow; verify exact original gene lists and row ordering.

suppressPackageStartupMessages({
  library(pheatmap)
})

dir.create("results/heatmaps", recursive = TRUE, showWarnings = FALSE)

sample_order <- c(
  "WT_1", "WT_2", "WT_3",
  "D320Y_1", "D320Y_2", "D320Y_3",
  "Q29P_1", "Q29P_2", "Q29P_3"
)

tpm_file <- "data/rsem_tpm.tsv"
if (!file.exists(tpm_file)) stop("Missing input: ", tpm_file)

tpm_df <- read.delim(tpm_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(sample_order %in% colnames(tpm_df))) {
  stop("TPM matrix must contain all nine sample columns with the expected names.")
}

gene_ids <- tpm_df[[1]]
tpm <- as.matrix(tpm_df[, sample_order, drop = FALSE])
rownames(tpm) <- gene_ids
storage.mode(tpm) <- "numeric"

if (anyNA(tpm) || any(tpm < 0)) {
  stop("TPM matrix contains missing or negative values.")
}
if (anyDuplicated(rownames(tpm))) {
  stop("Duplicated gene identifiers detected; resolve explicitly.")
}

log_tpm <- log2(tpm + 1)

gene_sets <- c(
  NFkB = "data/gene_sets/NFkB_64_genes.txt",
  MAPK = "data/gene_sets/MAPK_64_genes.txt",
  Interferon = "data/gene_sets/Interferon_30_genes.txt"
)

for (set_name in names(gene_sets)) {
  f <- gene_sets[[set_name]]
  if (!file.exists(f)) {
    warning("Skipping ", set_name, ": missing gene list ", f)
    next
  }

  genes <- unique(trimws(readLines(f, warn = FALSE)))
  genes <- genes[nzchar(genes)]
  present <- genes[genes %in% rownames(log_tpm)]
  missing <- setdiff(genes, present)

  writeLines(missing, paste0("results/heatmaps/", set_name,
                             "_genes_not_found.txt"))

  if (length(present) == 0) {
    warning("No genes from ", set_name, " matched the TPM matrix.")
    next
  }

  mat_unscaled <- log_tpm[present, , drop = FALSE]

  # Supplementary heatmap: log2(TPM + 1), without row scaling.
  pdf(paste0("results/heatmaps/", set_name,
             "_log2TPM_unscaled.pdf"), width = 8, height = 10)
  pheatmap(
    mat_unscaled,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = paste0(set_name, ": log2(TPM + 1)")
  )
  dev.off()

  # Fig. 6D-F-style heatmap: gene-wise Z-score across the nine samples.
  row_sd <- apply(mat_unscaled, 1, sd)
  keep <- is.finite(row_sd) & row_sd > 0
  mat_for_z <- mat_unscaled[keep, , drop = FALSE]
  if (nrow(mat_for_z) > 0) {
    mat_z <- t(scale(t(mat_for_z)))
    pdf(paste0("results/heatmaps/", set_name,
               "_log2TPM_geneWise_Zscore.pdf"), width = 8, height = 10)
    pheatmap(
      mat_z,
      scale = "none", # Already Z-scored explicitly.
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      main = paste0(set_name, ": gene-wise Z-score")
    )
    dev.off()
  }
}

capture.output(sessionInfo(), file = "results/heatmaps/sessionInfo.txt")
cat("Heatmap generation completed.\n")
