# notebooks/02_make_signatures.R (or any R session)

library(readr)
library(dplyr)
library(tidyr)
library(ashr)         # shrinkage for DESeq2
library(VennDiagram)  # simple Venn diagrams
library(tidyverse)
library(DESeq2)

# 1) Load
counts <- read_tsv("datasets/salmon.merged.gene_counts.tsv")
# Columns: gene_id, gene_name, then samples (e.g., d15_Tcm_REP1 ... Naive_REP2)

# 2) Basic shape and a peek
dim(counts)           # rows/cols
head(counts, 3)       # first genes
glimpse(counts)       # column types

# 3) Identify the sample columns
sample_cols <- setdiff(names(counts), c("gene_id", "gene_name"))



# Load full counts
counts <- read_tsv("datasets/salmon.merged.gene_counts.tsv")

# Keep d8 subsets plus Naive
tcm_cols <- c("d8_Tcm_REP1","d8_Tcm_REP2","d8_Tcm_REP3")
tem_cols <- c("d8_Tem_REP1","d8_Tem_REP2","d8_Tem_REP3")
tpm_cols <- c("d8_Tpm_REP1","d8_Tpm_REP2","d8_Tpm_REP3")
naive_cols <- c("Naive_REP1", "Naive_REP2")
use_cols <- c(tcm_cols, tem_cols, tpm_cols, naive_cols)
analysis_counts <- counts |> select(gene_id, gene_name, all_of(use_cols))

# Build count matrix and metadata
count_mat <- analysis_counts |>
  select(all_of(use_cols)) |>
  mutate(across(everything(), ~ round(.))) |>
  as.matrix()
storage.mode(count_mat) <- "integer"  # DESeq2 requires integer counts
rownames(count_mat) <- analysis_counts$gene_id

sample_info <- tibble(
  sample = use_cols,
  condition = c(rep("Tcm", length(tcm_cols)),
                rep("Tem", length(tem_cols)),
                rep("Tpm", length(tpm_cols)),
                rep("Naive", length(naive_cols)))
) |>
  mutate(condition = factor(condition, levels = c("Naive", "Tcm", "Tem", "Tpm"))) |>  # Naive as reference
  column_to_rownames("sample")

# DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ condition
)

# Filter low counts (keep genes with ≥10 total)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DE analysis
dds <- DESeq(dds)

# Retrieve contrasts vs Naive with ashr shrinkage
res_Tcm_vs_naive <- lfcShrink(dds, contrast = c("condition", "Tcm", "Naive"), type = "ashr")
res_Tem_vs_naive <- lfcShrink(dds, contrast = c("condition", "Tem", "Naive"), type = "ashr")
res_Tpm_vs_naive <- lfcShrink(dds, contrast = c("condition", "Tpm", "Naive"), type = "ashr")

# Quick summaries
summary(res_Tcm_vs_naive)
summary(res_Tem_vs_naive)
summary(res_Tpm_vs_naive)

# Top hits example
res_Tcm_vs_naive |>
  as_tibble(rownames = "gene_id") |>
  left_join(analysis_counts |> select(gene_id, gene_name), by = "gene_id") |>
  arrange(padj) |>
  slice_head(n = 20)

# Euler diagram of up-regulated genes (padj < 0.05, log2FC > 0)
sig_Tcm <- res_Tcm_vs_naive |> as_tibble(rownames = "gene_id") |>
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 0)
sig_Tem <- res_Tem_vs_naive |> as_tibble(rownames = "gene_id") |>
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 0)
sig_Tpm <- res_Tpm_vs_naive |> as_tibble(rownames = "gene_id") |>
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 0)

sig_sets <- list(
  Tcm = sig_Tcm$gene_id,
  Tem = sig_Tem$gene_id,
  Tpm = sig_Tpm$gene_id
)
# Simple Venn diagram (upregulated gene overlaps)
grid.newpage()
venn.plot <- venn.diagram(
  x = sig_sets,
  filename = NULL,
  fill = c("#4C78A8", "#F58518", "#54A24B"),
  alpha = 0.6,
  cex = 1.2,
  cat.cex = 1.2,
  lwd = 1.5
)
grid.draw(venn.plot)

# Exclusive signatures: genes up in only one subset (padj < 0.05, log2FC > 0)
tcm_only <- setdiff(sig_sets$Tcm, union(sig_sets$Tem, sig_sets$Tpm))
tem_only <- setdiff(sig_sets$Tem, union(sig_sets$Tcm, sig_sets$Tpm))
tpm_only <- setdiff(sig_sets$Tpm, union(sig_sets$Tcm, sig_sets$Tem))

make_sig <- function(ids, label) {
  tibble(gs_name = label, EnsemblID = ids) |>
    left_join(analysis_counts |> select(gene_id, gene_name), by = c("EnsemblID" = "gene_id")) |>
    mutate(gene_symbol = gene_name) |>
    select(gs_name, gene_symbol, EnsemblID)
}

exclusive_sigs <- bind_rows(
  make_sig(tcm_only, "Tcm_only"),
  make_sig(tem_only, "Tem_only"),
  make_sig(tpm_only, "Tpm_only")
)

write_csv(exclusive_sigs, "signatures_exclusive.csv")
print("Wrote exclusive signatures to signatures_exclusive.csv")
