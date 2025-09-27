# ============================================================
# R script: Unstranded Single-End RNA-Seq Analysis: PRJNA728240 Dataset
# ============================================================

# 1. Design matrix:
#    - Sample names: 3T3L1-D0-rep1, 3T3L1-D0-rep2, 3T3L1-D0-rep3,
#                    3T3L1-D2-rep1, 3T3L1-D2-rep2, 3T3L1-D2-rep3,
#                    3T3L1-D6-rep1, 3T3L1-D6-rep2, 3T3L1-D6-rep3
#    - Factor: condition (levels = d0, d2, d6), with d0 as reference.

# 2. Import counts:
#    - Input file: raw_counts.txt from featureCounts.
#    - Skip program header line.
#    - Remove first 5 annotation columns (Geneid, Chr, Start, End, Strand, Length).
#    - Keep Geneid as rownames.
#    - Rename columns to clean sample names.

# 3. QC: Pre-normalization
#    - Library size barplot (total counts per sample).
#    - Boxplot of log2(raw counts + 1).
#    - Sample correlation heatmap.

# 4. Filtering:
#    - Remove low-count genes (e.g., keep genes with ≥10 counts in at least 3 samples).

# 5. Normalization:
#    - Estimate size factors (DESeq2::estimateSizeFactors()).
#    - Apply vst() or rlog() transformation for visualization (PCA/heatmap).

# 6. QC: Post-normalization
#    - PCA plot (plotPCA).
#    - Sample-to-sample distance heatmap (pheatmap).
#    - MA plot for quick overview.

# 7. Model fitting:
#    - Design formula: ~ condition.
#    - Run DESeq2 model (no batch correction unless PCA shows batch structure).

# 8. DE analysis + Post-DE filtering:
#    - Contrasts: d6 vs d0, d2 vs d0.
#    - Shrink LFC with apeglm.
#    - Post-DE filtering with baseMean >=10.
#    - Exclude genes with very low average expression (padj=0)
#    - Save results (log2FC, p-value, padj).

# 9. Map Ensembl IDs to gene symbols/descriptive names

# 10. Volcano plots:
#    - Plot log2FC vs -log10(padj).
#    - Add thresholds (padj < 0.05, |log2FC| > 1).
#    - Label top 10 up and top 10 down genes.
#    - Export plots (PNG + PDF).

# 11. Export results:
#    - Normalized counts (counts(dds, normalized=TRUE)).
#    - DE results (raw + shrunken LFC).
#    - Save sessionInfo() for reproducibility.

# ----------------------------- Install and load dependencies -----------------------------

setwd("C:/Users/Dr. Shahgaldi/Desktop/PRJNA728240_analysis/3. Sent to prof Raj/ZBP1 Included")

# CRAN packages for plotting, data manipulation, and QC
req.cran <- c("ggplot2", "dplyr", "readr", "pheatmap", "RColorBrewer", "ggrepel")

for (p in req.cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

# Ensure BiocManager is installed for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Bioconductor packages required for DESeq2 pipeline
req.bioc <- c(
  "DESeq2",        # Differential expression analysis
  "apeglm",        # Log2 fold change shrinkage
  "sva",           # Batch correction (optional if needed)
  "EnhancedVolcano", # Volcano plots (optional, can also use ggplot2)
  "org.Mm.eg.db"   # Mouse gene annotation
)

for (p in req.bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

# Load all packages, suppressing startup messages
pkgs <- c(req.cran, req.bioc)

invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

# Global options / reproducibility
options(stringsAsFactors = FALSE)
set.seed(123)

# ----------------------------- 1. Design matrix and contrasts -----------------------------

# Define sample names exactly as in your counts file
sample_names <- c(
  "3T3L1-D0-rep1", "3T3L1-D0-rep2", "3T3L1-D0-rep3",
  "3T3L1-D2-rep1", "3T3L1-D2-rep2", "3T3L1-D2-rep3",
  "3T3L1-D6-rep1", "3T3L1-D6-rep2", "3T3L1-D6-rep3"
)

# Create colData dataframe
coldata <- data.frame(
  row.names = sample_names,
  condition = factor(
    rep(c("d0", "d2", "d6"), each = 3),
    levels = c("d0", "d2", "d6")  # d0 as reference
  )
)

# Display design matrix for verification
coldata

# Define contrasts for DESeq2
contrasts_list <- list(
  c("condition", "d6", "d0"),  # d6 vs d0
  c("condition", "d2", "d0")   # d2 vs d0
)

# Optional: print contrasts
contrasts_list

# ----------------------------- 2. Import raw counts -----------------------------

# Path to featureCounts output
counts_file <- "raw_counts.txt"

# Read counts, skipping the featureCounts program info line
counts_raw <- read.table(
  counts_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,       # Use Geneid as rownames
  check.names = FALSE, # Preserve column names exactly
  skip = 1,            # Skip program info line
  stringsAsFactors = FALSE
)

# Inspect original header / column names from featureCounts (helpful to verify order)
cat("Original header columns (first 15):\n")
print(head(colnames(counts_raw), 15))

# Remove first 5 annotation columns (Chr, Start, End, Strand, Length)
counts_matrix <- counts_raw[, -(1:5), drop = FALSE]

# Inspect the columns we've kept as count columns
cat("Columns kept as potential count columns:\n")
print(colnames(counts_matrix))

# Try to map the file's count columns to your `sample_names` (pattern match)
orig_cols_counts <- colnames(counts_matrix)

map_idx <- sapply(sample_names, function(sn) {
  which_idx <- which(grepl(sn, orig_cols_counts, fixed = TRUE))
  if (length(which_idx) == 1) return(which_idx)
  if (length(which_idx) > 1) return(which_idx[1])  # fallback if ambiguous
  return(NA_integer_)
})

mapping_df <- data.frame(expected = sample_names,
                         mapped_index = map_idx,
                         orig_name = orig_cols_counts[map_idx],
                         stringsAsFactors = FALSE)
cat("Mapping summary (expected -> mapped_index -> original column name):\n")
print(mapping_df)

# If mapping failed for any sample, but column counts match, fall back to positional rename;
# otherwise stop and ask for manual inspection.
if (any(is.na(map_idx))) {
  warning("Some sample_names did not auto-map by pattern.\n")
  if (ncol(counts_matrix) == length(sample_names)) {
    cat("Falling back to renaming by column order (counts columns == sample_names length).\n")
    colnames(counts_matrix) <- sample_names
  } else {
    stop("Automatic mapping failed AND counts columns != sample_names length. Inspect 'colnames(counts_matrix)' and 'sample_names'.")
  }
} else {
  # Reorder and rename columns to match sample_names
  counts_matrix <- counts_matrix[, map_idx, drop = FALSE]
  colnames(counts_matrix) <- sample_names
  cat("Columns successfully matched and reordered to sample_names using pattern matching.\n")
}

# Coerce count columns to integer numeric safely
counts_matrix[] <- lapply(counts_matrix, function(x) as.integer(as.character(x)))
counts_matrix <- as.matrix(counts_matrix)

# Check for any NAs introduced by coercion (should be 0)
na_count <- sum(is.na(counts_matrix))

cat("NA values introduced by coercion (should be 0):", na_count, "\n")

if (na_count > 0) {
  cat("Example NA positions (row, col):\n")
  print(head(which(is.na(counts_matrix), arr.ind = TRUE)))
  stop("There are NA values in counts_matrix after coercion — inspect the source file for non-numeric entries.")
}

# Final verification
cat("Import finished. counts_matrix dimensions (genes x samples):", dim(counts_matrix), "\n")
View(head(counts_matrix[, 1:min(9, ncol(counts_matrix))]))

# ----------------------------- 3. QC: Pre-normalization -----------------------------
# The idea is to check for library size differences, extreme outliers, or any sample issues before any normalization.
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ------------- Library size barplot (total counts per sample)
lib_sizes <- colSums(counts_matrix)

lib_sizes_df <- data.frame(
  sample = names(lib_sizes),
  total_counts = lib_sizes
)

ggplot(lib_sizes_df, aes(x = sample, y = total_counts)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Sizes (Raw Counts)", y = "Total Counts", x = "Sample")

ggsave("preQC_library_sizes.png", width = 7, height = 5)
ggsave("preQC_library_sizes.pdf", width = 7, height = 5)

# ------------- Boxplot of log2(counts + 1)

# Adding +1 avoids log2(0) errors.
logcounts <- log2(counts_matrix + 1)

# Convert to data.frame and preserve Gene IDs
logcounts_df <- as.data.frame(logcounts)
logcounts_df$Gene <- rownames(logcounts_df)

# Melt using Gene as ID variable
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
library(reshape2)

logcounts_df <- melt(logcounts_df, id.vars = "Gene")
colnames(logcounts_df) <- c("Gene", "Sample", "log2Count")

# Plot
p_boxplot <- ggplot(logcounts_df, aes(x = Sample, y = log2Count)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of log2(raw counts + 1)", y = "log2(count + 1)", x = "Sample")

# Display plot
print(p_boxplot)

# Save plot as PNG and PDF
ggsave("preQC_boxplot_log2counts.png", plot = p_boxplot, width = 7, height = 5)
ggsave("preQC_boxplot_log2counts.pdf", plot = p_boxplot, width = 7, height = 5)

# ------------- Sample correlation heatmap
sample_cor <- cor(counts_matrix, method = "pearson")

# Save PNG directly
pheatmap(
  sample_cor,
  color = colorRampPalette(brewer.pal(9, "Blues"))(100),
  main = "Sample Correlation (Raw Counts)",
  display_numbers = TRUE,
  filename = "preQC_sample_correlation_heatmap.png",
  width = 6,
  height = 5
)

# Save PDF by opening/closing device manually (no filename arg inside pheatmap)
pdf("preQC_sample_correlation_heatmap.pdf", width = 6, height = 5)
pheatmap(
  sample_cor,
  color = colorRampPalette(brewer.pal(9, "Blues"))(100),
  main = "Sample Correlation (Raw Counts)",
  display_numbers = TRUE
)
dev.off()

# ----------------------------- 4. Filtering: Low-count genes -----------------------------

# ------------- Principled filtering using edgeR

# Quick ad-hoc thresholds to see how many genes would be kept by simple rules:
thresholds <- c(2, 5, 10, 20)
threshold_keep_counts <- sapply(thresholds, function(t) sum(rowSums(counts_matrix >= t) >= 3))

# Print a small table: for each threshold t, how many genes have >= t counts in at least 3 samples
# "threshold 2" will keep the most genes (least stringent). "threshold 20" the fewest (most stringent).
# Large drops between thresholds (e.g. 2 -> 10) indicate many very-low-count genes.
# Ad-hoc thresholds are simple but not adaptive to library size; prefer filterByExpr for a principled choice.
thresh_tbl <- data.frame(threshold = thresholds, genes_kept = threshold_keep_counts)
print(thresh_tbl)

# Now use edgeR::filterByExpr (recommended, adaptive to library sizes & group design)
if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}
library(edgeR)

# filterByExpr expects a counts matrix (genes x samples) and a group factor
keep_edgeR <- filterByExpr(counts_matrix, group = coldata$condition)

# Diagnostics: how many genes are kept by filterByExpr
# - If filterByExpr retains a reasonable number (e.g. thousands), use it and proceed.
# - If it retains *very few* genes (< ~1000) AND you expect many expressed genes biologically, inspect library sizes and sample quality.
# - If you prefer a manual cutoff, pick the ad-hoc threshold that gives a similar number of genes as filterByExpr.
cat("Genes before filtering:", nrow(counts_matrix), "\n")
cat("Genes kept by filterByExpr:", sum(keep_edgeR), "\n")
cat("Percentage retained by filterByExpr:", round(100 * sum(keep_edgeR) / nrow(counts_matrix), 1), "%\n")

# Apply the filter to a DESeqDataSet (create dds if not already)
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = coldata, design = ~ condition)
dds_filtered <- dds[keep_edgeR, ]

# Save short report to disk (optional)
write.csv(thresh_tbl, file = "filter_threshold_summary.csv", row.names = FALSE)
cat("Filtered DESeqDataSet saved in 'dds_filtered' object (in memory).\n")

# ----------------------------- 5. Normalization -----------------------------
library(DESeq2)

# Estimate size factors for normalization.
# Normalizes for library size differences. 
# You don’t need log2 or CPM here, DESeq2 handles it internally.
dds_filtered <- estimateSizeFactors(dds_filtered)

# Access normalized counts 
norm_counts <- counts(dds_filtered, normalized = TRUE)
write.csv(norm_counts, "normalized_counts.csv", quote = FALSE)
saveRDS(norm_counts, "normalized_counts.rds")

# Variance-stabilizing transformation (VST) (for visualization/QC only)
# vst() or rlog() are mainly for visualization and QC (PCA, heatmaps). 
# Do not use these transformed counts for DE testing; DESeq2 uses the raw counts internally.
# Alternatively, you could use rlog transformation:
# rlog_counts <- rlog(dds_filtered, blind = TRUE)

vst_counts <- vst(dds_filtered, blind = TRUE)  # blind = TRUE ignores design; good for QC
saveRDS(vst_counts, file = "vst_dds_filtered.rds")

# Extract the transformed VST matrix for PCA/heatmaps
vst_mat <- assay(vst_counts)

# Quick check
dim(vst_mat)  # should be genes x samples
View(head(vst_mat[, 1:9]))

# ----------------------------- 6. QC: Post-normalization -----------------------------
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(reshape2)

# --------------------------- PCA plot
# Uses assay(vst_counts) explicitly and ensures sample/order matching with coldata.

# Compute PCA on VST matrix (samples as rows)
vst_mat <- assay(vst_counts)                     # genes x samples
pca_res <- prcomp(t(vst_mat), center = TRUE)     # samples as rows

# Percent variance explained for axis labels
percentVar <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
xlab <- sprintf("PC1 (%.1f%%)", percentVar[1])
ylab <- sprintf("PC2 (%.1f%%)", percentVar[2])

# Build PCA dataframe ensuring Condition is matched to sample order
sample_order <- colnames(vst_mat)

pca_df <- data.frame(
  Sample = sample_order,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Condition = coldata[sample_order, "condition", drop = TRUE]
)

pca_df$Condition <- factor(pca_df$Condition, levels = levels(coldata$condition))

# Plot with ggplot2 + ggrepel for readable labels

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(min.segment.length = 0.01, size = 3) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA of VST-normalized counts", x = xlab, y = ylab) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")

# Display and save
print(p_pca)
ggsave("postQC_PCA_plot.png", plot = p_pca, width = 7, height = 5, dpi = 300)
ggsave("postQC_PCA_plot.pdf", plot = p_pca, width = 7, height = 5)


# --------------------------- Sample-to-sample distance heatmap

sample_dist <- dist(t(vst_mat))
sample_dist_mat <- as.matrix(sample_dist)
rownames(sample_dist_mat) <- colnames(vst_mat)
colnames(sample_dist_mat) <- colnames(vst_mat)

# Consistent color scale
heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(100)

# PNG
pheatmap(sample_dist_mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = heatmap_colors,
         main = "Sample-to-Sample Distance (VST)",
         filename = "postQC_sample_distance_heatmap.png")

# Also save PDF
pdf("postQC_sample_distance_heatmap.pdf", width = 6, height = 5)
pheatmap(sample_dist_mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = heatmap_colors,
         main = "Sample-to-Sample Distance (VST)")
dev.off()

# --------------------------- Library sizes after normalization

lib_sizes_norm <- colSums(norm_counts)

lib_sizes_norm_df <- data.frame(
  sample = names(lib_sizes_norm),
  total_counts = lib_sizes_norm
)

ggplot(lib_sizes_norm_df, aes(x = sample, y = total_counts)) +
  geom_bar(stat = "identity", fill = "darkseagreen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Sizes (Normalized Counts)", 
       y = "Total Normalized Counts", x = "Sample")

ggsave("postQC_library_sizes.png", width = 7, height = 5)
ggsave("postQC_library_sizes.pdf", width = 7, height = 5)

# --------------------------- Boxplot of VST-transformed counts
# VST already stabilizes variance, so distributions should be aligned

vst_df <- as.data.frame(vst_mat)
vst_df$Gene <- rownames(vst_df)

# Melt to long format
library(reshape2)
vst_melt <- melt(vst_df, id.vars = "Gene")
colnames(vst_melt) <- c("Gene", "Sample", "VST_Count")

# Plot boxplot
p_boxplot_vst <- ggplot(vst_melt, aes(x = Sample, y = VST_Count)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgreen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of VST-normalized counts", 
       y = "VST-normalized expression", x = "Sample")

print(p_boxplot_vst)

# Save plots
ggsave("postQC_boxplot_vst_counts.png", plot = p_boxplot_vst, width = 7, height = 5, dpi = 300)
ggsave("postQC_boxplot_vst_counts.pdf", plot = p_boxplot_vst, width = 7, height = 5)

# --------------------------- Sample correlation heatmap (VST-normalized counts)

# Compute Pearson correlation across samples using VST matrix
sample_cor_vst <- cor(vst_mat, method = "pearson")

# PNG version
pheatmap(
  sample_cor_vst,
  color = colorRampPalette(brewer.pal(9, "Greens"))(100),
  main = "Sample Correlation (VST-normalized)",
  display_numbers = TRUE,
  filename = "postQC_sample_correlation_heatmap.png",
  width = 6,
  height = 5
)

# PDF version
pdf("postQC_sample_correlation_heatmap.pdf", width = 6, height = 5)
pheatmap(
  sample_cor_vst,
  color = colorRampPalette(brewer.pal(9, "Greens"))(100),
  main = "Sample Correlation (VST-normalized)",
  display_numbers = TRUE
)
dev.off()

# ----------------------------- 7. Model fitting: DESeq2 -----------------------------

library(DESeq2)

# Use the filtered and normalized DESeqDataSet (dds_filtered)
# Design formula: ~ condition
# If batch effects are detected in PCA, you could add them here, e.g., ~ batch + condition

# Ensure reference level (d0)
dds_filtered$condition <- relevel(dds_filtered$condition, ref = "d0")

# Run DESeq model
dds_filtered <- DESeq(dds_filtered, quiet = FALSE)

# Plot dispersion estimates
plotDispEsts(dds_filtered, main = "Dispersion estimates")

# Extract fitted normalized counts (optional)
fitted_counts <- counts(dds_filtered, normalized = TRUE)

# Quick check of size factors and conditions
print(sizeFactors(dds_filtered))
print(colData(dds_filtered))

# Save the DESeqDataSet object for downstream steps
saveRDS(dds_filtered, file = "dds_filtered_DESeq2_model.rds")

# ----------------------------- 8. DE analysis + Post-DE filtering -----------------------------

library(DESeq2)
library(apeglm)

# Ensure reference level is correct (already done in model fitting)
dds_filtered$condition <- relevel(dds_filtered$condition, ref = "d0")

# Map contrasts to coefficients in the DESeq2 model
# Typically: "Intercept (baseline expression for d0)", "condition_d2_vs_d0", "condition_d6_vs_d0"
coef_map <- c(
  d2_vs_d0 = "condition_d2_vs_d0",
  d6_vs_d0 = "condition_d6_vs_d0"
)

coef(dds_filtered)

# Initialize list to store DE results
de_results <- list()

for (contrast_name in names(coef_map)) {
  
  coef_name <- coef_map[[contrast_name]]
  
  # Shrink log2 fold changes using apeglm
  res_shrunk <- lfcShrink(dds_filtered, coef = coef_name, type = "apeglm")
  
  # Add post-DE filter: keep only genes with baseMean >=10
  keep_baseMean <- res_shrunk$baseMean >= 10
  res_shrunk <- res_shrunk[keep_baseMean, ]
  
  # Convert to data frame
  res_df <- as.data.frame(res_shrunk)
  
  # Some columns may be missing if no genes pass filtering; handle safely
  expected_cols <- c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  missing_cols <- setdiff(expected_cols, colnames(res_df))
  if(length(missing_cols) > 0) {
    for(col in missing_cols) res_df[[col]] <- NA
  }
  
  # Add gene IDs
  res_df$GeneID <- rownames(res_df)
  
  # Reorder columns safely
  res_df <- res_df[, c("GeneID", expected_cols)]
  
  # Save to CSV
  write.csv(res_df, paste0("DE_results_", contrast_name, ".csv"), row.names = FALSE, quote = FALSE)
  
  # Store in list
  de_results[[contrast_name]] <- res_df
}

# Quick check: number of genes after filtering for each contrast
for (contrast_name in names(de_results)) {
  cat("Genes retained after baseMean >=10 for", contrast_name, ":", nrow(de_results[[contrast_name]]), "\n")
}

# Verify by seeing: DE_results_d2_vs_d0.csv; DE_results_d6_vs_d0.csv; de_results_all_contrasts.rds

head(de_results$d2_vs_d0)
head(de_results$d6_vs_d0)

# Save all results together
saveRDS(de_results, file = "DE_results_all_contrasts.rds")

cat("Post-DE filtering complete. Filtered results saved for each contrast.\n")

# ----------------------------- 9. Map Ensembl IDs to gene symbols -----------------------------
library(AnnotationDbi)
library(org.Mm.eg.db)

# Function to strip version numbers from Ensembl IDs
strip_version <- function(ensembl_ids) {
  sub("\\..*$", "", ensembl_ids)
}

# Loop through DE results to map Ensembl IDs to gene symbols
for (contrast_name in names(de_results)) {
  
  res_df <- de_results[[contrast_name]]
  
  # Strip version numbers from Ensembl IDs
  res_df$ENSEMBL_nover <- strip_version(res_df$GeneID)
  
  # Map to gene symbols
  gene_mapping <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = res_df$ENSEMBL_nover,
    columns = c("SYMBOL"),
    keytype = "ENSEMBL"
  )
  
  # Remove duplicates (keep first mapping)
  gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENSEMBL), ]
  
  # Merge mapping back to DE results
  res_df <- merge(
    res_df,
    gene_mapping[, c("ENSEMBL", "SYMBOL")],
    by.x = "ENSEMBL_nover",
    by.y = "ENSEMBL",
    all.x = TRUE
  )
  
  # Fallback: use Ensembl ID if SYMBOL is NA
  res_df$GeneSymbol <- ifelse(is.na(res_df$SYMBOL), res_df$GeneID, res_df$SYMBOL)
  
  # Optional: reorder columns for clarity
  res_df <- res_df[, c("GeneID", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  
  # Save updated DE results with gene symbols
  de_results[[contrast_name]] <- res_df
  write.csv(res_df, paste0("DE_results_", contrast_name, "_withSymbols.csv"), 
            row.names = FALSE, quote = FALSE)
}

cat("Ensembl IDs successfully mapped to gene symbols for all contrasts.\n")

# ----------------------------- 10. Volcano plots (padj < 0.05, |log2FC| > 1) -----------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)

# Volcano plot thresholds
padj_thresh <- 0.05
log2FC_thresh <- 1

set.seed(123)  # For reproducible jitter

for (contrast_name in names(de_results)) {
  
  res_df <- de_results[[contrast_name]]
  
  # Remove genes with NA in log2FC or padj to avoid plotting warnings
  res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$padj), ]
  
  # Handle very small padj values to avoid Inf in -log10
  res_df$padj <- pmax(res_df$padj, 1e-320)
  
  # Compute -log10 adjusted p-value
  res_df$negLog10Padj <- -log10(res_df$padj)
  
  # Flag significance based on thresholds
  res_df$Significant <- "NotSig"
  res_df$Significant[res_df$padj < padj_thresh & res_df$log2FoldChange > log2FC_thresh] <- "Up"
  res_df$Significant[res_df$padj < padj_thresh & res_df$log2FoldChange < -log2FC_thresh] <- "Down"
  
  # Filter genes that have proper symbols (avoid ENSMUSG IDs)
  labeled_genes <- res_df$GeneSymbol[!grepl("^ENSMUSG", res_df$GeneSymbol)]
  
  # Select top 10 up and top 10 down genes by smallest padj for labeling (reduced to thin clustering)
  top_up <- res_df %>% 
    filter(Significant == "Up", GeneSymbol %in% labeled_genes) %>% 
    arrange(padj) %>% 
    slice_head(n = 10)
  
  top_down <- res_df %>% 
    filter(Significant == "Down", GeneSymbol %in% labeled_genes) %>% 
    arrange(padj) %>% 
    slice_head(n = 10)
  
  top_genes <- bind_rows(top_up, top_down) %>% 
    arrange(desc(negLog10Padj))  # Sort by descending significance for better label priority
  
  # Force label specific genes like Zbp1
  genes_to_label <- c("Zbp1")  # Add more if needed
  specific_labeled <- res_df %>% filter(GeneSymbol %in% genes_to_label & Significant != "NotSig")
  top_genes <- bind_rows(top_genes, specific_labeled) %>% distinct()
  
  # Build volcano plot
  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Padj)) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 1.5, position = position_jitter(width = 0.2, height = 7, seed = 123), na.rm = TRUE) +  # Reduced jitter; na.rm to suppress warnings
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "black")) +
    geom_vline(xintercept = c(-log2FC_thresh, log2FC_thresh), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = top_genes,
      aes(label = GeneSymbol),
      size = 3,  # Slightly increased size for visibility
      max.overlaps = 50,  # Increased to handle more labels
      segment.color = "grey50",  # Add connecting lines
      segment.size = 0.5,  # Thin line size
      segment.alpha = 0.8,  # Visibility for lines
      arrow = arrow(length = unit(0.1, "cm"))  # Directed arrows for clarity
    ) +
    theme_classic(base_size = 12) +  # Switch to classic theme for subtle grids
    labs(
      title = paste0("Volcano Plot: ", contrast_name),
      x = "log2 Fold Change",
      y = "-log10(adj. p-value)"
    ) +
    theme(legend.position = "right") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +  # Enforce y >=0; add extra space at top
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 300), clip = "off") +  # Keep increased y-cap; cap x; allow labels outside
    expand_limits(y = 0)  # Start y-axis at 0, auto scale top
  
  # Display plot
  print(p_volcano)
  
  # Save plots as PNG and PDF with increased height for tall y-axis
  ggsave(paste0("Volcano_", contrast_name, ".png"), plot = p_volcano, width = 8, height = 8, dpi = 300)
  ggsave(paste0("Volcano_", contrast_name, ".pdf"), plot = p_volcano, width = 11, height = 9)
  
  cat("Volcano plot generated for contrast:", contrast_name, "\n")
}
