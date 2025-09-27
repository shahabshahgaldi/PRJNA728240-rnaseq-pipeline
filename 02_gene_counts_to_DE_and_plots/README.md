# 02_gene_counts_to_DE_and_plots

This folder contains scripts and templates for performing **differential expression (DE) analysis** and generating plots from gene-level counts obtained in the previous step (`01_sra_to_gene_counts`).  
The analysis is performed using DESeq2 and includes quality control (QC), normalization, DE testing, and volcano plot generation.

## Files

- `01_DE_analysis_pipeline.R` — R script implementing the full DE analysis pipeline (design matrix, QC, filtering, normalization, DE testing, volcano plots, mapping gene symbols).
- `raw_counts.rar` — Compressed raw counts output from the previous step (`featureCounts` results).
- `README.md` — This file.
- `LICENSE` — licensing information for this folder.
- `results/` — Folder containing analysis outputs.

## Important Notes

- The pipeline assumes raw counts produced from **featureCounts** with unstranded single-end reads.
- Ensure R and required packages are installed before running the pipeline.
- QC and DE results will be saved in the `results/` folder.
- Large intermediate files are not included; only essential results are stored.

## Results

After running the analysis pipeline, the following key results will be available in `results/`:

- `DE_results_d2_vs_d0_withSymbols.csv` — Differential expression results (d2 vs d0) with mapped gene symbols.
- `DE_results_d6_vs_d0_withSymbols.csv` — Differential expression results (d6 vs d0) with mapped gene symbols.
- `filter_threshold_summary.csv` — Summary of filtering thresholds and number of genes retained.
- `normalized_counts.csv` — Normalized gene counts for all samples.
- `Volcano_d2_vs_d0.pdf` — Volcano plot for d2 vs d0.
- `Volcano_d6_vs_d0.pdf` — Volcano plot for d6 vs d0.

Additional results and plots may be stored in this folder depending on analysis settings.

## Extended Description of Key Files

- **`01_DE_analysis_pipeline.R`** — Main R script implementing the DE analysis workflow. Includes QC, filtering, normalization, DE testing, gene annotation, and volcano plots.
- **`raw_counts.rar`** — Compressed raw gene count file from the previous step (featureCounts output), which will be used to analyze DE using RStudio.
- **`results/`** — Folder containing final DE analysis results and plots for downstream interpretation.
