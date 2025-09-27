# PRJNA728240 RNA-Seq Analysis Pipeline

This repository contains a complete RNA-Seq analysis pipeline for the unstranded single-end PRJNA728240 dataset.  
The pipeline processes raw SRA data into gene-level counts, followed by differential expression analysis and visualization of results.

## Repository Structure

- **01_sra_to_gene_counts/**  
  Contains scripts and resources to process raw SRA data into gene-level counts. This includes:
  - `01_setup_miniconda.sh` — script to set up a conda environment for reproducibility
  - `02_make_counts.sh` — script to convert SRA→FASTQ→trim→align→featureCounts
  - `mapping.csv` — file containing sample mapping information
  - `runs.txt` — list of SRR accession numbers for processing
  - `README.md` — detailed instructions for running this stage of the pipeline
  - `step-by-step_manual_instructions.docx` — comprehensive workflow guide
  - `results/` — contains raw and cleaned gene-level counts after processing

- **02_gene_counts_to_DE_and_plots/**  
  Contains the R pipeline for differential expression (DE) analysis and visualization. This includes:
  - `01_DE_analysis_pipeline.R` — DESeq2-based DE analysis pipeline
  - `raw_counts.rar` — compressed raw counts file
  - `README.md` — description of DE analysis steps
  - `results/` — contains DE results, normalized counts, QC plots, and volcano plots

## Getting Started

1. Clone the repository:
https://github.com/shahabshahgaldi/PRJNA728240-rnaseq-pipeline.git
2. Follow the instructions in `01_sra_to_gene_counts/README.md` to obtain gene counts from raw SRA files.
3. Follow the instructions in `02_gene_counts_to_DE_and_plots/README.md` to perform DE analysis and generate plots.

## Requirements

- Conda (Miniconda or Anaconda)
- Bash shell
- SRA Toolkit
- HISAT2, featureCounts
- R (≥4.1.0) with DESeq2, apeglm, and plotting packages

## Citation

If you use this pipeline for your research, please cite the PRJNA728240 dataset and reference this repository in your methods.
