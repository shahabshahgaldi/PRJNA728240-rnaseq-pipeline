#!/usr/bin/env bash
set -euo pipefail

# 01_setup_miniconda.sh
# Create the `rnaseq_env` conda environment for this RNA-seq pipeline.
# Usage: bash 01_setup_miniconda.sh
# Note: If conda is not installed this script will download the Miniconda installer
#       but you will need to run the installer interactively.

echo "== RNASEQ ENV SETUP =="
if ! command -v conda >/dev/null 2>&1; then
  echo "Conda not found on PATH. Downloading Miniconda installer to ~/Desktop..."
  wget -O ~/Desktop/Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  chmod +x ~/Desktop/Miniconda3-latest-Linux-x86_64.sh
  echo
  echo "Miniconda installer downloaded to ~/Desktop."
  echo "Run: bash ~/Desktop/Miniconda3-latest-Linux-x86_64.sh"
  echo "Follow interactive prompts, then reopen the terminal and re-run this script."
  exit 0
fi

echo "Conda found: $(conda --version)"

# Accept TOS for core channels if supported (non-fatal)
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r    2>/dev/null || true

echo "Creating conda environment 'rnaseq_env' (this may take some time)..."
conda create -y -n rnaseq_env \
  -c conda-forge -c bioconda --strict-channel-priority \
  python=3.10 \
  fastqc multiqc trimmomatic cutadapt sra-tools samtools \
  star hisat2 subread rsem salmon rseqc bwa bowtie2 pigz parallel \
  R r-base r-essentials \
  bioconductor-deseq2 bioconductor-tximport bioconductor-biomaRt bioconductor-apeglm || {
    echo "Conda env create failed. Inspect errors above."
    exit 1
  }

echo "Environment created. Activate with: conda activate rnaseq_env"
echo "Quick tests (run after activating the env):"
echo "  fastqc --version"
echo "  samtools --version"
echo "  R --version"
