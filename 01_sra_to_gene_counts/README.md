````markdown
# 01_sra_to_gene_counts

This folder contains scripts and templates to produce **raw gene-level counts** (featureCounts output) from SRA for PRJNA728240.  
It covers the entire process: downloading SRA data, converting to FASTQ, trimming, aligning, and generating gene-level counts.

## Files
- `01_setup_miniconda.sh` — Creates the `rnaseq_env` conda environment (see comments inside).
- `02_make_counts.sh` — Main SRA → FASTQ → trim → align → featureCounts pipeline. Edit before running.
- `runs.txt` — List of SRR accessions (example provided).
- `mapping.csv` — Optional sample rename file (format: `SRR,sample_name`).
- `step-by-step_manual_instructions.docx` — Step-by-step tutorial with executable codes included (highlighted in yellow).
- `README.md` — This file.

## Important Notes
- These scripts are written for Linux and assume you run them on a machine with network access and enough disk space (~30+ GB for indices and FASTQ).
- `mapping.csv` should be comma-delimited with the first line as `SRR,sample_name`.
- The pipeline assumes single-end, unstranded reads (`featureCounts -s 0`). Change `-s` if your library is stranded.
- Don’t upload large FASTQ or alignment files to GitHub directly — use Git LFS or external storage (e.g., SRA, Zenodo, Dropbox).
- After cloning locally, make scripts executable if desired:
  ```bash
  chmod +x 01_setup_miniconda.sh 02_make_counts.sh
````

## Quick Usage (on a Linux machine)

1. Clone the repo and change to the folder:

   ```bash
   git clone <repo_url>
   cd PRJNA728240-rnaseq-pipeline/01_sra_to_gene_counts
   ```

2. Create and activate the environment:

   ```bash
   bash 01_setup_miniconda.sh
   conda activate rnaseq_env
   ```

3. Run the pipeline:

   ```bash
   bash 02_make_counts.sh
   ```

4. Check results in:

   ```
   counts/featureCounts/
   ```

## Extended Description of Key Files

* **`runs.txt`** — A simple text file listing SRR accession numbers to be downloaded and processed. Example:

  ```
  SRR14470855
  SRR14470856
  SRR14470857
  SRR14470861
  SRR14470862
  SRR14470863
  SRR14739702
  SRR14739703
  SRR14739704
  ```

* **`mapping.csv`** — Optional file used for renaming samples to meaningful names in downstream analysis. Format:

  ```
  SRR,sample_name
  SRR14739702,3T3L1-D2-rep3
  SRR14739703,3T3L1-D2-rep2
  ```

* **`step-by-step_manual_instructions.docx`** — Provides a detailed description of each pipeline step with highlighted executable commands for easy reference.

## Results

After running the pipeline, results will be available in:

```
counts/featureCounts/
```

Key result files:

* `raw_counts.txt` — Raw gene-level counts from featureCounts.
* `clean_raw_counts.txt` — Cleaned gene count table (tab-delimited).
* `clean_raw_counts.csv` — Cleaned gene count table (comma-delimited).

---

```
