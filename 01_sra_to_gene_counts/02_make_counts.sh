#!/usr/bin/env bash
set -euo pipefail

# 02_make_counts.sh
# Create raw gene-level counts for PRJNA728240.
# Usage: bash 02_make_counts.sh /path/to/PRJNA728240
# Default: ~/Desktop/PRJNA728240

PROJECT_DIR="${1:-$HOME/Desktop/PRJNA728240}"
echo "Project dir: $PROJECT_DIR"
mkdir -p "$PROJECT_DIR"
cd "$PROJECT_DIR"

# 1) Create directory structure
mkdir -p adapters alignments counts/featureCounts hisat2_index logs qc qc/qc_pre qc/qc_post \
         reference_genome fastq sra_cache trimmed sra_tmp

# 2) Create runs.txt (edit if needed)
if [ ! -f runs.txt ]; then
cat <<'EOF' > runs.txt
SRR14470855
SRR14470856
SRR14470857
SRR14470861
SRR14470862
SRR14470863
SRR14739702
SRR14739703
SRR14739704
EOF
fi
echo "runs.txt:"
cat runs.txt

# 3) Download SRA with prefetch into sra_cache
while read -r SRR; do
  echo "Prefetching $SRR ..."
  prefetch "$SRR" -O ./sra_cache || echo "prefetch failed for $SRR"
done < runs.txt

# 4) Validate SRA files (if vdb-validate present)
for d in ./sra_cache/*; do
  if [ -d "$d" ]; then
    for f in "$d"/*.sra; do
      [ -f "$f" ] || continue
      echo "Validating $f"
      vdb-validate "$f" || echo "Validation failed for $f"
    done
  fi
done

# 5) Convert SRA to FASTQ (single-end assumed) with fasterq-dump and compress with pigz
while read -r SRR; do
  mkdir -p fastq sra_tmp
  if [ -f "fastq/${SRR}.fastq.gz" ]; then
    echo "fastq/${SRR}.fastq.gz already exists, skipping."
    continue
  fi
  echo "Converting $SRR -> FASTQ..."
  fasterq-dump "$SRR" -O fastq --threads 4 --temp ./sra_tmp || { echo "fasterq-dump failed for $SRR"; continue; }
  if [ -f "fastq/${SRR}.fastq" ]; then
    pigz -p 4 "fastq/${SRR}.fastq" || true
  fi
  rm -rf ./sra_tmp/*
done < runs.txt

# 6) Rename FASTQ using mapping.csv (comma-delimited) if present
if [ -f mapping.csv ]; then
  echo "Renaming FASTQ files using mapping.csv (comma-delimited: SRR,sample_name)"
  while IFS=, read -r SRR sample; do
    [[ $SRR == "SRR" ]] && continue
    src="fastq/${SRR}.fastq.gz"
    dest="fastq/${sample}.fastq.gz"
    if [[ -f "$src" ]]; then
      mv "$src" "$dest"
      echo "Renamed: $src -> $dest"
    else
      echo "Warning: $src not found, skipping rename"
    fi
  done < mapping.csv
else
  echo "mapping.csv not found. Keep default SRR filenames or create mapping.csv (SRR,sample_name)."
fi

# 7) Pre-trim QC (FastQC + MultiQC)
mkdir -p qc/qc_pre
fastqc fastq/*.fastq.gz -o qc/qc_pre -t 4 || echo "FastQC failed (pre-trim)"
multiqc -o qc/qc_pre qc/qc_pre || echo "MultiQC failed (pre-trim)"

# 8) Trimming with Trimmomatic (single-end)
mkdir -p adapters trimmed logs
if [ ! -f adapters/BGI_adapters.fa ]; then
  cat > adapters/BGI_adapters.fa <<'AD'
>BGI_3prime
AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG
AD
fi

# Try to locate trimmomatic jar inside conda env
TRIMMOMATIC_JAR=$(find "$(conda info --base 2>/dev/null)" -type f -name "trimmomatic*.jar" 2>/dev/null | head -n 1 || true)
if [ -z "$TRIMMOMATIC_JAR" ]; then
  echo "Warning: trimmomatic jar not found automatically. Edit the script with correct path if needed."
fi

for fq in fastq/*.fastq.gz; do
  sample=$(basename "$fq" .fastq.gz)
  out="trimmed/${sample}_trim.fastq.gz"
  echo "Trimming $sample -> $out"
  if [ -n "$TRIMMOMATIC_JAR" ]; then
    java -jar "$TRIMMOMATIC_JAR" SE -phred33 -threads 4 \
      "$fq" "$out" \
      ILLUMINACLIP:adapters/BGI_adapters.fa:2:30:10:3:true \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30 \
      2> logs/trim_${sample}.log || echo "Trimmomatic failed for $sample"
  else
    echo "Skipping trimming because TRIMMOMATIC_JAR was not found."
  fi
done

# 9) Post-trim QC
mkdir -p qc/qc_post
fastqc trimmed/*.fastq.gz -o qc/qc_post -t 4 || echo "FastQC failed (post-trim)"
multiqc -o qc/qc_post qc/qc_post || echo "MultiQC failed (post-trim)"

# 10) Download reference genome + GTF (GRCm39) (if not present)
REFDIR="reference_genome"
GENOME_FA="$REFDIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF="$REFDIR/Mus_musculus.GRCm39.115.gtf"

if [ ! -f "$GENOME_FA" ]; then
  mkdir -p "$REFDIR" && cd "$REFDIR"
  wget -c https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
  wget -c https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz
  gunzip -k *.gz
  cd ..
fi

# 11) HISAT2 index build
mkdir -p hisat2_index
hisat2-build -p 4 "$GENOME_FA" hisat2_index/GRCm39_index > hisat2_index/build.log 2>&1 || echo "hisat2-build may have failed"

# 12) Align trimmed reads with HISAT2 -> sorted BAMs
mkdir -p alignments logs
for fq in trimmed/*_trim.fastq.gz; do
  sample=$(basename "$fq" _trim.fastq.gz)
  echo "Aligning $sample"
  hisat2 -p 4 --dta -x hisat2_index/GRCm39_index -U "$fq" \
    --summary-file logs/hisat2_${sample}.log -S alignments/${sample}_aligned.sam > logs/hisat2_${sample}_full.log 2>&1 || echo "hisat2 failed for $sample"
  samtools view -bS alignments/${sample}_aligned.sam | samtools sort -@ 4 -o alignments/${sample}_aligned.bam
  samtools index alignments/${sample}_aligned.bam
  rm -f alignments/${sample}_aligned.sam
done

# Optional: infer strandness with RSeQC (if installed)
if command -v infer_experiment.py >/dev/null 2>&1; then
  echo "Running RSeQC infer_experiment.py on first BAM (to check strandness)"
  firstbam=$(ls alignments/*_aligned.bam | head -n 1 || true)
  if [ -n "$firstbam" ]; then
    gtf2bed < "$GTF" > "$REFDIR/$(basename "$GTF" .gtf).bed" 2>/dev/null || true
    infer_experiment.py -i "$firstbam" -r "$REFDIR/$(basename "$GTF" .gtf).bed" -s 1000000 > logs/strandedness_check.txt || echo "infer_experiment.py failed"
  fi
fi

# 13) Quantify gene-level counts with featureCounts (unstranded: -s 0)
mkdir -p counts/featureCounts
featureCounts -T 4 -s 0 -a "$GTF" -o counts/featureCounts/raw_counts.txt -g gene_id -t exon alignments/*_aligned.bam 2> counts/featureCounts/featurecounts.log || echo "featureCounts may have failed"

# 14) Clean raw_counts.txt into simple table (geneID + sample counts)
head -n 2 counts/featureCounts/raw_counts.txt | tail -n 1 | cut -f1,7- | sed 's|alignments/||g' > counts/featureCounts/clean_raw_counts.txt
tail -n +3 counts/featureCounts/raw_counts.txt | cut -f1,7- >> counts/featureCounts/clean_raw_counts.txt
cat counts/featureCounts/clean_raw_counts.txt | tr '\t' ',' > counts/featureCounts/clean_raw_counts.csv

echo "Done. Check counts/featureCounts/clean_raw_counts.txt and .csv and logs/ for details."
