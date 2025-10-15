#!/bin/bash
set -euo pipefail

# === PREREQS (one-time) ===
# conda create -n gatk4 -c bioconda -c conda-forge gatk4 bwa samtools fastqc -y
# conda activate gatk4

# === ROOT DIRS ===
BASE="$HOME/Desktop/demo"
SUP="$BASE/supporting_files/hg38"
VC_DIR="$BASE/VC"
READS="$VC_DIR/reads"
ALN="$VC_DIR/aligned_reads"
RES="$VC_DIR/results"
DATA="$VC_DIR/data"

mkdir -p "$SUP" "$READS" "$ALN" "$RES" "$DATA"

echo "== Download example reads =="
cd "$READS"
# If FTP is blocked on your network, use 'prefetch/fastq-dump' instead (SRA toolkit).
wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

echo "== Prep: download MATCHING hg38 reference & known-sites =="
cd "$SUP"
# Use the GATK bundle hg38 (no 'chr' prefixes), consistent with dbSNP v0 files.
# Reference (bgzipped) + dict + fai:
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

# Known sites (dbSNP); index is provided:
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

REF="$SUP/Homo_sapiens_assembly38.fasta"
KNOWN="$SUP/Homo_sapiens_assembly38.dbsnp138.vcf"

echo "== Step 1: FastQC =="
fastqc "$READS/SRR062634_1.filt.fastq.gz" -o "$READS"
fastqc "$READS/SRR062634_2.filt.fastq.gz" -o "$READS"

echo "== Step 2: BWA index (one-time) & align =="
# BWA index is persistent; guarded to save time on reruns.
[ -f "$REF.bwt" ] || bwa index "$REF"

bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" \
  "$REF" \
  "$READS/SRR062634_1.filt.fastq.gz" "$READS/SRR062634_2.filt.fastq.gz" \
  | samtools view -Sb - \
  | samtools sort -@ 4 -o "$ALN/SRR062634.sorted.bam"

samtools index "$ALN/SRR062634.sorted.bam"

echo "== Step 3: Mark duplicates (sorted output) =="
gatk MarkDuplicatesSpark \
  -I "$ALN/SRR062634.sorted.bam" \
  -O "$ALN/SRR062634.sorted.dedup.bam" \
  --create-output-bam-index true

echo "== Step 4: Base Quality Score Recalibration (BQSR) =="
gatk BaseRecalibrator \
  -R "$REF" \
  -I "$ALN/SRR062634.sorted.dedup.bam" \
  --known-sites "$KNOWN" \
  -O "$DATA/recal_data.table"

gatk ApplyBQSR \
  -R "$REF" \
  -I "$ALN/SRR062634.sorted.dedup.bam" \
  --bqsr-recal-file "$DATA/recal_data.table" \
  -O "$ALN/SRR062634.sorted.dedup.bqsr.bam"

samtools index "$ALN/SRR062634.sorted.dedup.bqsr.bam"

echo "== Step 5: Alignment & insert-size metrics =="
gatk CollectAlignmentSummaryMetrics \
  -R "$REF" \
  -I "$ALN/SRR062634.sorted.dedup.bqsr.bam" \
  -O "$ALN/alignment_metrics.txt"

gatk CollectInsertSizeMetrics \
  -I "$ALN/SRR062634.sorted.dedup.bqsr.bam" \
  -O "$ALN/insert_size_metrics.txt" \
  -H "$ALN/insert_size_histogram.pdf"

echo "== Step 6: Variant calling (HaplotypeCaller) =="
# For multi-sample best practices, use -ERC GVCF and combine later.
gatk HaplotypeCaller \
  -R "$REF" \
  -I "$ALN/SRR062634.sorted.dedup.bqsr.bam" \
  -O "$RES/raw_variants.vcf"

echo "== Step 7: Split SNPs and INDELs =="
gatk SelectVariants -R "$REF" -V "$RES/raw_variants.vcf" --select-type-to-include SNP   -O "$RES/raw_snps.vcf"
gatk SelectVariants -R "$REF" -V "$RES/raw_variants.vcf" --select-type-to-include INDEL -O "$RES/raw_indels.vcf"

echo "All done. Outputs in: $RES"