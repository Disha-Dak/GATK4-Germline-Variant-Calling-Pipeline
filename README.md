# GATK4-Germline-Variant-Calling-Pipeline

This repository contains a **bash script** (`gatk_germline.sh`) implementing the **GATK4 Best Practices workflow** for **germline short variant discovery (SNPs & INDELs)** in human whole-genome sequencing (WGS) data.

The script automates every major step — from downloading reference files to variant calling — following Broad Institute’s recommendations.

---

## Overview

**Input:**
- Paired-end FASTQ files (2 × 100 bp) from 1000 Genomes Project (HG00096 example)

**Output:**
- Quality reports (`FastQC`)
- Deduplicated, recalibrated BAM file
- Raw variant calls (`raw_variants.vcf`)
- Separate SNP and INDEL VCFs

**Reference genome:**
- Homo sapiens assembly GRCh38 (Broad GATK bundle)



##  Workflow Steps

| Step | Tool | Description |
|:--|:--|:--|
| 1 | `fastqc` | Quality control on raw reads |
| 2 | `bwa mem` | Alignment to hg38 reference |
| 3 | `samtools` | Convert SAM → BAM and sort |
| 4 | `gatk MarkDuplicatesSpark` | Remove PCR duplicates |
| 5 | `gatk BaseRecalibrator` & `ApplyBQSR` | Base quality score recalibration using dbSNP |
| 6 | `gatk CollectMetrics` | Collect alignment & insert-size statistics |
| 7 | `gatk HaplotypeCaller` | Call SNPs & INDELs |
| 8 | `gatk SelectVariants` | Split SNPs and INDELs into separate VCFs |

---

## 
Requirements

Install dependencies (preferably using Conda):

```bash
conda create -n gatk4 -c bioconda -c conda-forge gatk4 bwa samtools fastqc -y
conda activate gatk4
