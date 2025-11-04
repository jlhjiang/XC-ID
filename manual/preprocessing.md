# XC-ID Preprocessing
Date Created: 2025-10-31

Here, we provide simple preprocessing pipeline for generating suitable input for the XC-ID algorithm by Jiang and Gillis (2025).

You will need to install **STAR** via pip or bioconda (or a package manager of your choice).

---

## What should I do?
1. I have no variant information (VCF).

We will need to call variants prior to using STAR+WASP. We will align the FASTQ sequences using base STAR, and call variants using bcftools. Start with step 1.

**For human data**: this step can be bypassed by using common SNP positions found in databases such as 1000 Genomes or gnomAD.
**For mice data**: this step can be bypassed by using tools to compare genotypes between two strains and extracting the heterozygous positions.

2. I have variant information and/or a list of SNP positions to consider.

We will move directly to unbiased alignment using STAR+WASP.

---

## Optional Step 0: Build reference genome
It can be beneficial to align to a genome that has the Y and ALT contigs removed. Here is a safe way of doing so:

---

## Step 1: Initial Alignment
(reference STAR manual to adjust according to your library prep method)

```bash
STAR --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READ2_CSV" "$READ1_CSV" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --outSAMtype BAM SortedByCoordinate \
    --threads "$THREADS"
samtools index Aligned.sortedByCoord.out.bam
```

## Step 2: Call variants
You can use any variant callers such as bcftools, GATK, FreeBayes, etc... We will use bcftools for high sensitivity in single-cell data.

```bash
bcftools mpileup \
  -a AD \
  -f "$GENOME_FASTA" \
  -q 30 \
  -r X "Aligned.sortedByCoord.out.bam" \
  -Ou \
| bcftools call -mv \
  -Ou \
| bcftools norm \
  -m-any \
  --check-ref w \
  -f "$GENOME_FASTA" \
  -Ou \
| bcftools view \
-i 'GT="0/1" && QUAL>=20'
-Oz -o "variants.vcf.gz"
&& bcftools index "variants.vcf.gz"

gunzip variants.vcf.gz
```

You can then provide the resulting variants.vcf file to STAR+WASP and XC-ID.


## Step 3: STAR+WASP allele-specific alignment
### Existing variant information.

Assuming you have an existing VCF file with SNPs of interest, which can be common variants (i.e. population frequency >= 0.01) taken from databases such as gnomAD or 1000 Genomes or called de novo from FASTA/FASTQ or BAM files (variant calling).

The important parameters to include for STAR is
```bash
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG NH HI AS nM CB UB vW
```
Every other parameter can be adjusted according to the STAR manual.

Example:
```bash
STAR --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READ2_CSV" "$READ1_CSV" \
    --varVCFfile variants.vcf \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix WASP_ \
    --threads "$THREADS"\
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG NH HI AS nM CB UB vW 
samtools index WASP_Aligned.sortedByCoord.out.bam
```