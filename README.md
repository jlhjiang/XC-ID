# XC-ID
**X Chromosome Inactivation IDentifier**

XC-ID is a computational tool for identifying the per-cell active X lineage in single-cell RNA-seq (scRNA-seq) data from female mammals.  
It leverages allele-specific expression patterns at heterozygous X-linked SNPs to assign each cell to its active X chromosome (X₀ or X₁).

---

## Overview
XC-ID performs:
1. **Allele-specific read parsing** from STAR+WASP–aligned BAMs --> Cells X SNPs counts matrices
2. **Haplotype inference** using a simulated annealing–based optimization algorithm
3. **Bootstrap-based confidence estimation** for robust X lineage assignment  
Finally, XC-ID generates a results table that can be integrated with cell metadata in single-cell toolkits such as Scanpy or seurat.

Each sample/individual should be processed independently.

---

## Installation
Clone the repository and install via `pip`:

```bash
git clone https://github.com/jlhjiang/XC-ID.git
cd XC-ID
pip install -e .
```

## Usage

XC-ID can be used both in Python API or as a command-line tool. It can be conveniently implemented downstream of QC steps, see notebooks/tutorial.ipynb.
The most basic usage is:

        Python API
        >>> x = XCID(
        ...     vcf_files=['x.vcf.gz'],
        ...     bam_files=['Aligned.sortedByCoord.out.bam'],
                    )
        >>> res = x.run_all()
        >>> res.head()


```bash
xcid \
  --vcf PATH/TO/VCF \
  --bam PATH/TO/BAM \
  --out-results results.tsv
```

### Useful parameters

It is often helpful to increase the number of bootstraps (default 100) to 500 or 1000 for better sensitivity in "lower quality" datasets. Try this if you are getting lots of `unknown` cell labels in the results.

This can achieved by changing the `n_boot` parameter in Python API and `--n-boot` in CLI.

You can also finetune the number of jobs (default -1 or all available threads) and set a random seed.
Python API: `n_jobs`, `rand_seed`
CLI: `--n-jobs`, `--rand-seed`


## Basic preprocessing

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

An optional .txt file with a list of filtered (high-quality) cell barcodes can be provided to denoise. (see ./data/filtered_cells.txt for format).

## Variant calling

You can use any variant callers such as bcftools, GATK, FreeBayes... We provide a simple procedure using bcftools for high sensitivity.

### Starting with FASTA/FASTQ 
(reference STAR manual to adjust according to your data format)

```bash
STAR --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READ2_CSV" "$READ1_CSV" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --outSAMtype BAM SortedByCoordinate \
    --threads "$THREADS"
samtools index Aligned.sortedByCoord.out.bam
```

### Starting with BAM 
make sure that reference genomes are compatible with the one used in STAR+WASP

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

