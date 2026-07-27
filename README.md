# XC-ID
**XC-ID: De novo identification of the active X chromosome in single-cell RNA-seq**
[Jennifer Jiang, Jesse Gillis (2025)](https://doi.org/10.64898/2025.12.05.692676)

## Overview
**XC-ID** is a computational tool for identifying the *active X chromosome lineage* (maternal or paternal haplotype) in **single-cell RNA-seq (scRNA-seq)** data from female mammals.  
It infers per-cell X-chromosome activation states using allele-specific expression patterns at heterozygous X-linked SNPs.

XC-ID is available both as:
- a **Python API** (for integration with Jupyter or Scanpy workflows), and  
- a **command-line interface** for batch processing.

---

## Features

XC-ID performs the following key steps:

1. **Allele-specific read parsing** – extracts per-cell, per-SNP counts from STAR+WASP–aligned BAM files.  
2. **Haplotype inference and phasing** – assigns each cell’s active X haplotype using a simulated annealing optimization.  
3. **Bootstrap-based confidence estimation** – quantifies assignment robustness through iterative resampling.

The output table integrates with standard single-cell analysis frameworks (e.g., **Scanpy**, **Seurat**) metadata.  
XC-ID has been tested across 10X Genomics (v2, v3, GEM-X) and SmartSeq platforms, and in both human and non-human mammalian datasets.

---

## Installation
Clone the repository and install via `pip`:

```bash
git clone https://github.com/jlhjiang/XC-ID.git
cd XC-ID
pip install -e .
```

---

## Usage

### 1. Preprocessing

XC-ID requires:
- **BAM files** aligned with STAR using WASP tagging (`--waspOutputMode SAMtag`)  
- **VCF files** containing heterozygous SNP positions for phasing

Preprocessing pipelines can vary depending on your data type. A basic setup only requires a single STAR+WASP alignment step.  
Detailed preprocessing instructions are available in the [documentation](/manual/preprocessing.md).

### 2. Running XC-ID

You can run XC-ID via the **Python API** or **CLI**. The most basic use example is:

#### **Python API**
```python
from xcid import XCID

xc = XCID(
    vcf=["/path/to/variants.vcf"],
    bam=["/path/to/aligned.bam"],
)
results = xc.run_all()
results.head()
```

See the Jupyter Notebook tutorial [here](/manual/tutorial.ipynb)

#### **Command-line**
```bash
xcid --vcf /path/to/variants.vcf --bam /path/to/aligned.bam --out-results /path/to/results.tsv
```

> You may pass multiple BAM or VCF files (space-separated) for the same sample, e.g. `--bam rep1.bam rep2.bam`. Multiple BAM files are read in parallel.

### 3. Useful parameters

Python API arguments now match the CLI flag names directly (e.g. `chrom` <-> `--chrom`, `seed` <-> `--seed`), so options translate one-to-one between the two interfaces.

It is often helpful to increase the number of bootstraps (default 100) to 500 or 1000 for better sensitivity in "lower quality" datasets. Try this if you are getting lots of `unknown` cell labels in the results.

You can also finetune the haplotype counts matrices by providing a filtered cells file, which should be a .txt or .tsv file with one barcode per line. This can be the filtered `features.tsv` file from the STAR+WASP output.

| Parameter | API argument | CLI flag | Description |
|------------|---------------|-----------|--------------|
| VCF file(s) | `vcf` | `--vcf` | One or more VCF files with heterozygous SNP positions |
| BAM file(s) | `bam` | `--bam` | One or more STAR+WASP aligned BAM files |
| Chromosome | `chrom` | `--chrom` | Chromosome name as it appears in the BAM/VCF (default: `chrX`) |
| Filtered cells | `cells` | `--cells` | Allowed cell barcodes for constructing the counts matrices |
| Number of bootstraps | `n_boot` | `--n-boot` | Increase to 500–1000 for low-quality datasets |
| Number of parallel jobs | `n_jobs` | `--n-jobs` | Default = -1 (use all available cores); also parallelizes reading multiple BAM files |
| Random seed | `seed` | `--seed` | Ensures reproducibility |
| Min UMIs per SNP | `min_per_snp` | `--min-per-snp` | Minimum UMIs covering a SNP to keep it for phasing (default: 5) |
| Min counts per cell | `min_per_cell` | `--min-per-cell` | Minimum SNP counts in a cell to keep it (default: 2) |
| Min minor allele frequency | `min_maf` | `--min-maf` | Minimum minor allele frequency to keep a SNP (default: 0.1) |

See full parameter list with:
```bash
xcid --help
```

### Output

XC-ID outputs a results table with columns:
1. **cell_id:** the cell barcodes
2. **score:** the haplotype ratio score
3. **p_value:** p-values from a binomial test on bootstrapped scores
4. **p_adj:** Benjamini-Hochberg adjusted p-values
5. **XCI_status:** the final X chromosome calls from the algorithm (`X0`, `X1`, or `unknown`)

Note that X₀ or X₁ are direction agnostic. For genotyping, refer to the full manual for determining the exact REF/ALT SNP assignment for each X haplotype.

---

## Full manual

Reference to the [manual](/manual/preprocessing.md) for more detailed explanation, including how to extract the exact REF/ALT genotypes, and get putative escape genes/proportions.

