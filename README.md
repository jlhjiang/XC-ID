# XC-ID
**X Chromosome inactivation IDentifier**

Jennifer Jiang, Jesse Gillis (2025)

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
    vcf_files=["/path/to/variants.vcf"],
    bam_files=["/path/to/aligned.bam"]
)
results = x.run_all()
results.head()
```

#### **Command-line**
```bash
xcid   --vcf /path/to/variants.vcf   --bam /path/to/aligned.bam   --out-results /path/to/results.tsv
```

> You may provide multiple BAM or VCF files (comma-separated) for the same sample.

### 3. Useful parameters

It is often helpful to increase the number of bootstraps (default 100) to 500 or 1000 for better sensitivity in "lower quality" datasets. Try this if you are getting lots of `unknown` cell labels in the results.

You can also finetune the haplotype counts matrices by providing a filtered cells file, which should be a .txt or .tsv file with one barcode per line. This can be the filtered `features.tsv` file from the STAR+WASP output.

| Parameter | API argument | CLI flag | Description |
|------------|---------------|-----------|--------------|
| Number of bootstraps | `n_boot` | `--n-boot` | Increase to 500–1000 for low-quality datasets |
| Number of parallel jobs | `n_jobs` | `--n-jobs` | Default = -1 (use all available cores) |
| Random seed | `rand_seed` | `--rand-seed` | Ensures reproducibility |
| Filtered cells | `cells_file` | `--cells` | Allowed cell barcodes for constructing the counts matrices |

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
5. **X_status: the final X chromosome calls from the algorithm**

Note that X₀ or X₁ are direction agnostic. For genotyping, refer to the full manual for determining the exact REF/ALT SNP assignment for each X haplotype.

---

## Full manual

Reference to the [manual] for more detailed explanation, including how to extract the exact REF/ALT genotypes, and get putative escape genes/proportions.

