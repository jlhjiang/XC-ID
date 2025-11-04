# XC-ID
**X Chromosome Inactivation IDentifier**

Jennifer Jiang, Jesse Gillis (2025)

XC-ID is a computational tool for identifying the per-cell active X lineage in single-cell RNA-seq (scRNA-seq) data from female mammals.  
It leverages allele-specific expression patterns at heterozygous X-linked SNPs to assign each cell to its active X chromosome (X₀ or X₁).

XC-ID can be used in both Python API and as a command-line tool

---

## Overview
XC-ID performs:
1. **Allele-specific read parsing** from STAR+WASP–aligned BAMs --> cells X SNPs counts matrices
2. **Haplotype inference** using a simulated annealing–based optimization algorithm
3. **Bootstrap-based confidence estimation** for robust X lineage assignment  
Finally, XC-ID generates a results table that can be readily integrated with cell metadata in single-cell toolkits such as Scanpy or seurat.

Currently, XC-ID has been tested on 10X v2, v3, GEM-X, and SmartSeq technology. It has also been successfully implemented in non-human/mice mammalian species.

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

### Preprocessing

XC-ID required the input of BAMs with WASP tags and the VCF files with variant positions to consider for the phasing algorithm.
Preprocessing steps can vary depending on what data you have on hand: at the simplest, it only requires one STAR+WASP alignment step. 

We provide simple preprocessing guidelines at [].

### Running the XC-ID

XC-ID can be used both in Python API or as a command-line tool. It can be conveniently implemented downstream of QC steps, see notebooks/tutorial.ipynb. Each sample should be processed independently.

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
  --out-results PATH/TO/RESULTS.tsv
```

### Useful parameters

It is often helpful to increase the number of bootstraps (default 100) to 500 or 1000 for better sensitivity in "lower quality" datasets. Try this if you are getting lots of `unknown` cell labels in the results. This can achieved by changing the `n_boot` parameter in Python API and `--n-boot` in CLI.

You can also finetune the number of jobs (default -1 or all available threads) and set a random seed.
Python API: `n_jobs`, `rand_seed`
CLI: `--n-jobs`, `--rand-seed`

---

## Full manual
### Genotyping & other uses

Reference to the manual here for more detailed explanation, including how to get the exact REF/ALT genotypes, and get putative escape genes/proportions.

