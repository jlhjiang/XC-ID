# XC-ID Preprocessing
Date Created: 2025-10-31

XC-ID needs two inputs per sample:
1. A **VCF** of heterozygous X-linked SNP positions to phase.
2. A **BAM** aligned with STAR using WASP allele-specific tagging (`--waspOutputMode SAMtag`).

This page walks through generating both. You will need **STAR** and **samtools** installed; if you need to call variants yourself, you will also need **bcftools**.

---

## Which path do I need?

Start here to figure out which steps below actually apply to you.

| Your situation | What to do |
|---|---|
| I don't have a VCF of SNP positions yet, and I'm working with **human** data | Skip variant calling — use the provided [gnomad_for_XCID.vcf](/data/gnomad_for_XCID.vcf) (common X chromosome variants, population AF ≥ 0.01), then go straight to **Step 3**. |
| I don't have a VCF, and I'm working with **mouse** (or other multi-strain) data | Skip variant calling — derive heterozygous positions by comparing genotypes between the two parental strains, then go to **Step 3**. |
| I don't have a VCF and want to call variants de novo from my own reads | Do **Step 0** (optional) → **Step 1** (alignment) → **Step 2** (variant calling) → **Step 3**. |
| I already have a VCF of SNP positions (from any source) | Go straight to **Step 3**. |

In short: most users only need **Step 3**. Steps 0–2 are only for building a VCF from scratch when you don't already have one and can't use a public/reference variant set.

---

## Step 0 (optional): Build a cleaned reference genome

Only needed if you're aligning from scratch (Steps 0–2 path above). It can be beneficial to align to a genome with the Y chromosome and ALT/random contigs removed, since reads that would otherwise map ambiguously to X/Y homologous regions can confound allele-specific counting:

```bash
awk '/^>/{keep = ($0 !~ /_alt|_random|chrUn|chrY/)} keep' \
    GRCh38.primary_assembly.genome.fa > GRCh38_noY_noALT.fa

grep -v -E "_alt|_random|chrUn|chrY" gencode.v47.annotation.gtf > gencode.v47_noY_noALT.gtf
```

---

## Step 1: Initial alignment (for de novo variant calling only)

A plain STAR alignment, used only to generate a BAM for variant calling in Step 2. Adjust parameters to your library prep method — see the [STAR manual](https://github.com/alexdobin/STAR) for details.

```bash
STAR --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READ2_CSV" "$READ1_CSV" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --outSAMtype BAM SortedByCoordinate \
    --threads "$THREADS"
samtools index Aligned.sortedByCoord.out.bam
```

## Step 2: Call variants (for de novo variant calling only)

Any variant caller works (bcftools, GATK, FreeBayes, ...). We use **bcftools** below for its sensitivity on single-cell data. This calls heterozygous SNPs on chromosome X and writes them to `variants.vcf`:

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
    -i 'GT="0/1" && QUAL>=20' \
    -Oz -o "variants.vcf.gz"

bcftools index "variants.vcf.gz"
gunzip variants.vcf.gz
```

The resulting `variants.vcf` is the SNP list to pass into Step 3 (as `--varVCFfile`) and, later, into XC-ID (as `--vcf`).

---

## Step 3: STAR+WASP allele-specific alignment (required for everyone)

This is the step that produces the BAM XC-ID actually reads. You need a VCF of SNP positions to phase — either the one you built in Step 2, the provided [gnomad_for_XCID.vcf](/data/gnomad_for_XCID.vcf) for human data, or one derived from strain comparison for mouse data.

The two STAR parameters that **must** be set for XC-ID to work are:
```bash
--waspOutputMode SAMtag \
--outSAMattributes vA vG NH HI AS nM CB UB vW
```
These add the `vA`/`vG`/`vW` read tags XC-ID uses for allele-specific counting, plus the `CB`/`UB` cell/UMI barcode tags. Every other STAR parameter can be adjusted freely to match your library prep — see the [STAR manual](https://github.com/alexdobin/STAR).

Full example:
```bash
STAR --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READ2_CSV" "$READ1_CSV" \
    --varVCFfile variants.vcf \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix WASP_ \
    --threads "$THREADS" \
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG NH HI AS nM CB UB vW
samtools index WASP_Aligned.sortedByCoord.out.bam
```

The resulting `WASP_Aligned.sortedByCoord.out.bam` and your `variants.vcf` are exactly the two inputs XC-ID needs:

```bash
xcid --vcf variants.vcf --bam WASP_Aligned.sortedByCoord.out.bam --out-results results.tsv
```
