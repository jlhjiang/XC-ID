import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def load_vcf(vcf, blacklist=frozenset(), chrom='chrX'):
    """Loads the list of VCF files and returns a dictionary of the format:
    {pos: {'rs':ID, 'ref':ref, 'alt':alt, 'type':'snp'}}."""
    nuc = [
        'A', 'C', 'G', 'T'
        ]
    snp_info = {}

    for path in vcf:
        with open(path, 'rt') as fh:
            for line in fh:
                if line[0] == '#':
                    continue

                c, pos, rs, ref, alt, *_ = line.split('\t', 5)
                if c != chrom:
                    continue

                pos = int(pos)
                if pos in blacklist:
                    continue

                var_type = 'snp' if (ref in nuc and alt in nuc) else 'indel'
                if var_type == 'snp':
                    snp_info[pos] = {'rs': rs, 'ref': ref, 'alt': alt,
                                     'type': 'snp'}

    n = len(snp_info)
    print(f'{n} SNPs loaded from {chrom}')
    if not n:
        sys.exit('No SNPs found, terminating.')
    return snp_info


def _read_bam_single(bam, snp_info, filtered_cells, chrom):
    '''Read a single BAM file and extract UMI and allele information for SNP positions.
    Returns (umis, cell2allele, bad_cells) for this file only.'''

    bad_cells = 0
    umis        = {}
    cell2allele = {}

    with pysam.AlignmentFile(bam, "rb") as samfile:
        # fetch reads on chromosome (or region)
        for read in samfile.fetch(chrom):
            # Check for the presence of vG and vW tags.
            try:
                vG = read.get_tag("vG")
                vW = read.get_tag("vW")
            except KeyError:
                continue

            # vW is 1 if passed WASP filtering
            if vW != 1:
                continue

            # get vA (allele information): 0 for ref, 1 for alt
            try:
                vA = read.get_tag("vA")
            except KeyError:
                continue
            # fall back
            if not vG or not vA:
                continue

            # assume vG and vA are lists (since they are B tags)
            try:
                vGpos = list(vG)
                vAvals = list(vA)
            except (TypeError, ValueError, AttributeError):
                continue

            # Parse cell id from CB or RG tags.
            cell_id = None
            umi_id  = None
            if read.has_tag("CB"):
                cell_id = read.get_tag("CB")
            elif read.has_tag("RG"):
                cell_id = read.get_tag("RG")
            if read.has_tag("UB"):
                umi_id = read.get_tag("UB")

            if not cell_id:
                continue  # skip read if no cell id

            if filtered_cells is not None:
                if cell_id not in filtered_cells: # filter out cells not in the filtered list
                    bad_cells += 1
                    continue

            if not umi_id:
                # Fallback: construct a UMI from available read data
                umi_id = f"{cell_id}_{read.reference_start}_{read.flag}_{read.cigarstring}_{read.template_length}"

            # For each variant position (and corresponding allele) stored in vG/vA:
            for g, a in zip(vGpos, vAvals):
                try:
                    g = int(g)
                except ValueError:
                    continue
                # vG positions are 0-based, add 1 to convert to 1-based
                posx = g + 1

                if posx not in snp_info:
                    continue

                try:
                    allele = 'ref' if int(a)==1 else 'alt'
                except ValueError:
                    continue
                if int(a) not in (1, 2):
                    continue

                umikey = (umi_id, cell_id)
                if posx not in umis:
                    umis[posx] = {'ref': {}, 'alt': {}}
                if umikey not in umis[posx][allele]:
                    umis[posx][allele][umikey] = 0
                umis[posx][allele][umikey] += 1

                if cell_id not in cell2allele:
                    cell2allele[cell_id] = {}
                if posx not in cell2allele[cell_id]:
                    cell2allele[cell_id][posx] = {}
                cell2allele[cell_id][posx][umi_id] = allele

    return umis, cell2allele, bad_cells


def _merge_umis(umis, other):
    '''Merge a per-file umis dict into the running total, summing per-UMI counts.'''
    for pos, alleles in other.items():
        dest = umis.setdefault(pos, {'ref': {}, 'alt': {}})
        for allele in ('ref', 'alt'):
            dest_allele = dest[allele]
            for umikey, cnt in alleles[allele].items():
                dest_allele[umikey] = dest_allele.get(umikey, 0) + cnt


def _merge_cell2allele(cell2allele, other):
    '''Merge a per-file cell2allele dict into the running total.'''
    for cell, posmap in other.items():
        dest = cell2allele.setdefault(cell, {})
        for pos, umimap in posmap.items():
            dest.setdefault(pos, {}).update(umimap)


def read_bam(bams = [],
             snp_info = {},
             filtered_cells = None,
             chrom = 'chrX',
             n_jobs = 1):
    '''Read BAM files and extract UMI and allele information for SNP positions.
    Returns a dictionary of UMIs and a mapping of cells to alleles.
    When multiple BAM files are provided, they are read in parallel (see n_jobs).'''

    if filtered_cells is not None:
        print(f"Number of filtered cells provided: {len(filtered_cells)}")
    else:
        filtered_cells = None
        print("No filtered cells file provided, all cells will be processed.")

    if len(bams) > 1:
        partials = Parallel(n_jobs=n_jobs, backend='loky')(
            delayed(_read_bam_single)(bam, snp_info, filtered_cells, chrom)
            for bam in bams
        )
    else:
        partials = [_read_bam_single(bam, snp_info, filtered_cells, chrom) for bam in bams]

    umis = {}
    cell2allele = {}
    bad_cells = 0
    for file_umis, file_cell2allele, file_bad_cells in partials:
        _merge_umis(umis, file_umis)
        _merge_cell2allele(cell2allele, file_cell2allele)
        bad_cells += file_bad_cells

    print(f"Loaded {len(umis)} SNP positions with coverage.")
    print(f"Loaded {len(cell2allele)} cells with informative reads.")
    if not umis:
        sys.exit("No SNPs after reading, check the chromosome name.")
    if not cell2allele:
        sys.exit("No cells found, check the chromosome name.")
    if filtered_cells is not None:
        print(f"Filtered out {bad_cells} cells not in the filtered list.")

    return umis, cell2allele


def make_results(cell_ids, score, score_array):
    '''Create a results DataFrame with cell IDs, scores, p-values, adjusted p-values, and XCI status.'''

    def binom_pval(score, score_array):
        '''Calculate binomial p-values for each cell.'''
        p_value = []
        for i, s in enumerate(score):
            if s > 0:
                n = (score_array[:, i] > 0).sum()
            elif s < 0:
                n = (score_array[:, i] < 0).sum()
            elif s == 0:
                n = 0
            p_value.append(binomtest(n, score_array.shape[0], 0.5, alternative='greater').pvalue)

        p_value = np.array(p_value)
        _, p_adj, _, _ = multipletests(p_value, alpha=0.05, method='fdr_bh')
        return p_value, p_adj

    p_value, p_adj = binom_pval(score, score_array)
    conf = (p_adj < 0.05) & (score!=0)

    def _label(i: int, s: float) -> str:
        if not conf[i]:
            return "unknown"
        return "X0" if s > 0 else "X1"

    results = pd.DataFrame({
        "cell_id": cell_ids,
        "score": score,
        "p_value": p_value,
        "p_adj": p_adj,
        "XCI_status": [_label(i, s) for i, s in enumerate(score)],
    })

    return results


def load_blacklist(path) -> set:
    '''Load a set of blacklisted SNP positions from a file (one 1-based position per line).'''
    blacklist = set()
    if not path:
        return blacklist
    with open(path, "rt") as fh:
        for line in fh:
            s = line.strip()
            if s:
                blacklist.add(int(s))
    return blacklist


def write_results_table(results: pd.DataFrame, path) -> None:
    '''Write the results table to a TSV file, creating parent directories as needed.'''
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(path, sep="\t", index=False)


def write_escape_table(escape_df: pd.DataFrame, path) -> None:
    '''Write the escape candidate table to a TSV file, creating parent directories as needed.'''
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    escape_df.to_csv(path, sep="\t", index=False)


def write_haplotype_matrices(hap0_counts: pd.DataFrame, hap1_counts: pd.DataFrame, out_dir) -> None:
    '''Write haplotype count matrices to hap0_counts.tsv / hap1_counts.tsv in out_dir.'''
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hap0_counts.to_csv(out_dir / "hap0_counts.tsv", sep="\t")
    hap1_counts.to_csv(out_dir / "hap1_counts.tsv", sep="\t")
