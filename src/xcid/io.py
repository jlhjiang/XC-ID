import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed, effective_n_jobs
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


def _open_bam(bam, threads=1):
    '''Open a BAM, silencing htslib's "index file is older than the data file" warning, which
    would otherwise print once per chunk. Verbosity is only lowered while the index loads, so
    warnings raised while reading still show.'''
    old_verbosity = pysam.set_verbosity(0)
    try:
        return pysam.AlignmentFile(bam, "rb", threads=threads)
    finally:
        pysam.set_verbosity(old_verbosity)


# size of the chromosome chunks read in parallel by read_bam
BAM_CHUNK_SIZE = 5_000_000
# SNP lookups past the chunk end within this distance (e.g. spliced reads) use a set;
# anything further falls back to a binary search over all positions
SNP_WINDOW = 1_000_000


def _read_bam_single(bam, snp_pos, filtered_cells, chrom, start, end, threads=1):
    '''Read one region of a single BAM file and extract UMI and allele information for SNP positions.
    snp_pos is a sorted array of 1-based SNP positions. Only reads starting in [start, end) are
    kept, so adjacent chunks never count a read twice.
    Returns (umis, cell2allele, bad_cells) for this region only.'''

    bad_cells = 0
    umis        = {}
    cell2allele = {}

    # reads starting at 0-based `start` only cover 1-based positions > start
    lo, hi = np.searchsorted(snp_pos, [start + 1, end + SNP_WINDOW + 1])
    local_snps = set(snp_pos[lo:hi].tolist())
    window_end = end + SNP_WINDOW

    def is_snp(posx):
        if posx <= window_end:
            return posx in local_snps
        i = np.searchsorted(snp_pos, posx)
        return i < len(snp_pos) and snp_pos[i] == posx

    with _open_bam(bam, threads=threads) as samfile:
        # fetch reads on chromosome (or region)
        for read in samfile.fetch(chrom, start, end):
            # reads spanning the chunk start belong to the previous chunk
            if read.reference_start < start:
                continue

            # most reads carry no WASP tags: reject them with has_tag rather than
            # a raised KeyError, which is much slower on this hot path
            if not read.has_tag("vW"):
                continue
            # vW is 1 if passed WASP filtering
            if read.get_tag("vW") != 1:
                continue
            if not (read.has_tag("vG") and read.has_tag("vA")):
                continue
            vG = read.get_tag("vG")
            # get vA (allele information): 1 for ref, 2 for alt
            vA = read.get_tag("vA")
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
                # Fallback: construct a UMI from the fragment's position. Both mates of a proper
                # pair get the same key (leftmost start + fragment length), so a fragment is counted
                # once; the flag is left out so duplicate-marked (0x400) copies collapse too.
                if read.is_proper_pair:
                    frag_start = min(read.reference_start, read.next_reference_start)
                    umi_id = f"{cell_id}_{frag_start}_{abs(read.template_length)}"
                else:
                    umi_id = f"{cell_id}_{read.reference_start}_{read.is_reverse}_{read.cigarstring}"

            # For each variant position (and corresponding allele) stored in vG/vA:
            for g, a in zip(vGpos, vAvals):
                try:
                    g = int(g)
                except ValueError:
                    continue
                # vG positions are 0-based, add 1 to convert to 1-based
                posx = g + 1

                if not is_snp(posx):
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


def _bam_chunks(bams, chrom, chunk_size=BAM_CHUNK_SIZE):
    '''Split `chrom` in each BAM into (bam, start, end) regions, in file then genomic order.'''
    chunks = []
    for bam in bams:
        with _open_bam(bam) as samfile:
            length = samfile.get_reference_length(chrom)
        for start in range(0, length, chunk_size):
            chunks.append((bam, start, min(start + chunk_size, length)))
    return chunks


def read_bam(bams = [],
             snp_info = {},
             filtered_cells = None,
             chrom = 'chrX',
             n_jobs = 1):
    '''Read BAM files and extract UMI and allele information for SNP positions.
    Returns a dictionary of UMIs and a mapping of cells to alleles.
    The chromosome in each BAM is split into chunks that are read in parallel. n_jobs is the
    total CPU budget: it is shared between parallel workers and per-worker BAM decompression
    threads, so the run never uses more than n_jobs cores.'''

    if filtered_cells is not None:
        filtered_cells = set(filtered_cells) # O(1) membership checks per read
        print(f"Number of filtered cells provided: {len(filtered_cells)}")
    else:
        filtered_cells = None
        print("No filtered cells file provided, all cells will be processed.")

    chunks = _bam_chunks(bams, chrom)
    n_cpus = effective_n_jobs(n_jobs)
    n_workers = min(n_cpus, len(chunks))
    threads = max(1, n_cpus // n_workers)
    print(f"Reading {len(bams)} BAM file(s) in {len(chunks)} chunks "
          f"({n_workers} workers x {threads} decompression thread(s))")

    # only membership is checked in the workers, so ship a sorted position array rather than the
    # full dict: joblib memory-maps large arrays instead of pickling a copy for every chunk
    snp_pos = np.fromiter(snp_info, dtype=np.int64, count=len(snp_info))
    snp_pos.sort()
    # chunks come back in submission order, so merging below matches a serial read
    partials = Parallel(n_jobs=n_workers, backend='loky')(
        delayed(_read_bam_single)(bam, snp_pos, filtered_cells, chrom, start, end, threads)
        for bam, start, end in chunks
    )

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
            # only count bootstraps where the cell passed filtering (non-NaN score)
            n_trials = int((~np.isnan(score_array[:, i])).sum())
            if s > 0:
                n = (score_array[:, i] > 0).sum()
            elif s < 0:
                n = (score_array[:, i] < 0).sum()
            else:
                n = 0
            if n_trials == 0:
                p_value.append(1.0)
            else:
                p_value.append(binomtest(int(n), n_trials, 0.5, alternative='greater').pvalue)

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


def get_system_specs() -> dict:
    '''Collect basic machine specs (hostname, platform, CPU count, total
    memory) using only the standard library -- no new dependency. Memory is
    read from /proc/meminfo and left as None on non-Linux systems or if that
    file isn't readable.'''
    specs = {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "total_memory_gb": None,
    }
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    specs["total_memory_gb"] = round(kb / 1024 / 1024, 2)
                    break
    except (FileNotFoundError, OSError):
        pass
    return specs


def write_run_summary(summary: dict, path) -> None:
    '''Write a run summary (parameters, system specs, runtime, data/output
    counts -- see XCID.summary_stats()) to a JSON file, creating parent
    directories as needed. `default=str` is a pragmatic fallback for any
    non-JSON-native values (numpy scalars, sets, Path objects) that end up
    in the parameters/specs dicts rather than requiring every caller to
    convert them individually.'''
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
