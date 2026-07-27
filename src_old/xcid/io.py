import sys
import numpy as np
import pandas as pd
import pysam
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def load_vcf(vcf_files, blacklist=frozenset(), chromosome='chrX'):
    """Loads the list of VCF files and returns a dictionary of the format: 
    {pos: {'rs':ID, 'ref':ref, 'alt':alt, 'type':'snp'}}."""
    nuc = [
        'A', 'C', 'G', 'T'
        ]
    snp_info = {}

    for path in vcf_files:
        with open(path, 'rt') as fh:   
            for line in fh:
                if line[0] == '#':
                    continue

                chrom, pos, rs, ref, alt, *_ = line.split('\t', 5)
                if chrom != chromosome:
                    continue

                pos = int(pos)
                if pos in blacklist:
                    continue

                var_type = 'snp' if (ref in nuc and alt in nuc) else 'indel'
                if var_type == 'snp':
                    snp_info[pos] = {'rs': rs, 'ref': ref, 'alt': alt,
                                     'type': 'snp'}

    n = len(snp_info)
    print(f'{n} SNPs loaded from {chromosome}')
    if not n:
        sys.exit('No SNPs found, terminating.')
    return snp_info


def read_bam(bams = [],
             snp_info = {},
             filtered_cells = None,
             chromosome = 'chrX'):
    '''Read BAM files and extract UMI and allele information for SNP positions.
    Returns a dictionary of UMIs and a mapping of cells to alleles.'''

    bad_cells = 0
    umis        = {}
    cell2allele = {}
    inf_reads   = 0
    inf_alleles = 0

    if filtered_cells is not None:
        print(f"Number of filtered cells provided: {len(filtered_cells)}")
    else:
        filtered_cells = None
        print("No filtered cells file provided, all cells will be processed.")

    for bam in bams:
        with pysam.AlignmentFile(bam, "rb") as samfile:
            # fetch reads on chromosome (or region)
            for read in samfile.fetch(chromosome):
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

                inf_reads += 1

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

                    inf_alleles += 1
                    
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

    if filtered_cells is not None:
        n_filtered_cells = len(filtered_cells)
    else:
        n_filtered_cells = 0
    print(f"Loaded {len(umis)} SNP positions with coverage.")
    print(f"Loaded {len(cell2allele)} cells with informative reads.")
    if not umis:
        sys.exit("No SNPs after reading, check the chromosome name.")
    if not cell2allele:
        sys.exit("No cells found, check the chromosome name.")
    if filtered_cells is not None:
        print(f"Filtered out {bad_cells} cells not in the filtered list.")

    return umis, cell2allele


def make_results(cells, score, score_array):
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
        "cell_id": cells,
        "score": score,
        "p_value": p_value,
        "p_adj": p_adj,
        "XCI_status": [_label(i, s) for i, s in enumerate(score)],
    })

    return results