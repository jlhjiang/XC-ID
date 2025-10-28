import sys
import numpy as np

def get_pos_for_phasing(umis=None,
                        min_maf=0.1, #minimum minor allele frequency to consider for phasing 
                        min_per_snp=5, #minimum UMIs per position to consider for phasing
                        blacklist=frozenset() #optional addtional blacklisting
                        ):
    """
    Filter snps based on UMIs, MAF, and blacklist.
    Returns a list of positions to phase and their corresponding UMI counts.
    """
    if not umis:
        sys.exit("No UMIs provided for filtering.")

    print(f"Filtering SNPs based on minimum coverage={min_per_snp}, MAF={min_maf}")
    
    pos4phasing = []
    snp_counts = []

    low_umis   = 0
    monoallelic = 0
    low_maf = 0

    for posx in sorted(umis.keys()):
        count_ref = len(umis[posx]['ref'])
        count_alt = len(umis[posx]['alt'])

        if posx in blacklist:
            continue
        if count_ref < 1 or count_alt < 1:
            monoallelic +=1
            continue
        elif count_ref+count_ref < min_per_snp: #(count_ref < min_per_snp) or (count_alt < min_per_snp):
            low_umis +=1
            continue
        elif (np.minimum(count_ref, count_alt)/(count_ref+count_alt)) < min_maf:
            low_maf += 1
            continue
        else:
            pos4phasing.append(posx)
            snp_counts.append(count_ref+count_alt)

    if not pos4phasing:
        sys.exit("No SNPs to phase after filtering.")
    print(f"Final SNPs to phase: {len(pos4phasing)}")

    return pos4phasing


def build_counts_matrix(cell2allele=None,
                        pos4phasing=None,
                        min_per_cell=2):
    '''Build REF and ALT counts matrices (n_cells x n_snps) for the given cell2allele dictionary and positions for phasing.'''

    cells = list(cell2allele.keys())
    n_cells = len(cells)
    n_positions = len(pos4phasing)
    counts_ref = np.zeros((n_cells, n_positions), dtype=int)
    counts_alt = np.zeros((n_cells, n_positions), dtype=int)

    pos_index = {pos: j for j, pos in enumerate(pos4phasing)}

    for i, cell in enumerate(cells):
        for pos, umi_dict in cell2allele[cell].items():
            # Only include positions in pos4phasing
            if pos in pos_index:
                j = pos_index[pos]
                counts_ref[i, j] = np.array([1 for allele in umi_dict.values() if allele == 'ref']).sum()
                counts_alt[i, j] = np.array([1 for allele in umi_dict.values() if allele == 'alt']).sum()
    filter_mask = ((counts_ref+counts_alt).sum(axis=1) >= min_per_cell) # filter for minimum counts in a cell
    print(f"Remaining informative cells: {filter_mask.sum()} out of {n_cells}.")

    cells = np.array(cells)[filter_mask]
    counts_ref = counts_ref[filter_mask,:]
    counts_alt = counts_alt[filter_mask,:]

    return cells, np.asarray(counts_ref, dtype=np.int32), np.asarray(counts_alt, dtype=np.int32)