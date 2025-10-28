import numpy as np
import pandas as pd

def make_haplotype_matrices(assignment_vector: np.ndarray, 
                         counts_ref: np.ndarray,
                         counts_alt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    '''Given an assignment vector (0/1 for each SNP), return haplotype-specific count matrices.'''

    is_hap1 = assignment_vector.astype(bool)
    is_hap0 = ~is_hap1

    hap1_counts = counts_ref.copy()
    hap0_counts = counts_alt.copy()

    hap1_counts[:, is_hap0] = counts_alt[:, is_hap0]
    hap0_counts[:, is_hap0] = counts_ref[:, is_hap0]

    return hap0_counts, hap1_counts


def make_escape_df(results: pd.DataFrame,
                    hap0_counts: np.ndarray,
                  hap1_counts: np.ndarray,
                  pos4phasing: np.ndarray,
                  score: np.ndarray,
                  min_cells=1,
                  ) -> pd.DataFrame:
    '''Make a dataframe of escape candidates based on results and haplotype count matrices.'''

    hap0_cells = score > 0 & (results['XCI_status'] != 'unknown')
    hap1_cells = score < 0 & (results['XCI_status'] != 'unknown')

    conf_cells = results['XCI_status'] != 'unknown'

    active_counts = hap0_counts.copy()
    active_counts[hap0_cells, :] = hap0_counts[hap0_cells, :]
    active_counts[hap1_cells, :] = hap1_counts[hap1_cells, :]

    inactive_counts = hap0_counts.copy()
    inactive_counts[hap0_cells, :] = hap1_counts[hap0_cells, :]
    inactive_counts[hap1_cells, :] = hap0_counts[hap1_cells, :]

    inactive_counts = inactive_counts[conf_cells, :]
    active_counts = active_counts[conf_cells, :]

    total_nonzero = (inactive_counts > 0) & (active_counts > 0)

    escape_df = pd.DataFrame({'position': pos4phasing,
                            'n_escape_cells': total_nonzero.sum(axis=0),
                             'active_counts': active_counts.sum(axis=0),
                             'inactive_counts': inactive_counts.sum(axis=0),
                             })
    escape_df = escape_df.loc[escape_df['n_escape_cells']>=min_cells]
    escape_df.sort_values('n_escape_cells', ascending=False, inplace=True)

    return escape_df