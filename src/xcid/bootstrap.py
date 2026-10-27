import numpy as np
from joblib import Parallel, delayed
from .simulated_annealing import simulated_annealing
from .scoring import get_score


def bootstrap_once(k: int,
                   pos4phasing: list,
                   sample_idx: np.ndarray,
                   counts_ref: np.ndarray,
                   counts_alt: np.ndarray,
                   sa_kwargs=None,
                   min_per_cell=2,
                   rng=None,
                   ):
    '''Perform a single bootstrap iteration for haplotype confidence estimation.'''

    # extract counts for the sample, applying filter
    count_ref_sample = counts_ref[:, sample_idx]
    count_alt_sample = counts_alt[:, sample_idx]
    total_counts = (count_ref_sample + count_alt_sample)
    filter_mask = total_counts.sum(axis=1) >= min_per_cell
    count_ref_sample = count_ref_sample[filter_mask, :]
    count_alt_sample = count_alt_sample[filter_mask, :]

    # Run simulated annealing on the sampled counts
    best_assignment_sample, _, _ = simulated_annealing(
        counts_ref=count_ref_sample, counts_alt=count_alt_sample,
        pos4phasing=pos4phasing, print_interval=False,
        rng=rng, 
        **sa_kwargs)  # use the filtered counts

    # but get the values for all cells
    score_sample = get_score(count_ref_sample, count_alt_sample, best_assignment_sample)

    bootstrap_tray = {'k': k,  # store the iteration number
                        'cells_mask': filter_mask,
                        'best_a': best_assignment_sample,
                        'score': score_sample,
                        }
    return bootstrap_tray  # return the result for this bootstrap iteration


def bootstrap_variants(pos4phasing: list,
                        counts_ref: np.ndarray, counts_alt: np.ndarray,
                        obs_assignment: np.ndarray,
                        n_boot=100,
                        min_per_cell=2,
                        n_jobs=-1,  # total - 1 number of jobs
                        verbose=2,
                        min_score=0,
                        rng=None,
                        **sa_kwargs):
    """
    Parallel version of bootstrap_variants using joblib.ProcessPool (loky).
    counts_ref / counts_alt may be large NumPy arrays; they are mem-mapped
    read-only so workers don't copy gigabytes of data.
    """
    print(f'Performing {n_boot} bootstraps to determine cell confidence...')

    obs_score = get_score(counts_ref, counts_alt, np.array(obs_assignment))
    nonzero_score = (obs_score != 0) & (np.abs(obs_score) >= min_score)

    # boolean fancy-indexing already returns a fresh array, no extra .copy() needed
    counts_ref_nz = counts_ref[nonzero_score, :]
    counts_alt_nz = counts_alt[nonzero_score, :]

    n_pos = len(pos4phasing)
    sample_idx = rng.choice(n_pos, (n_boot, n_pos), replace=True) # sample index number of snps with replacement
    child_seeds = rng.integers(2**32, size=n_boot)

    results = Parallel(
        n_jobs=n_jobs,
        backend='loky',
        mmap_mode='r',
        verbose=verbose, # progress bar
    )(
        delayed(bootstrap_once)(
            k=k,
            pos4phasing=pos4phasing,
            sample_idx=sample_idx[k, :],
            counts_ref=counts_ref_nz,
            counts_alt=counts_alt_nz,
            min_per_cell=min_per_cell,
            sa_kwargs=sa_kwargs,
            rng=np.random.default_rng(int(child_seeds[k])),
        )
        for k in range(n_boot)
    )

    w = (counts_ref_nz + counts_alt_nz).sum(axis=0)
    obs_assignment = np.asarray(obs_assignment)

    score_array = np.full((n_boot, counts_ref.shape[0]), np.nan, dtype=float)
    for k, r in enumerate(results):
        # sample_idx[k, :] are the same integer SNP indices used to build this
        # bootstrap sample, so we can reuse them directly instead of re-deriving
        # them with an O(n_pos) pos4phasing.index() lookup per sampled position
        # (which made this loop O(n_boot * n_pos^2) for no reason).
        pos_indices = sample_idx[k, :]
        assignment_from_obs = obs_assignment[pos_indices]
        weight = w[pos_indices]

        best_a = np.asarray(r['best_a'])
        cell_mask = np.zeros(counts_ref.shape[0], dtype=bool)
        cell_mask[nonzero_score] = r['cells_mask']

        if ((assignment_from_obs != best_a) * weight).sum() > ((assignment_from_obs == best_a) * weight).sum():
            score_array[r['k'], cell_mask] = -r['score']  # flip
        else:
            score_array[r['k'], cell_mask] = r['score']

    return score_array
