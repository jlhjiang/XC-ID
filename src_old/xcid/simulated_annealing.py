import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from .scoring import get_score


def cell_cost_function_multi(assignment_matrix, counts_ref, counts_alt):
    '''cell wise cost function (maximizing haplotype coverage) for vectorized multiple initializations'''
    # no need for mask, just dot products of counts and assignment_matrix
    # assignment_matrix is n_pos x n_init

    hap0_cov = np.dot(counts_ref, (1 - assignment_matrix).T) + np.dot(counts_alt, assignment_matrix.T) # n_cells x n_init
    # sum of hap on ref counts and hap on alt counts
    hap1_cov = np.dot(counts_alt, (1 - assignment_matrix).T) + np.dot(counts_ref, assignment_matrix.T)

    cell_cost = hap0_cov - hap1_cov

    return cell_cost


def simulated_annealing(counts_ref, counts_alt, 
                        max_iter=None,
                        pos4phasing=None,
                        temp=None, alpha=0.999, 
                        print_interval=True, 
                        n_init=1,
                        rng=None,
                        ):
    '''Simulated annealing to assign haplotypes to SNPs based on allele counts across cells.'''

    n_pos = len(pos4phasing)

    if max_iter is None:
        max_iter = np.maximum(15000, n_pos * 20).astype(int) 
        # ensure max_iter is at least 15000 or 1000 times the number of positions
    if temp is None:
        temp = n_pos * 1000
    
    if n_pos < 2500:
        alpha = 0.995
    else:
        alpha = 0.999

    t_schedule = temp * alpha ** np.arange(max_iter)

    if print_interval:
        print(f"Running simulated annealing with {n_init} initializations, {max_iter} max iterations.")
        print(f"Initial temperature: {temp}, cooling rate: {alpha}.")

    print_schedule = max_iter // 4

    pos4phasing = np.asarray(pos4phasing)

    assignment_matrix = rng.choice([0, 1], size=(n_init, n_pos))

    # Compute initial cost
    init_cell_cost = cell_cost_function_multi(assignment_matrix, counts_ref, counts_alt).astype(np.float32)
    init_cost = -np.abs(init_cell_cost).sum(axis=0)

    best_assignment = assignment_matrix.copy()
    best_cost = init_cost.copy()
    best_cell_cost = init_cell_cost.copy()

    # precompute the difference in counts
    diff = (counts_ref - counts_alt).astype(np.float32)
    
    for iteration in range(max_iter):

        temp = t_schedule[iteration]
        idx = rng.integers(0, n_pos, size=n_init)
        flip = 1 - 2*best_assignment[np.arange(n_init), idx] # directly + or - the cost difference

        delta_vec = 2 * diff[:, idx] * flip

        cur_abs = np.abs(best_cell_cost)
        new_abs = np.abs(best_cell_cost + delta_vec)
        delta_abs = new_abs - cur_abs
        delta_E = - delta_abs.sum(axis=0)
        new_cost = best_cost + delta_E

        # accept this best neighbour using the simulated annealing acceptance criterion
        eps = 1e-10
        worse  = delta_E >  eps
        better = delta_E < -eps
        if temp > 1e-16:
            cond_accept = np.log(rng.random(n_init)) < (-delta_E / temp)
        else:
            cond_accept = np.zeros(n_init, dtype=bool)
        
        accept = better | (worse & cond_accept)

        if accept.any():
            new_cost = best_cost + delta_E
            accepting = np.where(accept)[0]
            best_cost[accepting] = new_cost[accepting]
            best_cell_cost[:, accepting] += delta_vec[:, accepting]
            best_assignment[accepting, idx[accepting]] ^= 1

        if print_interval:
            if iteration % print_schedule == 0:
                print(f'Iteration {iteration}: cost={best_cost}')
                print(f'Current temperature: {t_schedule[iteration].round(2)}')

    best_run = np.argmin(best_cost)
    final_best_assignment = best_assignment[best_run, :]
    final_best_cost = best_cost[best_run]

    final_score = get_score(counts_ref, counts_alt, final_best_assignment)

    if print_interval:
        print(f'Best cost found after {iteration} iterations: {best_cost[best_run]}. Initial cost: {init_cost[best_run]}.')

    return final_best_assignment, final_best_cost, final_score