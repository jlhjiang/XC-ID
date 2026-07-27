import numpy as np

def get_score(counts_ref, counts_alt, assignment_vector):
    '''get the score for each cell given the assignment vector'''

    hap0_cov = np.dot(counts_ref, (1 - assignment_vector)) + np.dot(counts_alt, assignment_vector)
    hap1_cov = np.dot(counts_alt, (1 - assignment_vector)) + np.dot(counts_ref, assignment_vector)
   
    score = np.log2(hap0_cov + 1) - np.log2(hap1_cov + 1)

    inf_mask = np.isinf(score)
    score[inf_mask] = np.nan
    return score