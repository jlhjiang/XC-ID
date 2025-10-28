from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Tuple
import logging
if __package__ is None or __package__ == "":
    import sys, os
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "xcid"

from .io import load_vcf, read_bam, make_results
from .preprocess import get_pos_for_phasing, build_counts_matrix
from .simulated_annealing import simulated_annealing
from .bootstrap import bootstrap_variants
from .escape import make_haplotype_matrices, make_escape_df

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

@dataclass
class XCID:
    '''Class to perform X chromosome inactivation identification (XC-ID) analysis.'''

    vcf_files: list[str] # list of VCF file paths, can provide multiple
    bam_files: list[str] # list of BAM file paths, can provide multiple
    cells_file: Optional[str] = None # path to filtered cells file (one cell barcode per line)
    chromosome: str = "X" # chromosome to analyze, watch out for naming conventions ("e.g. chrX" vs "X")
    blacklist: set[int] = field(default_factory=set) # set of positions to blacklist (e.g. PARs)
    rand_seed: Optional[int] = None

    min_per_snp: int = 5 # minimum UMIs per SNP to consider for phasing
    min_per_cell: int = 2 # minimum SNP counts per cell to consider for phasing
    n_boot: int = 100
    n_jobs: int = -1
    verbose: int = 2

    # simulated annealing parameters
    temp: Optional[int] = None # initial temperature
    alpha: Optional[float] = None # cooling rate
    max_iter: Optional[int] = None # maximum iterations
    n_init: int = 1 # number of initializations

    # internal state
    rng: np.random.Generator = field(init=False)

    filtered_cells: Optional[np.ndarray] = field(init=False, default=None) # array of filtered cell barcodes
    snp_info: Optional[dict] = field(init=False, default=None) # SNP information dictionary
    umis: Optional[dict] = field(init=False, default=None) # UMI information dictionary
    cell2allele: Optional[dict] = field(init=False, default=None) # cell to allele mapping
    pos4phasing: Optional[list] = field(init=False, default=None) # valid positions to phase
    counts_ref: Optional[np.ndarray] = field(init=False, default=None) # reference counts (n_cells x n_snps)
    counts_alt: Optional[np.ndarray] = field(init=False, default=None) # alternate counts (n_cells x n_snps)
    cells: Optional[np.ndarray] = field(init=False, default=None) # array of cell barcodes

    sa_score: Optional[np.ndarray] = field(init=False, default=None)
    best_cost: Optional[float] = field(init=False, default=None)
    best_assignment: Optional[np.ndarray] = field(init=False, default=None)
    score_array: Optional[np.ndarray] = field(init=False, default=None)
    results: Optional[pd.DataFrame] = field(init=False, default=None)
    hap0_counts: Optional[np.ndarray] = field(init=False, default=None)
    hap1_counts: Optional[np.ndarray] = field(init=False, default=None)
    escape_df: Optional[pd.DataFrame] = field(init=False, default=None)

    def __post_init__(self):
        self.rng = np.random.default_rng(self.rand_seed) if self.rand_seed is not None else np.random.default_rng()

    # ---------- compute / chainable ----------
    def load_data(self) -> "XCID":
        '''Load and preprocess data from VCF and BAM files.'''
        if self.cells_file is not None:
            self.filtered_cells = np.loadtxt(self.cells_file, dtype=str)

        self.snp_info = load_vcf(self.vcf_files, 
                                 chromosome=self.chromosome, 
                                 blacklist=self.blacklist,
                                 )
        self.umis, self.cell2allele = read_bam(
            bams=self.bam_files,
            snp_info=self.snp_info,
            filtered_cells=self.filtered_cells,
            chromosome=self.chromosome,
        )
        self.pos4phasing = get_pos_for_phasing(umis=self.umis, 
                                               min_per_snp=self.min_per_snp, 
                                               blacklist=self.blacklist,
                                               )
        self.cells, self.counts_ref, self.counts_alt = build_counts_matrix(
            cell2allele=self.cell2allele,
            pos4phasing=self.pos4phasing,
            min_per_cell=self.min_per_cell,
        )
        return self

    def assign_haplotypes(self) -> "XCID":
        '''Assign haplotypes using simulated annealing.'''

        self._require(self.counts_ref is not None, "Call load_data() first.")
        # Derive a child seed for this stage to keep streams independent

        self.best_assignment, self.best_cost, self.sa_score = simulated_annealing(
            self.counts_ref,
            self.counts_alt,
            pos4phasing=self.pos4phasing,
            temp=self.temp,
            alpha=self.alpha,
            max_iter=self.max_iter,
            n_init=self.n_init,
            rng=self.rng,
        )
        return self

    def estimate_confidence(self) -> "XCID":
        '''Estimate confidence of haplotype assignments via bootstrapping.'''

        self._require(self.best_assignment is not None, "Call assign_haplotypes() first.")

        self.score_array = bootstrap_variants(
            pos4phasing=self.pos4phasing,
            obs_assignment=self.best_assignment,
            counts_ref=self.counts_ref,
            counts_alt=self.counts_alt,
            n_boot=self.n_boot,
            min_per_cell=self.min_per_cell,
            n_jobs=self.n_jobs,
            verbose=self.verbose,
            temp=self.temp,
            alpha=self.alpha,
            max_iter=self.max_iter,
            n_init=self.n_init,
            rng=self.rng,
        )
        return self

    # ---------- accessors / return data ----------
    def results_table(self) -> pd.DataFrame:
        '''Return a table of results with columns: cell, score, p_value, p_adj, XCI_status.'''
        self._require(self.sa_score is not None, "Call assign_haplotypes() first.")
        self.results = make_results(cells=self.cells,
                                    score=self.sa_score,
                                    score_array=self.score_array
                                    )
        return self.results

    def haplotype_matrices(self) -> Tuple[np.ndarray, np.ndarray]:
        '''Return haplotype count matrices (n_cells x n_snps) for reference and alternate alleles.'''
        self._require(self.best_assignment is not None, "Call assign_haplotypes() first.")
        self.hap0_counts, self.hap1_counts = make_haplotype_matrices(
            self.best_assignment,
            self.counts_ref,
            self.counts_alt,
            )
        return self.hap0_counts, self.hap1_counts

    def escape_table(self, min_cells: int = 1) -> pd.DataFrame:
        '''Return a table of escape candidates with columns: position, n_escape_cells, active_counts, inactive_counts.'''
        if self.hap0_counts is None or self.hap1_counts is None:
            self.haplotype_matrices()
        self._require(self.results is not None, "Call results_df() first.")
        self.escape_df = make_escape_df(
            results=self.results,
            hap0_counts=self.hap0_counts,
            hap1_counts=self.hap1_counts,
            pos4phasing=self.pos4phasing,
            score=self.sa_score,
            min_cells=min_cells,
            )
        return self.escape_df

    # ---------- convenience ----------
    def run_all(self) -> pd.DataFrame:
        '''Run the full XCID pipeline and return the results table.'''
        log.info("Running XCID: load → SA → bootstrap → results")

        print("1. Loading data and constructing count matrices...")
        self.load_data()
        print('==================================')

        print("2. Assigning haplotypes via simulated annealing...")
        self.assign_haplotypes()
        print('==================================')

        print("3. Estimating confidence via bootstrapping...")
        self.estimate_confidence()
        print('==================================')

        print("4. Compiling results table...")
        out = self.results_table()
        print("Done!")
        return out

    # ---------- guards ----------
    @staticmethod
    def _require(cond: bool, msg: str):
        if not cond:
            raise RuntimeError(msg)
