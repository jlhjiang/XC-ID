from __future__ import annotations
from dataclasses import dataclass, field, fields
from typing import Optional, Tuple
import logging
import time
if __package__ is None or __package__ == "":
    import sys, os
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "xcid"

from .io import load_vcf, read_bam, make_results, get_system_specs
from .preprocess import get_pos_for_phasing, build_counts_matrix
from .simulated_annealing import simulated_annealing
from .bootstrap import bootstrap_variants
from .escape import make_haplotype_matrices, make_escape_df

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

@dataclass
class XCID:
    '''Class to perform X chromosome inactivation identification (XC-ID) analysis.

    Field names match the CLI flags (e.g. `chrom` <-> `--chrom`, `seed` <-> `--seed`)
    so options translate directly between the Python API and the command line.'''

    vcf: list[str] # list of VCF file paths, can provide multiple
    bam: list[str] # list of BAM file paths, can provide multiple
    cells: Optional[str] = None # path to filtered cells file (one cell barcode per line)
    chrom: str = "chrX" # chromosome to analyze, watch out for naming conventions ("chrX" vs "X")
    blacklist: set[int] = field(default_factory=set) # set of positions to blacklist (e.g. PARs)
    seed: Optional[int] = None

    min_per_snp: int = 5 # minimum UMIs per SNP to consider for phasing
    min_per_cell: int = 2 # minimum SNP counts per cell to consider for phasing
    min_maf: float = 0.1 # minimum minor allele frequency to consider SNPs
    n_boot: int = 1000
    n_jobs: int = -1
    verbose: int = 2

    # simulated annealing parameters
    temp: Optional[float] = None # initial temperature
    alpha: Optional[float] = None # cooling rate
    max_iter: Optional[int] = None # maximum iterations
    n_init: int = 10 # number of initializations
    n_init_bootstrap: int = 1 # number of initializations for bootstrap

    # internal state
    rng: np.random.Generator = field(init=False)
    _start_time: float = field(init=False, default=0.0) # wall-clock start, set in __post_init__; used by summary_stats()

    filtered_cells: Optional[np.ndarray] = field(init=False, default=None) # array of allowed cell barcodes loaded from `cells`
    snp_info: Optional[dict] = field(init=False, default=None) # SNP information dictionary
    umis: Optional[dict] = field(init=False, default=None) # UMI information dictionary
    cell2allele: Optional[dict] = field(init=False, default=None) # cell to allele mapping
    pos4phasing: Optional[list] = field(init=False, default=None) # valid positions to phase
    counts_ref: Optional[np.ndarray] = field(init=False, default=None) # reference counts (n_cells x n_snps)
    counts_alt: Optional[np.ndarray] = field(init=False, default=None) # alternate counts (n_cells x n_snps)
    cell_ids: Optional[np.ndarray] = field(init=False, default=None) # array of cell barcodes surviving filtering

    sa_score: Optional[np.ndarray] = field(init=False, default=None)
    best_cost: Optional[float] = field(init=False, default=None)
    best_assignment: Optional[np.ndarray] = field(init=False, default=None)
    score_array: Optional[np.ndarray] = field(init=False, default=None)
    results: Optional[pd.DataFrame] = field(init=False, default=None)
    hap0_counts: Optional[np.ndarray] = field(init=False, default=None)
    hap1_counts: Optional[np.ndarray] = field(init=False, default=None)
    filt_hap0_counts: Optional[np.ndarray] = field(init=False, default=None)
    filt_hap1_counts: Optional[np.ndarray] = field(init=False, default=None)
    escape_df: Optional[pd.DataFrame] = field(init=False, default=None)

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed) if self.seed is not None else np.random.default_rng()
        self._start_time = time.time()

    # ---------- compute / chainable ----------
    def load_data(self) -> "XCID":
        '''Load and preprocess data from VCF and BAM files.'''
        if self.cells is not None:
            self.filtered_cells = np.loadtxt(self.cells, dtype=str, ndmin=1)

        self.snp_info = load_vcf(self.vcf,
                                 chrom=self.chrom,
                                 blacklist=self.blacklist,
                                 )
        self.umis, self.cell2allele = read_bam(
            bams=self.bam,
            snp_info=self.snp_info,
            filtered_cells=self.filtered_cells,
            chrom=self.chrom,
            n_jobs=self.n_jobs,
        )
        self.pos4phasing = get_pos_for_phasing(umis=self.umis,
                                               min_per_snp=self.min_per_snp,
                                               min_maf=self.min_maf,
                                               blacklist=self.blacklist,
                                               )
        self.cell_ids, self.counts_ref, self.counts_alt = build_counts_matrix(
            cell2allele=self.cell2allele,
            pos4phasing=self.pos4phasing,
            min_per_cell=self.min_per_cell,
        )
        return self

    def assign_haplotypes(self) -> "XCID":
        '''Assign haplotypes using simulated annealing.'''

        self._require(self.counts_ref is not None, "Call load_data() first.")

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
            n_init=self.n_init_bootstrap,
            rng=self.rng,
        )
        return self

    # ---------- accessors / return data ----------
    def results_table(self) -> pd.DataFrame:
        '''Return a table of results with columns: cell, score, p_value, p_adj, XCI_status.'''
        self._require(self.sa_score is not None, "Call assign_haplotypes() first.")
        self._require(self.score_array is not None, "Call estimate_confidence() first.")
        results = make_results(cell_ids=self.cell_ids,
                                    score=self.sa_score,
                                    score_array=self.score_array
                                    )

        if results['XCI_status'].isin(['X0', 'X1']).sum()/len(results) < 0.2:
            log.warning("Less than 20%% of cells were confidently assigned. "
                        "This may indicate low data quality: **consider increasing n_boot (currently %d)**",
                        self.n_boot)

        self.results = results
        return self.results

    def summary_stats(self) -> dict:
        '''Build a run-summary dict: the parameters this run used, basic
        system specs, wall-clock runtime so far, and counts describing the
        data and the final results (how many SNPs/cells/reads were
        informative, and how the final X0/X1/unknown calls broke down).

        Safe to call at any point in the pipeline -- sections that depend on
        steps not yet run (e.g. output_summary before results_table()) are
        simply left as empty dicts rather than raising.'''
        parameters = {}
        for f in fields(self):
            if not f.init:
                continue
            value = getattr(self, f.name)
            if isinstance(value, set):
                value = sorted(value)
            parameters[f.name] = value

        summary = {
            "parameters": parameters,
            "runtime": {
                "elapsed_seconds": round(time.time() - self._start_time, 2),
            },
            "system": get_system_specs(),
            "data_summary": {},
            "output_summary": {},
        }

        if self.umis is not None:
            summary["data_summary"]["n_snps_with_coverage"] = len(self.umis)
        if self.pos4phasing is not None:
            summary["data_summary"]["n_informative_snps"] = len(self.pos4phasing)
        if self.cell2allele is not None:
            summary["data_summary"]["n_cells_with_informative_reads"] = len(self.cell2allele)
        if self.cell_ids is not None:
            summary["data_summary"]["n_informative_cells"] = len(self.cell_ids)
        if self.counts_ref is not None and self.counts_alt is not None:
            summary["data_summary"]["n_informative_reads"] = int(
                self.counts_ref.sum() + self.counts_alt.sum()
            )

        if self.results is not None:
            status_counts = self.results["XCI_status"].value_counts()
            summary["output_summary"] = {
                "n_X0": int(status_counts.get("X0", 0)),
                "n_X1": int(status_counts.get("X1", 0)),
                "n_unknown": int(status_counts.get("unknown", 0)),
                "n_total_cells": int(len(self.results)),
            }

        return summary

# ---------- experimental haplotype matrices and escape ----------

    def haplotype_matrices(self) -> Tuple[np.ndarray, np.ndarray]:
        '''Return haplotype count matrices (n_cells x n_snps) for reference and alternate alleles.'''
        self._require(self.best_assignment is not None, "Call assign_haplotypes() first.")
        self.hap0_counts, self.hap1_counts = make_haplotype_matrices(
            self.best_assignment,
            self.counts_ref,
            self.counts_alt,
            )
        return self.hap0_counts, self.hap1_counts

    def filtered_haplotype_matrices(self) -> Tuple[np.ndarray, np.ndarray]:
        '''Return haplotype count matrices (n_filtered_cells x n_snps) for reference and alternate alleles.'''
        self._require(self.best_assignment is not None, "Call assign_haplotypes() first.")
        self._require(self.results is not None, "Call results_table() first.")
        filt_hap0_counts, filt_hap1_counts = make_haplotype_matrices(
            self.best_assignment,
            self.counts_ref,
            self.counts_alt,
            )

        conf_idx = self.results['XCI_status'] != 'unknown'
        conf_cells = self.results.loc[conf_idx, 'cell_id'].values
        filt_hap0_counts = filt_hap0_counts[conf_idx, :]
        filt_hap1_counts = filt_hap1_counts[conf_idx, :]

        filt_hap0_counts = pd.DataFrame(filt_hap0_counts, index=conf_cells, columns=self.pos4phasing)
        filt_hap1_counts = pd.DataFrame(filt_hap1_counts, index=conf_cells, columns=self.pos4phasing)

        self.filt_hap0_counts = filt_hap0_counts
        self.filt_hap1_counts = filt_hap1_counts

        return self.filt_hap0_counts, self.filt_hap1_counts

    def escape_table(self, min_cells: int = 1) -> pd.DataFrame:
        '''Return a table of escape candidates with columns: position, n_escape_cells, active_counts, inactive_counts.'''
        if self.hap0_counts is None or self.hap1_counts is None:
            self.haplotype_matrices()
        self._require(self.results is not None, "Call results_table() first.")
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
