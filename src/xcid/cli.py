# xcid/cli.py
from __future__ import annotations
import argparse
from pathlib import Path

from .xcid import XCID
from .io import load_blacklist, write_results_table, write_escape_table, write_haplotype_matrices, write_run_summary

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="xcid",
        description="XC-ID: X Chromosome inactivation IDentifier"
    )
    io = p.add_argument_group("I/O")
    io.add_argument("--vcf", nargs="+", required=True, help="One or more VCF files (chrX variants).")
    io.add_argument("--bam", nargs="+", required=True, help="One or more BAM files (STAR+WASP aligned).")
    io.add_argument("--cells", default=None, help="Optional file of allowed cell barcodes (one per line).")
    io.add_argument("--chrom", default="chrX", help='Chromosome name in BAM/VCF (default: "chrX").')
    io.add_argument("--blacklist", default=None, help="Optional file of SNP positions to exclude (one 1-based pos per line).")

    out = p.add_argument_group("Outputs")
    out.add_argument("--out-results", required=True, help="Path to write results table TSV.")

    filt = p.add_argument_group("Filtering")
    filt.add_argument("--min-per-snp", type=int, default=5, help="Min UMIs per SNP to keep (default: 5).")
    filt.add_argument("--min-per-cell", type=int, default=2, help="Min SNP counts per cell to keep (default: 2).")
    filt.add_argument("--min-maf", type=float, default=0.05, help="Min minor allele frequency to keep a SNP (default: 0.05).")

    sa = p.add_argument_group(
        "Simulated annealing (do not change these unless you know what you are doing)"
    )
    sa.add_argument("--temp", type=float, default=None, help="Initial temperature (default: n_pos * 1000).")
    sa.add_argument("--alpha", type=float, default=None, help="Cooling rate (default auto: 0.995/0.999).")
    sa.add_argument("--max-iter", type=int, default=None, help="Max iterations (default: max(15000, 20*n_pos)).")
    sa.add_argument("--n-init", type=int, default=10, help="Number of random initializations (default: 1).")

    boot = p.add_argument_group("Bootstrap")
    boot.add_argument("--n-boot", type=int, default=100, help="Number of bootstraps (default: 100).")
    boot.add_argument("--n-jobs", type=int, default=-1, help="Parallel jobs for BAM reading and bootstrap (default: -1, all cores).")
    boot.add_argument("--verbose", type=int, default=2, help="joblib verbosity (default: 2).")

    misc = p.add_argument_group("Reproducibility")
    misc.add_argument("--seed", type=int, default=None, help="Random seed for RNG streams.")

    exp = p.add_argument_group(
        "Experimental (for development/testing purposes only)"
    )
    exp.add_argument("--out-hap-counts", default=None, help="Optional output directory for haplotype count matrices.")
    exp.add_argument("--out-escape", default=None, help="Optional TSV of escape candidates.")

    return p

def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    blacklist = load_blacklist(args.blacklist)

    x = XCID(
        vcf=args.vcf,
        bam=args.bam,
        cells=args.cells,
        chrom=args.chrom,
        blacklist=blacklist,
        seed=args.seed,
        min_per_snp=args.min_per_snp,
        min_per_cell=args.min_per_cell,
        min_maf=args.min_maf,
        n_boot=args.n_boot,
        n_jobs=args.n_jobs,
        verbose=args.verbose,
        temp=args.temp,
        alpha=args.alpha,
        max_iter=args.max_iter,
        n_init=args.n_init,
    )

    # Run full pipeline
    results = x.run_all()

    # Save results
    write_results_table(results, args.out_results)

    # Run summary (parameters, system specs, runtime, data/output counts) --
    # saved alongside the results table, not behind a flag, since it's cheap
    # to produce and useful for every run.
    summary_path = Path(args.out_results).parent / "run_summary.json"
    write_run_summary(x.summary_stats(), summary_path)

    # Optional: save escape table
    if args.out_escape:
        write_escape_table(x.escape_table(min_cells=1), args.out_escape)

    # Optional: save haplotype count matrices
    if args.out_hap_counts:
        hap0_counts, hap1_counts = x.filtered_haplotype_matrices()
        write_haplotype_matrices(hap0_counts, hap1_counts, args.out_hap_counts)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
