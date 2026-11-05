#!/usr/bin/env python3
"""Generate per-reference coverage depth plots from samtools output files."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot depth across reference positions from samtools depth TSV files."
    )
    parser.add_argument(
        "--depth-file",
        required=True,
        help="samtools depth TSV file with columns: reference, position, depth.",
    )
    parser.add_argument(
        "--coverage-file",
        help="Optional samtools coverage file used to extract numreads per reference.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output PNG file path. Defaults to the input file stem with .png.",
    )
    return parser.parse_args()


def read_depth_file(path):
    """Read a samtools depth TSV file with reference, position, and depth columns."""
    refs = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for ref, pos, depth in reader:
            refs.setdefault(ref, {"positions": [], "depths": []})
            refs[ref]["positions"].append(int(pos))
            refs[ref]["depths"].append(int(depth))
    return refs


def read_coverage_file(path):
    """Read samtools coverage file.

    Expected columns: rname, startpos, endpos, numreads, covbases, coverage,
    meandepth, meanbaseq, meanmapq.
    """
    coverage_by_ref = {}
    if not path:
        return coverage_by_ref

    # Reuse the samtools coverage summary so each subplot title can include
    # the total number of mapped reads for the matching reference sequence.
    with open(path, newline="") as handle:
        header = handle.readline().strip().lstrip("#").split("\t")
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
        for row in reader:
            if row.get("rname"):
                coverage_by_ref[row["rname"]] = row
    return coverage_by_ref


def plot_depth_file(depth_file, output, coverage_file=None):
    """Plot depth by reference"""
    refs = read_depth_file(depth_file)
    coverage_by_ref = read_coverage_file(coverage_file)
    ref_summaries = []

    for ref, values in refs.items():
        positions = values["positions"]
        depths = values["depths"]
        coverage_row = coverage_by_ref.get(ref, {})
        covered_bases = sum(1 for x in depths if x > 0)
        coverage_pct = 100 * covered_bases / len(depths) if depths else 0
        max_depth = max(depths) if depths else 0
        ref_summaries.append(
            {
                "ref": ref,
                "plot_positions": positions,
                "plot_depths": depths,
                "coverage_pct": coverage_pct,
                "max_depth": max_depth,
                "mapped_reads": coverage_row.get("numreads"),
            }
        )

    # Show the most supported references first so the more informative panels
    # are easier to find when one plot contains many accessions.
    ref_summaries.sort(key=lambda item: (item["coverage_pct"], item["max_depth"]), reverse=True)

    n_refs = len(ref_summaries)
    ncols = 1 if n_refs == 1 else 2
    nrows = math.ceil(n_refs / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(14, max(3.2 * nrows, 4.5)),
        squeeze=False,
        sharey=False,
    )

    for ax, summary in zip(axes.flat, ref_summaries):
        ax.plot(
            summary["plot_positions"],
            summary["plot_depths"],
            linewidth=1.6,
            color="#3f6fb5",
            alpha=0.95,
            solid_capstyle="round",
            solid_joinstyle="round",
        )
        ax.set_title(
            f"{summary['ref']} | mapped_reads={summary['mapped_reads'] or 'NA'} | coverage={summary['coverage_pct']:.2f}% | max_depth={summary['max_depth']}",
            fontsize=9,
        )
        ax.set_xlabel("Position")
        ax.set_ylabel("Depth")
        ax.grid(color="#d9d9d9", linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)

    for ax in axes.flat[n_refs:]:
        # Hide any unused panels in the final grid.
        ax.axis("off")

    plt.tight_layout()
    outfile = Path(output)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outfile


def main():
    args = parse_args()
    output = args.output
    outfile = plot_depth_file(args.depth_file, output, args.coverage_file)
    print(outfile)


if __name__ == "__main__":
    main()
