#!/usr/bin/env python

import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Filter and summarize the BLAST results."
    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("--header", required=True, help="Path to a BLASTn header file.")
    parser.add_argument("--input", required=True, help="Path to the raw BLASTn result file.")
    parser.add_argument("--filtered-output", required=True, help="Output file for filtered BLASTn results.")
    parser.add_argument("--summary-output", required=True, help="Output file for summarized BLASTn results.")
    parser.add_argument("--min-qlen", type=int, default=50, help="Minimum query sequence length (default: 50)")
    parser.add_argument("--min-pident", type=int, default=50, help="Minimum percentage of identical matches (default: 50)")
    parser.add_argument("--min-length", type=int, default=50, help="Minimum alignment length (default: 50)")
    parser.add_argument("--max-evalue", type=float, default=0.001, help="Maximum expected value (default: 0.001)")

    return parser.parse_args(args)


def filter_summary_blast(blast_header, input_file, filtered_output, summary_output, min_qlen, min_pident, min_length, max_evalue):
    # Load the header of BLASTn results
    with open(blast_header) as f:
        header_line = f.readline().strip()
        column_names = header_line.split("\t")

    # Load BLASTn results
    raw_results = pd.read_csv(input_file, sep="\t", header=None, names=column_names, keep_default_na=False)

    # Remove entries with missing or invalid staxids
    raw_results = raw_results[
        raw_results["staxids"].notna() &
        (raw_results["staxids"].astype(str).str.strip() != "") &
        (~raw_results["staxids"].astype(str).isin(["NA", "N/A"]))
    ]

    # Apply filtering
    filtered = raw_results[
        (raw_results["qlen"] >= min_qlen) &
        (raw_results["pident"] >= min_pident) &
        (raw_results["length"] >= min_length) &
        (raw_results["evalue"] <= max_evalue)
    ]

    # Remove entries with missing scientific name
    filtered = filtered[
        filtered["sscinames"].notna() &
        (filtered["sscinames"].astype(str).str.strip() != "") &
        (~filtered["sscinames"].astype(str).isin(["NA", "N/A"]))
    ]

    # Check if any blast hits pass the filtering threshold
    if filtered.empty:
        return

    # Save filtered DataFrame
    filtered.to_csv(filtered_output, index=False, sep="\t")

    # Summarize filtered results
    summary = (
        filtered.groupby(["qseqid", "staxids", "sscinames"])
        .agg(
            count=("staxids", "count"),
            min_pident=("pident", "min"),
            max_pident=("pident", "max"),
            median_pident=("pident", "median"),
            min_length=("length", "min"),
            max_length=("length", "max"),
            median_length=("length", "median"),
            min_bitscore=("bitscore", "min"),
            max_bitscore=("bitscore", "max"),
            median_bitscore=("bitscore", "median")
        )
        .reset_index()
    )

    # Save summary DataFrame
    summary.to_csv(summary_output, index=False, sep="\t")


def main(args=None):
    args = parse_args(args)
    filter_summary_blast(
        blast_header=args.header,
        input_file=args.input,
        filtered_output=args.filtered_output,
        summary_output=args.summary_output,
        min_qlen=args.min_qlen,
        min_pident=args.min_pident,
        min_length=args.min_length,
        max_evalue=args.max_evalue
    )


if __name__ == "__main__":
    main()
