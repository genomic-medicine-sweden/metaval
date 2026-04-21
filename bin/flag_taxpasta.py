#!/usr/bin/env python3
"""
Compare one sample taxpasta profile against one optional NTC taxpasta profile.

Writes one TSV with columns:
  taxonomy_id, name, rank, lineage, <sample>
  or with NTC: ..., <sample>_vs_<ntc>, <ntc>
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

META_COLS = ["taxonomy_id", "name", "rank", "lineage"]
REQ_COLS = {"taxonomy_id", "name"}

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Flag taxpasta tables against NTC tables.")
    p.add_argument("--sample_taxpasta", required=True, help="Taxpasta TSV containing sample columns.")
    p.add_argument("--sample_name", required=True, help="Sample column name to select.")
    p.add_argument("--ntc_taxpasta", help="Optional NTC taxpasta TSV.")
    p.add_argument("--ntc_name", help="Optional NTC column name to select.")
    p.add_argument(
        "--prefix",
        help="Output filename stem. Defaults to '<sample_taxpasta_stem>_<sample_name>'.",
    )
    p.add_argument("-o", "--outdir", help="Output directory. Defaults to sample taxpasta directory.")
    return p.parse_args(argv)

#def safe_filename_part(name: str) -> str:
#    return re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip("-_.")

def load_taxpasta_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    if missing := REQ_COLS - set(df.columns):
        raise ValueError(f"'{path}' missing required columns: {', '.join(sorted(missing))}")

    for col in ("rank", "lineage"):
        if col not in df.columns:
            df[col] = ""

    profile_cols = [c for c in df.columns if c not in META_COLS]
    if not profile_cols:
        raise ValueError(f"'{path}' has no profile columns.")

    df = df[META_COLS + profile_cols].copy()
    df[["rank", "lineage"]] = df[["rank", "lineage"]].fillna("")
    df[profile_cols] = df[profile_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df


def flag(sample_reads_count: int, ntc_reads_count: int) -> str:
    if sample_reads_count > 0 and ntc_reads_count == 0:
        return "in_sample"
    if sample_reads_count == 0 and ntc_reads_count > 0:
        return "in_NTC"
    if sample_reads_count > ntc_reads_count:
        return "more"
    if sample_reads_count < ntc_reads_count:
        return "less"
    else:
        return "equal"

def find_profile_column(df: pd.DataFrame, name: str, path: str) -> str:
    profile_cols = [c for c in df.columns if c not in META_COLS]
    if name in profile_cols:
        return name

    pattern = re.compile(rf"^{re.escape(name)}(?:[._-].*)?$")
    matches = sorted([c for c in profile_cols if pattern.match(c)], key=len)

    if not matches:
        raise ValueError(
            f"No column matching '{name}' in '{path}'. Available: {', '.join(profile_cols)}"
        )

    # Return the unique shortest match, otherwise raise
    shortest = [c for c in matches if len(c) == len(matches[0])]
    if len(shortest) == 1:
        return shortest[0]

    raise ValueError(f"Ambiguous columns for '{name}' in '{path}': {', '.join(matches)}")


def build_output(
    sample_df: pd.DataFrame,
    sample_col: str,
    sample_label: str,
    ntc_df: pd.DataFrame | None = None,
    ntc_col: str | None = None,
    ntc_label: str | None = None,
) -> pd.DataFrame:
    base = sample_df[META_COLS + [sample_col]].copy()

    if ntc_df is None or ntc_col is None or ntc_label is None:
        out = base.loc[base[sample_col] > 0].rename(columns={sample_col: sample_label})
        return out[META_COLS + [sample_label]]

    ntc_slice = ntc_df[META_COLS + [ntc_col]].rename(columns={ntc_col: ntc_label})
    merged = base.merge(ntc_slice, on=META_COLS, how="outer")
    merged[sample_col] = merged[sample_col].fillna(0).astype(int)
    merged[ntc_label] = merged[ntc_label].fillna(0).astype(int)
    merged[["rank", "lineage"]] = merged[["rank", "lineage"]].fillna("")

    merged = merged.loc[(merged[sample_col] > 0) | (merged[ntc_label] > 0)].copy()
    cmp_label = f"{sample_label}_vs_{ntc_label}"
    merged[cmp_label] = [flag(sample_reads_count, ntc_reads_count) for sample_reads_count, ntc_reads_count in zip(merged[sample_col], merged[ntc_label])]
    merged = merged.rename(columns={sample_col: sample_label})

    return merged[META_COLS + [sample_label, cmp_label, ntc_label]]


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    sample_df = load_taxpasta_file(args.sample_taxpasta)
    sample_col = find_profile_column(sample_df, args.sample_name, args.sample_taxpasta)

    ntc_df, ntc_col = None, None
    if args.ntc_taxpasta and args.ntc_name:
        ntc_df = load_taxpasta_file(args.ntc_taxpasta)
        ntc_col = find_profile_column(ntc_df, args.ntc_name, args.ntc_taxpasta)

    outdir = Path(args.outdir) if args.outdir else Path(args.sample_taxpasta).parent
    outdir.mkdir(parents=True, exist_ok=True)

    output_df = build_output(sample_df, sample_col, args.sample_name, ntc_df, ntc_col, args.ntc_name)
    filename_stem = args.prefix or f"{Path(args.sample_taxpasta).stem}_{safe_filename_part(args.sample_name)}"
    out_path = outdir / f"{filename_stem}.tsv"
    output_df.to_csv(out_path, sep="\t", index=False)

    ntc_info = args.ntc_name if ntc_col else "none"
    print(f"Wrote {out_path} (sample={args.sample_name}, ntc={ntc_info})")


if __name__ == "__main__":
    main()
