#!/usr/bin/env python3
"""
- Parses a taxpasta TSV file with columns like:
    taxonomy_id, name, rank, lineage, <sample1>, <sample2>, ...
- Uses regex to extract sample name and molecule type (DNA/RNA) from cleaned sample names.
- Splits the report by molecule type.
- Compares each sample vs the NTC from the same wet-lab prep within the same subset.
- Writes one TSV per sample and molecule subset.

Expected cleaned sample name patterns:
  NTCdate-prep-RNA
  NTCdate-prep-DNA
  NTC-date-prep-RNA
  NTC-date-prep-DNA
  samplename-prep-RNA
  samplename-prep-DNA
  NTC-date-prepRNA
  samplename-prepDNA

Usage:
  python flag_taxpasta_report.py taxpasta_report.tsv
  python flag_taxpasta_report.py taxpasta_report.tsv --outdir results/taxpasta
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd

META_COLS = {"taxonomy_id", "name", "rank", "lineage"}
REQ_COLS = {"taxonomy_id", "name"}

SUFFIX_RE = re.compile(r"(DNA|RNA).*$", re.IGNORECASE)


def parse_args(args=None):
    description = "Split taxpasta report by molecule with NTC comparison."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("input", help="Taxpasta TSV file")
    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for flagged taxpasta TSV files (default: input file directory)",
    )
    return parser.parse_args(args)


def clean_name(name: str) -> str:
    return SUFFIX_RE.sub(r"\1", name)  # \1 = DNA or RNA


def safe_filename_part(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip("-_.")


def extract_prep_token(sample: str) -> str:
    # Infer the wet-lab prep from the last token before the molecule suffix,
    # e.g. NTC-260305-ELB -> ELB and sample-HLSAN -> HLSAN.
    tokens = [token for token in re.split(r"[-_]+", sample) if token]
    if len(tokens) < 2:
        return ""
    return tokens[-1].upper()


def parse_sample(name: str) -> dict:
    cleaned = clean_name(name)
    # Accept both ...-DNA / ...-RNA and compact forms such as ...-HLSANDNA.
    m = re.match(r"^(?P<sample>.+?)(?:[-_])?(?P<molecule>DNA|RNA)$", cleaned, re.IGNORECASE)
    if not m:
        raise ValueError(
            f"Invalid sample name format: {cleaned}\n"
            f"Expected something like: sample-DNA, sample-ELB-RNA, or NTC-260305-HLSAN-RNA"
        )
    d = m.groupdict()
    d["molecule"] = d["molecule"].upper()
    d["prep"] = extract_prep_token(d["sample"])
    return d


def is_ntc(sample: str) -> bool:
    return sample.upper().startswith("NTC")


def per_sample_flag(sample_series: pd.Series, ntc_series: pd.Series) -> list[str]:
    out = []
    for v, n in zip(sample_series, ntc_series):
        if v == 0:
            out.append("in_NTC")
        elif n == 0:
            out.append("in_Sample")
        elif v > n:
            out.append("More")
        else:
            out.append("Less")
    return out


def build_ntc_reference(df_sub: pd.DataFrame, clean_ntc_cols: list[str]) -> pd.Series:
    if not clean_ntc_cols:
        return pd.Series(0, index=df_sub.index, dtype=int)
    return df_sub[clean_ntc_cols].sum(axis=1)


def load_taxpasta_report(path: str) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path, sep="\t")

    missing = REQ_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    for col in ("rank", "lineage"):
        if col not in df.columns:
            df[col] = ""

    sample_cols = [c for c in df.columns if c not in META_COLS]
    if not sample_cols:
        raise ValueError("No sample columns found")

    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df, sample_cols


def group_samples(columns: list[str]) -> dict:
    groups = defaultdict(lambda: {"all": [], "samples": [], "ntc": [], "meta": {}})
    for col in columns:
        meta = parse_sample(col)
        mol = meta["molecule"]
        # Keep the original column plus parsed metadata so downstream output
        # can be built per molecule and per sample.
        groups[mol]["meta"][col] = {**meta, "clean": clean_name(col)}
        groups[mol]["all"].append(col)
        groups[mol]["ntc" if is_ntc(meta["sample"]) else "samples"].append(col)
    return groups


def ntc_prep_order(cols: dict) -> list[str]:
    preps = []
    for ntc_col in cols["ntc"]:
        prep = cols["meta"][ntc_col]["prep"]
        if prep and prep not in preps:
            preps.append(prep)
    return preps


def process_subset(df: pd.DataFrame, cols: dict, output_path: str) -> None:
    all_cols, sample_cols, ntc_cols = cols["all"], cols["samples"], cols["ntc"]
    clean_sample_cols = [clean_name(c) for c in sample_cols]
    clean_ntc_cols = [clean_name(c) for c in ntc_cols]
    prep_order = ntc_prep_order(cols)

    df_sub = df.loc[
        df[all_cols].sum(axis=1) > 0,
        ["taxonomy_id", "name", "rank"] + sample_cols + ntc_cols + ["lineage"]
    ].copy()

    rename_map = {c: clean_name(c) for c in sample_cols + ntc_cols}
    df_sub.rename(columns=rename_map, inplace=True)

    compare_cols = clean_sample_cols + clean_ntc_cols
    if compare_cols:
        # Keep rows that have signal in at least one sample or NTC column.
        df_sub = df_sub.loc[df_sub[compare_cols].sum(axis=1) > 0].copy()

    for sample_col in sample_cols:
        clean_sample = cols["meta"][sample_col]["clean"]
        if prep_order:
            # Add one comparison column per available prep, e.g. _vsNTC-ELB.
            for prep in prep_order:
                matched_clean_ntc = [
                    cols["meta"][c]["clean"]
                    for c in cols["ntc"]
                    if cols["meta"][c]["prep"] == prep
                ]
                sample_ntc_reference = build_ntc_reference(df_sub, matched_clean_ntc)
                df_sub[f"{clean_sample}_vsNTC-{prep}"] = per_sample_flag(
                    df_sub[clean_sample], sample_ntc_reference
                )
        # If there are no NTC columns, the sample file is still written, just
        # without any vsNTC comparison columns.

    ordered = ["taxonomy_id", "name", "rank", "lineage"]
    for c in clean_sample_cols:
        ordered.append(c)
        if prep_order:
            for prep in prep_order:
                compare_col = f"{c}_vsNTC-{prep}"
                if compare_col in df_sub.columns:
                    ordered.append(compare_col)
    ordered.extend(clean_ntc_cols)

    df_sub.to_csv(output_path, sep="\t", index=False, columns=ordered)
    print(f"Wrote output -> {output_path}")


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.outdir) if args.outdir else input_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    df, sample_cols = load_taxpasta_report(args.input)
    groups = group_samples(sample_cols)

    for molecule, cols in sorted(groups.items()):
        print(f"{molecule}: {len(cols['samples'])} samples, {len(cols['ntc'])} NTC")
        for sample_col in cols["samples"]:
            sample_name = cols["meta"][sample_col]["clean"]
            output = output_dir / f"{input_path.stem}_{safe_filename_part(sample_name)}.tsv"
            # Write one output file per sample within each molecule subset.
            sample_cols = {
                "all": [sample_col] + cols["ntc"],
                "samples": [sample_col],
                "ntc": cols["ntc"],
                "meta": cols["meta"],
            }
            process_subset(df, sample_cols, str(output))


if __name__ == "__main__":
    main()
