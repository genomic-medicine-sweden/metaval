#!/usr/bin/env python3
"""Generate a self-contained static HTML report from metaval pipeline outputs."""

import argparse
import base64
import csv
import gzip
import io
import json
import re
import zipfile
from datetime import date
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

SCRIPT_DIR = Path(__file__).resolve().parent
ASSETS_DIR = SCRIPT_DIR.parent / "assets"
STATIC_REPORT_DIR = ASSETS_DIR / "static_report"
REPORT_CSS = STATIC_REPORT_DIR / "report.css"
REPORT_JS = STATIC_REPORT_DIR / "report.js"
REPORT_TEMPLATE = STATIC_REPORT_DIR/ "report.html.j2"
FOOTER_LOGO = ASSETS_DIR / "metaval_logo_light.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a static HTML metaval report from a samplesheet.")
    parser.add_argument("--samplesheet", required=True, help="Path to the input samplesheet CSV.")
    parser.add_argument("--ticket", default="NA", help="Ticket number to show in the report header.")
    parser.add_argument("--version", required=True, help="Report or workflow version to show in the report header.")
    parser.add_argument("--output", required=True, help="Path to the output HTML report.")
    parser.add_argument("--flagged-dir", required=True, help="Path to flagged taxpasta TSV files.")
    parser.add_argument("--reads-dir", required=True, help="Path to extracted read files grouped by classifier.")
    parser.add_argument("--blastn-dir", required=True, help="Path to BLASTN result files.")
    parser.add_argument("--blastx-dir", required=True, help="Path to BLASTX result files.")
    parser.add_argument("--coverage-dir", required=True, help="Path to samtools coverage files.")
    parser.add_argument("--coverage-plots-dir", required=True, help="Path to coverage plot image files.")

    return parser.parse_args()

REPORT_TITLE = "Clinical Metagenomics Report"
CLASSIFIERS = ("kraken2", "centrifuge", "diamond")
CLASSIFIER_LABELS = {
    "kraken2": "Kraken2",
    "centrifuge": "Centrifuge",
    "diamond": "Diamond",
}

HEADER_HELP = {
    "qseqid": "Query sequence identifier.",
    "sseqid": "Subject or reference sequence identifier.",
    "slen": "Length of the subject or reference sequence.",
    "pident": "Percentage identity across the aligned region.",
    "qlen": "Length of the query sequence.",
    "length": "Alignment length.",
    "qcovs": "Query coverage per subject, as a percentage.",
    "nident": "Number of identical matches in the alignment.",
    "evalue": "Expected number of chance matches with this score or better.",
    "bitscore": "BLAST bit score for the alignment.",
    "staxid": "NCBI taxonomy identifier for the subject hit.",
    "ssciname": "Scientific name of the subject hit.",
    "#rname": "Reference sequence name.",
    "startpos": "Start position on the reference.",
    "endpos": "End position on the reference.",
    "numreads": "Number of reads mapped to the reference region.",
    "covbases": "Number of reference bases covered by reads.",
    "coverage": "Percentage of the reference covered by reads.",
    "meandepth": "Mean read depth across covered bases.",
    "meanbaseq": "Mean base quality of mapped reads.",
    "meanmapq": "Mean mapping quality of mapped reads.",
}


def file_to_data_uri(path: Path) -> str:
    """Encode a .png file as a data URI for embedding in the report."""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def read_samplesheet(samplesheet_path: Path) -> list[dict[str, str]]:
    """Read the samplesheet"""
    with samplesheet_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def infer_rows_from_flagged_dir(flagged_dir: Path) -> list[dict[str, str]]:
    """Infer sample rows from flagged taxpasta filenames when the samplesheet is empty."""
    samples = sorted({
        file_path.stem.rsplit("_", 1)[0]
        for file_path in flagged_dir.glob("*.tsv")
        if file_path.stem.rsplit("_", 1)[-1] in {"kraken2", "centrifuge", "diamond"}
    })
    return [
        {
            "sample": sample,
            "instrument_platform": "",
            "library_type": sample.rsplit("-", 1)[-1] if "-" in sample else "",
            "is_ntc": "false",
            "batch": "",
        }
        for sample in samples
    ]


def ntc_badge(value: str) -> str:
    """Return the CSS class used to style the sample/NTC badge."""
    normalized = value.strip().lower()
    return "ntc-true" if normalized == "true" else "ntc-false"


def ntc_label(value: str) -> str:
    """Return the display label shown inside the sample/NTC badge."""
    normalized = value.strip().lower()
    return "NTC" if normalized == "true" else "Sample"


def detect_superkingdom(name: str, rank: str, lineage: str) -> str:
    text = " ".join([name, rank, lineage]).lower()
    if "virus" in text or "viruses" in text:
        return "virus"
    if "bacteria" in text:
        return "bacteria"
    if "archaea" in text:
        return "archaea"
    if "fungi" in text:
        return "fungi"
    if "eukary" in text or "homo sapiens" in text:
        return "eukaryote"
    return "other"


def parse_read_count(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def load_assigned_reads(flagged_path: Path, sample_name: str, taxid: str) -> str:
    """Return the formatted read count assigned to a taxid in a flagged table."""
    if not flagged_path.exists():
        return ""

    with flagged_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if not rows or sample_name not in rows[0] or "taxonomy_id" not in rows[0]:
        return ""

    headers = rows[0]
    sample_index = headers.index(sample_name)
    taxid_index = headers.index("taxonomy_id")

    for row in rows[1:]:
        if sample_index < len(row) and taxid_index < len(row) and row[taxid_index] == taxid:
            return f"{parse_read_count(row[sample_index]):,}"

    return ""


def detect_flag(row: list[str]) -> str:
    """Return the first sample-vs-NTC flag found in a taxpasta row."""
    valid_flags = {"in_sample", "in_NTC", "> NTC", "< NTC", "equal"}
    for value in row:
        if value in valid_flags:
            return value
    return ""


def load_extracted_read_index(reads_dir: Path) -> dict[str, dict[str, list[dict[str, str]]]]:
    """Index extracted read and assembly files by sample and classifier."""
    index: dict[str, dict[str, list[dict[str, str]]]] = {}
    extracted_pattern = re.compile(
        r"^(?P<sample>.+?)_taxid_(?P<taxid>\d+)_(?P<organism>.+?)"
        r"\.extracted_(?P<classifier>kraken2|centrifuge|diamond)_"
    )
    assembly_pattern = re.compile(
        r"^(?P<sample>.+?)_taxid_(?P<taxid>\d+)_(?P<organism>.+?)_"
        r"(?P<classifier>kraken2|centrifuge|diamond)\.(?:scaffolds|contigs)\.fa(?:sta)?$"
    )

    if not reads_dir.exists():
        return index

    for file_path in reads_dir.iterdir():
        if not file_path.is_file():
            continue

        match = extracted_pattern.match(file_path.name) or assembly_pattern.match(file_path.name)
        if not match:
            continue

        sample = match.group("sample")
        taxid = match.group("taxid")
        organism = match.group("organism")
        classifier = match.group("classifier")
        index.setdefault(sample, {}).setdefault(classifier, [])
        entry = {"taxid": taxid, "organism": organism}
        if entry not in index[sample][classifier]:
            index[sample][classifier].append(entry)

    return index


def resolve_flagged_path(row: dict[str, str], flagged_dir: Path, classifier: str) -> Path:
    """Resolve the expected flagged taxpasta path for a sample and classifier."""
    candidates = [
        flagged_dir / f"{row['sample']}_{row['library_type']}_{row['batch']}_{classifier}.tsv",
        flagged_dir / f"{row['sample']}_{classifier}.tsv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def simplify_id(value: str) -> str:
    """Replaces characters that is not uppercase/lowercase letter, number,
    dot, underscore and hyphen with a single hyphen (-) """
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-").lower()


def read_tsv_table(path: Path) -> dict[str, list[list[str]]]:
    """Read a TSV file into report table headers and rows."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {"headers": [], "rows": []}
    return {"headers": rows[0], "rows": rows[1:]}


def parse_read_entries(path: Path, max_lines: int = 80) -> list[dict[str, str]]:
    """Parse FASTA or plain text read snippets for report display."""
    lines = path.read_text(encoding="utf-8").splitlines()[:max_lines]
    if not lines:
        return []

    entries: list[dict[str, str]] = []

    if lines[0].startswith(">"):
        current: list[str] = []
        for line in lines:
            if line.startswith(">") and current:
                entries.append({"header": current[0][1:], "content": "\n".join(current)})
                current = [line]
            else:
                current.append(line)
        if current:
            entries.append({"header": current[0][1:], "content": "\n".join(current)})
        return entries

    return [{"header": path.name, "content": "\n".join(lines)}]


def select_reads_source(
    reads_dir: Path,
    classifier: str,
    sample: str,
    taxid: str,
    organism: str,
) -> list[Path]:
    """Select read or assembly files matching a sample, taxid, organism, and classifier."""
    base_name = f"{sample}_taxid_{taxid}_{organism}"
    candidates = []
    candidates.extend(reads_dir.glob(f"{base_name}.extracted_{classifier}_*"))
    candidates.extend(reads_dir.glob(f"{base_name}_{classifier}.scaffolds.fa*"))
    candidates.extend(reads_dir.glob(f"{base_name}_{classifier}.contigs.fa*"))
    return sorted({path for path in candidates if path.is_file()})


def build_detail_context(
    sample: str,
    classifier: str,
    taxid: str,
    organism: str,
    assigned_reads: str,
    reads_files: list[Path],
    blastn_dir: Path,
    blastx_dir: Path,
    coverage_dir: Path,
    coverage_plots_dir: Path,
) -> dict[str, object]:
    """Build the lazy-loaded detail panel context for one extracted organism."""
    blastn_summary = blastn_dir / f"{sample}_taxid_{taxid}_{organism}_{classifier}_blast_filtered.txt"
    blastx_summary = blastx_dir / f"{sample}_taxid_{taxid}_{organism}_{classifier}_blastx_filtered.txt"

    blastn = read_tsv_table(blastn_summary) if blastn_summary.exists() else {"headers": [], "rows": []}
    blastx = read_tsv_table(blastx_summary) if blastx_summary.exists() else {"headers": [], "rows": []}

    coverage_files = sorted(coverage_dir.glob(f"{sample}_{classifier}_taxid_{taxid}_{organism}_mappingorganism_*.txt"))
    coverage_tables = [
        {"name": path.name, **read_tsv_table(path)}
        for path in coverage_files
    ]

    coverage_plot_files = sorted(
        path for path in coverage_plots_dir.glob(f"{sample}_{classifier}_taxid_{taxid}_{organism}_mappingorganism_*")
        if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg", ".webp"}
    )
    coverage_plot_sections = [
        {
            "name": path.name,
            "src": file_to_data_uri(path),
        }
        for path in coverage_plot_files
    ]

    reads_sections = [
        {
            "name": path.name,
            "entries": parse_read_entries(path),
        }
        for path in reads_files
        if "_read_2" not in path.name
    ]

    return {
        "sample": sample,
        "classifier": CLASSIFIER_LABELS[classifier],
        "taxid": taxid,
        "assigned_reads": assigned_reads,
        "organism": organism.replace("-", " "),
        "blastn": blastn,
        "blastx": blastx,
        "coverage_tables": coverage_tables,
        "coverage_plot_sections": coverage_plot_sections,
        "reads_sections": reads_sections,
    }


def collect_detail_sections(
    rows: list[dict[str, str]],
    flagged_dir: Path,
    extracted_reads_index: dict[str, dict[str, list[dict[str, str]]]],
    reads_dir: Path,
    blastn_dir: Path,
    blastx_dir: Path,
    coverage_dir: Path,
    coverage_plots_dir: Path,
) -> tuple[dict[tuple[str, str, str, str], str], dict[str, dict[str, object]]]:
    """Collect detail panel links and contexts for extracted organisms."""
    links: dict[tuple[str, str, str, str], str] = {}
    detail_sections: dict[str, dict[str, object]] = {}

    for row_index, row in enumerate(rows, start=1):
        sample = row["sample"]
        sample_id = f"sample-{row_index}"
        for classifier, entries in extracted_reads_index.get(sample, {}).items():
            for entry in entries:
                taxid = entry["taxid"]
                organism = entry["organism"]
                reads_files = select_reads_source(
                    reads_dir=reads_dir,
                    classifier=classifier,
                    sample=sample,
                    taxid=taxid,
                    organism=organism,
                )
                flagged_path = resolve_flagged_path(row, flagged_dir, classifier)
                panel_id = f"detail-{sample_id}-{classifier}-{taxid}-{simplify_id(organism)}"
                context = build_detail_context(
                    sample=sample,
                    classifier=classifier,
                    taxid=taxid,
                    organism=organism,
                    assigned_reads=load_assigned_reads(flagged_path, sample, taxid),
                    reads_files=reads_files,
                    blastn_dir=blastn_dir,
                    blastx_dir=blastx_dir,
                    coverage_dir=coverage_dir,
                    coverage_plots_dir=coverage_plots_dir,
                )
                detail_sections[panel_id] = {
                    "panel_id": panel_id,
                    "sample_id": sample_id,
                    "classifier_key": classifier,
                    "sample": context["sample"],
                    "classifier": context["classifier"],
                    "taxid": context["taxid"],
                    "assigned_reads": context["assigned_reads"],
                    "organism": context["organism"],
                    "blastn": context["blastn"],
                    "blastx": context["blastx"],
                    "coverage_tables": context["coverage_tables"],
                    "coverage_plot_sections": context["coverage_plot_sections"],
                    "reads_sections": context["reads_sections"],
                }
                links[(sample, classifier, taxid, organism)] = f"#{panel_id}"

    return links, detail_sections


def read_taxpasta_table(path: Path, sample_name: str) -> dict[str, list[list[str]]]:
    """Read a flagged taxpasta table and derive report display metadata."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {
            "headers": [],
            "rows": [],
            "lineages": [],
            "superkingdoms": [],
            "flags": [],
            "metaval_checked": [],
            "total_classified_reads": "",
            "total_unclassified_reads": "",
            "host_reads": "",
        }
    headers = rows[0]
    data_rows = rows[1:]
    lineages = [""] * len(data_rows)
    superkingdoms = ["other"] * len(data_rows)

    if "lineage" in headers:
        lineage_index = headers.index("lineage")
        reordered_indices = [
            index for index, header in enumerate(headers) if header != "lineage"
        ]
        headers = [headers[index] for index in reordered_indices]
        lineages = [
            row[lineage_index] if lineage_index < len(row) else ""
            for row in data_rows
        ]
        data_rows = [
            [row[index] if index < len(row) else "" for index in reordered_indices]
            for row in data_rows
        ]
    else:
        reordered_indices = list(range(len(headers)))

    name_index = headers.index("name") if "name" in headers else None
    rank_index = headers.index("rank") if "rank" in headers else None
    superkingdoms = [
        detect_superkingdom(
            row[name_index] if name_index is not None and name_index < len(row) else "",
            row[rank_index] if rank_index is not None and rank_index < len(row) else "",
            lineages[index],
        )
        for index, row in enumerate(data_rows)
    ]
    flags = [detect_flag(row) for row in data_rows]

    total_classified_reads = ""
    total_unclassified_reads = ""
    host_reads = ""
    sample_index = headers.index(sample_name) if sample_name in headers else None
    taxid_index = headers.index("taxonomy_id") if "taxonomy_id" in headers else None
    name_index = headers.index("name") if "name" in headers else None
    if sample_index is not None:
        total_count = sum(
            parse_read_count(row[sample_index])
            for row in data_rows
            if sample_index < len(row) and not (taxid_index is not None and taxid_index < len(row) and row[taxid_index] == '0')
        )
        host_count = sum(
            parse_read_count(row[sample_index])
            for row in data_rows
            if sample_index < len(row) and (
                (taxid_index is not None and taxid_index < len(row) and row[taxid_index] == '9606')
                or (name_index is not None and name_index < len(row) and row[name_index] == 'Homo sapiens')
            )
        )
        unclassified_count = sum(
            parse_read_count(row[sample_index])
            for row in data_rows
            if sample_index < len(row)
            and taxid_index is not None
            and taxid_index < len(row)
            and row[taxid_index] == '0'
        )
        total_classified_reads = f"{total_count:,}"
        total_unclassified_reads = f"{unclassified_count:,}"
        host_reads = f"{host_count:,}"

    return {
        "headers": headers,
        "rows": data_rows,
        "lineages": lineages,
        "superkingdoms": superkingdoms,
        "flags": flags,
        "metaval_checked": [False for _ in data_rows],
        "total_classified_reads": total_classified_reads,
        "total_unclassified_reads": total_unclassified_reads,
        "host_reads": host_reads,
    }


def add_cross_classifier_counts(classifiers: dict[str, dict[str, object]], sample_name: str) -> None:
    """Add Kraken2/Centrifuge/Diamond counts beside each classifier row."""
    count_maps = {key: {} for key in CLASSIFIERS}

    for key in CLASSIFIERS:
        table = classifiers[key]["table"]
        if not table or sample_name not in table["headers"] or "taxonomy_id" not in table["headers"]:
            continue
        sample_index = table["headers"].index(sample_name)
        taxid_index = table["headers"].index("taxonomy_id")
        count_maps[key] = {
            row[taxid_index]: row[sample_index]
            for row in table["rows"]
            if taxid_index < len(row) and sample_index < len(row)
        }

    for key in CLASSIFIERS:
        table = classifiers[key]["table"]
        if not table or sample_name not in table["headers"] or "taxonomy_id" not in table["headers"]:
            continue
        sample_index = table["headers"].index(sample_name)
        taxid_index = table["headers"].index("taxonomy_id")
        table["headers"] = (
            table["headers"][: sample_index + 1]
            + ["classifiers"]
            + table["headers"][sample_index + 1 :]
        )
        table["rows"] = [
            row[: sample_index + 1]
            + [[
                ("K", count_maps["kraken2"].get(row[taxid_index], "0")),
                ("C", count_maps["centrifuge"].get(row[taxid_index], "0")),
                ("D", count_maps["diamond"].get(row[taxid_index], "0")),
            ]]
            + row[sample_index + 1 :]
            for row in table["rows"]
            if taxid_index < len(row)
        ]


def prepare_rows(
    rows: list[dict[str, str]],
    flagged_dir: Path,
    extracted_reads_index: dict[str, dict[str, list[dict[str, str]]]],
    detail_links: dict[tuple[str, str, str, str], str],
) -> list[dict[str, str]]:
    """Prepare sample rows with classifier tables and links for rendering."""
    prepared_rows = []
    for index, row in enumerate(rows, start=1):
        prepared = dict(row)
        prepared["sample_id"] = f"sample-{index}"
        prepared["ntc_class"] = ntc_badge(row["is_ntc"])
        prepared["ntc_label"] = ntc_label(row["is_ntc"])
        prepared["classifiers"] = {}
        for classifier in CLASSIFIERS:
            file_path = resolve_flagged_path(row, flagged_dir, classifier)
            table = read_taxpasta_table(file_path, row["sample"]) if file_path.exists() else None
            extracted_matches = extracted_reads_index.get(row["sample"], {}).get(classifier, [])
            if table and "taxonomy_id" in table["headers"]:
                taxid_index = table["headers"].index("taxonomy_id")
                table_taxids = {
                    table_row[taxid_index]
                    for table_row in table["rows"]
                    if taxid_index < len(table_row)
                }
                seen_organisms = set()
                extracted_organisms = []
                panel_ids_by_taxid: dict[str, str] = {}
                for match in extracted_matches:
                    if match["taxid"] not in table_taxids:
                        continue
                    key = (match["taxid"], match["organism"])
                    if key in seen_organisms:
                        continue
                    seen_organisms.add(key)
                    href = detail_links.get((row["sample"], classifier, match["taxid"], match["organism"]), "")
                    extracted_organisms.append({
                        "name": match["organism"],
                        "href": href,
                        "panel_id": href[1:] if href.startswith("#") else "",
                    })
                    if href.startswith("#") and match["taxid"] not in panel_ids_by_taxid:
                        panel_ids_by_taxid[match["taxid"]] = href[1:]
                checked_taxids = {match["taxid"] for match in extracted_matches if match["taxid"] in table_taxids}
                table["metaval_checked"] = [
                    taxid_index < len(table_row) and table_row[taxid_index] in checked_taxids
                    for table_row in table["rows"]
                ]
                table["metaval_panel_ids"] = [
                    panel_ids_by_taxid.get(table_row[taxid_index], "")
                    if taxid_index < len(table_row)
                    else ""
                    for table_row in table["rows"]
                ]
            else:
                extracted_organisms = []
                if table is not None:
                    table["metaval_checked"] = [False for _ in table["rows"]]
                    table["metaval_panel_ids"] = ["" for _ in table["rows"]]
            prepared["classifiers"][classifier] = {
                "label": CLASSIFIER_LABELS[classifier],
                "table": table,
                "extracted_organisms": extracted_organisms,
            }
        add_cross_classifier_counts(prepared["classifiers"], row["sample"])
        prepared_rows.append(prepared)
    return prepared_rows


def build_html(
    ticket: str,
    version: str,
    rows: list[dict[str, str]],
    flagged_dir: Path,
    reads_dir: Path,
    blastn_dir: Path,
    blastx_dir: Path,
    coverage_dir: Path,
    coverage_plots_dir: Path,
) -> str:
    """Render the static report HTML and embed lazy-loaded detail data."""
    created_date = date.today().isoformat()
    report_css = REPORT_CSS.read_text(encoding="utf-8")
    report_js = REPORT_JS.read_text(encoding="utf-8")
    extracted_reads_index = load_extracted_read_index(reads_dir)
    detail_links, detail_sections = collect_detail_sections(
        rows=rows,
        flagged_dir=flagged_dir,
        extracted_reads_index=extracted_reads_index,
        reads_dir=reads_dir,
        blastn_dir=blastn_dir,
        blastx_dir=blastx_dir,
        coverage_dir=coverage_dir,
        coverage_plots_dir=coverage_plots_dir,
    )
    environment = Environment(
        loader=FileSystemLoader(REPORT_TEMPLATE.parent),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = environment.get_template(REPORT_TEMPLATE.name)
    archive_buffer = io.BytesIO()
    with zipfile.ZipFile(archive_buffer, "w", compression=zipfile.ZIP_STORED) as archive:
        for panel_id, detail in detail_sections.items():
            detail_json = json.dumps(
                detail,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
            archive.writestr(
                f"details/{panel_id}.json.gz",
                gzip.compress(detail_json),
            )
    detail_sections_zip_b64 = base64.b64encode(archive_buffer.getvalue()).decode("ascii")
    return template.render(
        report_title=REPORT_TITLE,
        header_help_map=HEADER_HELP,
        report_css=report_css,
        report_js=report_js,
        footer_logo=file_to_data_uri(FOOTER_LOGO),
        ticket=ticket,
        version=version,
        created_date=created_date,
        rows=prepare_rows(rows, flagged_dir, extracted_reads_index, detail_links),
        detail_sections_zip_b64=detail_sections_zip_b64,
    )


def main() -> None:
    args = parse_args()
    samplesheet_path = Path(args.samplesheet)
    output_path = Path(args.output)
    flagged_dir = Path(args.flagged_dir)
    reads_dir = Path(args.reads_dir)
    blastn_dir = Path(args.blastn_dir)
    blastx_dir = Path(args.blastx_dir)
    coverage_dir = Path(args.coverage_dir)
    coverage_plots_dir = Path(args.coverage_plots_dir)

    rows = read_samplesheet(samplesheet_path)
    if not rows:
        rows = infer_rows_from_flagged_dir(flagged_dir)
    report_html = build_html(
        args.ticket,
        args.version,
        rows,
        flagged_dir,
        reads_dir,
        blastn_dir,
        blastx_dir,
        coverage_dir,
        coverage_plots_dir,
    )
    output_path.write_text(report_html, encoding="utf-8")

    print(f"Wrote report to {output_path}")


if __name__ == "__main__":
    main()
