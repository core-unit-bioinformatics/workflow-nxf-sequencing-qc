#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) Core Unit Bioinformatics (CUBI), Medical Faculty, HHU.
# See the LICENSE file in the repository root.
"""Build the compact SequencingQC summary section for the top of the MultiQC report.

Collects per case (DRAGEN output prefix):

* TMB   from ``<case>.tmb.metrics.csv``
* MSI   from ``<case>.microsat_output.json``
* HRD   from ``<case>.hrdscore.csv``
* CNV QC status, purity and ploidy from the ``cnv-metrics.csv`` table of
  ``analyze_cnv_metrics.R`` (optional)
* run-level Picard CrosscheckFingerprints results (optional)

Biomarker parsing follows ``05-process-tsv-biomarker.py`` of the diagnostic
somatic-exome workflow, so the full biomarker table has the same columns as
its ``biomarker.csv``. One difference: HRD rows are keyed by the file prefix
(like TMB and MSI) instead of the ``Sample`` column inside the file, so all
three biomarkers of one DRAGEN run end up in one row.

Outputs:

* ``sequencingqc_summary_multiqc_config.yml``  MultiQC config setting ``intro_text``
  to the summary, so it is shown first (MultiQC always renders General
  Statistics before any custom-content section)
* ``sequencingqc_summary.html``  the same HTML fragment, standalone
* ``sequencingqc_summary.tsv``   visible summary columns, one row per case
* ``biomarker.tsv``         all biomarker values, one row per case

No existing MultiQC module reads DRAGEN TMB/MSI/HRD files or combines them
with CNV and fingerprint QC, hence this script.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import sys
from collections.abc import Iterable
from pathlib import Path

__version__ = "1.0.0"

NA = "n/a"
CONFIG_NAME = "sequencingqc_summary_multiqc_config.yml"
TMB_SUFFIX = ".tmb.metrics.csv"
MSI_SUFFIX = ".microsat_output.json"
HRD_SUFFIX = ".hrdscore.csv"

# DRAGEN keys shown in the compact table
TMB_KEY = "TMB"
MSI_KEY = "PercentageUnstableSites"
MSI_VALID_KEY = "ResultIsValid"
HRD_KEY = "HRD_sum"
# analyze_cnv_metrics.R columns
CNV_SAMPLE = "sample"
CNV_STATUS = "QC_STATUS"
CNV_REASONS = "QC_FAIL_REASONS"
CNV_PURITY = "EstimatedTumorPurity"
CNV_PLOIDY = "OverallPloidy"

FP_EXPECTED = {"EXPECTED_MATCH", "EXPECTED_MISMATCH"}

Table = dict[str, dict[str, str]]


class InputError(Exception):
    """Invalid or inconsistent input; message names file, value and fix."""


# --------------------------------------------------------------------------
# File discovery
# --------------------------------------------------------------------------
def find_files(roots: Iterable[Path], suffix: str) -> dict[str, Path]:
    """Map case -> file for all files ending in ``suffix`` below ``roots``.

    Symlinks are followed (Nextflow stages inputs as symlinks). The same case
    found in two different files is an error.
    """
    found: dict[str, Path] = {}
    for root in roots:
        if not root.is_dir():
            raise InputError(
                f"DRAGEN directory does not exist or is not a directory: {root}. "
                "Check --dragen_dirs."
            )
        for dirpath, _dirs, files in os.walk(root, followlinks=True):
            for name in sorted(files):
                if not name.endswith(suffix):
                    continue
                case = name[: -len(suffix)]
                path = Path(dirpath) / name
                if case in found and found[case].resolve() != path.resolve():
                    raise InputError(
                        f"Case '{case}' found twice for '*{suffix}': "
                        f"{found[case]} and {path}. Each DRAGEN output prefix must "
                        "be unique across --dragen_dirs."
                    )
                found[case] = path
    return found


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------
def read_tmb(path: Path) -> dict[str, str]:
    """DRAGEN ``*.tmb.metrics.csv``: ``SECTION,,metric,value`` lines."""
    metrics: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.reader(handle):
            if len(row) < 4:
                continue
            metrics[row[2]] = row[3]
    if TMB_KEY not in metrics:
        raise InputError(
            f"TMB file {path} has no '{TMB_KEY}' metric (expected a line like "
            f"'TMB SUMMARY,,{TMB_KEY},<value>'). Is this a DRAGEN TMB metrics file?"
        )
    return metrics


def flatten(obj, prefix: str = "") -> dict[str, str]:
    """Flatten nested JSON like R ``unlist()`` (as in the diag script)."""
    flat: dict[str, str] = {}
    if isinstance(obj, dict):
        for key, value in obj.items():
            flat.update(flatten(value, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(obj, list):
        for idx, value in enumerate(obj, start=1):
            flat.update(flatten(value, f"{prefix}{idx}"))
    elif isinstance(obj, bool):
        flat[prefix] = "TRUE" if obj else "FALSE"
    else:
        flat[prefix] = "" if obj is None else str(obj)
    return flat


def read_msi(path: Path) -> dict[str, str]:
    """DRAGEN ``*.microsat_output.json`` without the ``Settings`` block."""
    try:
        with path.open() as handle:
            data = json.load(handle)
    except json.JSONDecodeError as err:
        raise InputError(f"MSI file {path} is not valid JSON: {err}") from err
    if not isinstance(data, dict) or MSI_KEY not in data:
        raise InputError(
            f"MSI file {path} has no '{MSI_KEY}' field. "
            "Is this a DRAGEN microsat_output.json?"
        )
    data.pop("Settings", None)
    return flatten(data)


def read_hrd(path: Path) -> dict[str, str]:
    """DRAGEN ``*.hrdscore.csv``: header + exactly one data row."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise InputError(
            f"HRD file {path} has {len(rows)} data rows, expected exactly 1."
        )
    row = rows[0]
    if HRD_KEY not in row:
        raise InputError(
            f"HRD file {path} has no '{HRD_KEY}' column "
            f"(observed: {', '.join(row)}). Is this a DRAGEN hrdscore.csv?"
        )
    return {k: v for k, v in row.items() if k != "Sample"}


def read_cnv(paths: Iterable[Path]) -> Table:
    """``cnv-metrics.csv`` from analyze_cnv_metrics.R; empty files are skipped."""
    cnv: Table = {}
    for path in paths:
        if path.stat().st_size == 0:
            log(f"CNV metrics file is empty, skipped: {path}")
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            missing = {CNV_SAMPLE, CNV_STATUS} - set(reader.fieldnames or [])
            if missing:
                raise InputError(
                    f"CNV metrics file {path} lacks column(s) {sorted(missing)}. "
                    "Expected output of analyze_cnv_metrics.R."
                )
            for row in reader:
                case = row[CNV_SAMPLE]
                if case in cnv:
                    raise InputError(
                        f"Case '{case}' appears twice in CNV metrics ({path})."
                    )
                cnv[case] = row
    return cnv


def read_crosscheck(paths: Iterable[Path]) -> list[dict[str, str]]:
    """Picard ``*.crosscheck_metrics.txt``; self-comparisons are dropped."""
    pairs: list[dict[str, str]] = []
    seen = set()
    for path in paths:
        with path.open() as handle:
            lines = [ln for ln in handle if ln.strip() and not ln.startswith("#")]
        if not lines:
            log(f"Crosscheck file has no data, skipped: {path}")
            continue
        reader = csv.DictReader(lines, delimiter="\t")
        need = {"LEFT_GROUP_VALUE", "RIGHT_GROUP_VALUE", "RESULT"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise InputError(
                f"Crosscheck file {path} lacks column(s) {sorted(missing)}. "
                "Expected Picard CrosscheckFingerprints output."
            )
        for row in reader:
            left = Path(row["LEFT_GROUP_VALUE"]).name
            right = Path(row["RIGHT_GROUP_VALUE"]).name
            if left == right:
                continue
            key = tuple(sorted((left, right)))
            if key in seen:
                continue
            seen.add(key)
            pairs.append(
                {
                    "left": key[0],
                    "right": key[1],
                    "result": row["RESULT"],
                    "lod": row.get("LOD_SCORE", ""),
                }
            )
    return pairs


# --------------------------------------------------------------------------
# Combination
# --------------------------------------------------------------------------
def combine(tmb: Table, msi: Table, hrd: Table) -> tuple[list[str], Table]:
    """Merge biomarkers per case; column order as in the diag biomarker.csv."""
    columns: list[str] = []
    for table in (tmb, msi, hrd):
        for row in table.values():
            for col in row:
                if col not in columns:
                    columns.append(col)
    merged: Table = {}
    for case in sorted(set(tmb) | set(msi) | set(hrd)):
        row: dict[str, str] = {}
        for table in (tmb, msi, hrd):
            row.update(table.get(case, {}))
        merged[case] = row
    return columns, merged


def to_float(value: str | None) -> float | None:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def pct(value: str | None) -> str:
    """Fraction (0.62) as whole percent (62%)."""
    num = to_float(value)
    return NA if num is None else f"{num * 100:.0f}%"


def fmt(value: str | None, digits: int) -> str:
    num = to_float(value)
    if num is None:
        return NA
    return f"{num:.{digits}f}"


def summary_rows(
    biomarker: Table, cnv: Table, thresholds: dict[str, float | None]
) -> list[dict[str, str]]:
    """One compact row per case with display values and flags."""
    rows = []
    for case in sorted(set(biomarker) | set(cnv)):
        bio = biomarker.get(case, {})
        cnv_row = cnv.get(case, {})
        msi_valid = bio.get(MSI_VALID_KEY, "").upper()
        msi = fmt(bio.get(MSI_KEY), 1)
        if msi != NA and msi_valid == "FALSE":
            msi = f"{msi} (invalid)"
        status = cnv_row.get(CNV_STATUS) or NA
        row = {
            "case": case,
            "cnv_qc": status,
            "purity": pct(cnv_row.get(CNV_PURITY)),
            "ploidy": fmt(cnv_row.get(CNV_PLOIDY), 2),
            "tmb": fmt(bio.get(TMB_KEY), 1),
            "msi": msi,
            "hrd": fmt(bio.get(HRD_KEY), 0),
            "cnv_reasons": cnv_row.get(CNV_REASONS, "") if status == "FAIL" else "",
        }
        for key, src in (("tmb", TMB_KEY), ("msi", MSI_KEY), ("hrd", HRD_KEY)):
            cut = thresholds.get(key)
            val = to_float(bio.get(src))
            invalid = key == "msi" and msi_valid == "FALSE"
            row[f"{key}_high"] = (
                "1"
                if cut is not None and val is not None and val >= cut and not invalid
                else ""
            )
        rows.append(row)
    return rows


def missing_inputs(
    cases: Iterable[str], sources: dict[str, Table]
) -> list[tuple[str, list[str]]]:
    out = []
    for case in sorted(cases):
        miss = [name for name, table in sources.items() if case not in table]
        if miss:
            out.append((case, miss))
    return out


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
CSS = """
<style>
#exsum{font-size:.9rem;margin-top:.5em}
#exsum .exs-title{margin:.2em 0 .3em}
#exsum .exs-chips{display:flex;flex-wrap:wrap;gap:.4em;margin:.2em 0 .6em}
#exsum .exs-chip{padding:.15em .6em;border-radius:1em;font-size:90%;
 border:1px solid rgba(128,128,128,.4)}
#exsum .exs-ok{background:rgba(40,167,69,.18)}
#exsum .exs-bad{background:rgba(220,53,69,.22)}
#exsum .exs-warn{background:rgba(255,193,7,.25)}
#exsum table{width:auto;margin-bottom:.4em}
#exsum td,#exsum th{padding:.15em .7em;white-space:nowrap}
#exsum td.num,#exsum th.num{text-align:right;font-variant-numeric:tabular-nums}
#exsum .exs-hi{font-weight:bold;color:#c0392b}
#exsum .exs-na{opacity:.5}
#exsum details{margin:.3em 0}
#exsum summary{cursor:pointer;opacity:.8}
#exsum .exs-scroll{overflow-x:auto;max-width:100%}
</style>
"""


def esc(value: str) -> str:
    return html.escape(str(value), quote=True)


def cell(value: str, cls: str = "", title: str = "") -> str:
    classes = [c for c in (cls, "exs-na" if value == NA else "") if c]
    attr = f' class="{" ".join(classes)}"' if classes else ""
    if title:
        attr += f' title="{esc(title)}"'
    return f"<td{attr}>{esc(value)}</td>"


def chip(text: str, state: str) -> str:
    return f'<span class="exs-chip exs-{state}">{esc(text)}</span>'


def render_html(
    rows: list[dict[str, str]],
    columns: list[str],
    biomarker: Table,
    pairs: list[dict[str, str]],
    missing: list[tuple[str, list[str]]],
    thresholds: dict[str, float | None],
    have_cnv: bool,
    have_fp: bool,
) -> str:
    parts = [CSS, '<div id="exsum">', '<h3 class="exs-title">Summary</h3>']

    # --- run-level chips ---------------------------------------------------
    chips = [chip(f"{len(rows)} case{'' if len(rows) == 1 else 's'}", "")]
    if have_cnv:
        n_fail = sum(r["cnv_qc"] == "FAIL" for r in rows)
        n_pass = sum(r["cnv_qc"] == "PASS" for r in rows)
        chips.append(
            chip(f"CNV QC {n_pass}/{n_pass + n_fail} pass", "bad" if n_fail else "ok")
        )
    if have_fp:
        unexp = [p for p in pairs if p["result"].startswith("UNEXPECTED")]
        inconc = [p for p in pairs if p["result"] == "INCONCLUSIVE"]
        state = "bad" if unexp else ("warn" if inconc else "ok")
        text = f"Fingerprint {len(unexp)} unexpected"
        if inconc:
            text += f", {len(inconc)} inconclusive"
        chips.append(chip(text, state))
    if missing:
        chips.append(chip(f"{len(missing)} incomplete", "warn"))
    parts.append(f'<div class="exs-chips">{"".join(chips)}</div>')

    # --- compact table -----------------------------------------------------
    def th(label: str, key: str = "", tip: str = "") -> str:
        cut = thresholds.get(key) if key else None
        if cut is not None:
            tip = f"{tip}; bold if ≥ {cut:g}" if tip else f"bold if ≥ {cut:g}"
        title = f' title="{esc(tip)}"' if tip else ""
        cls = ' class="num"' if key or label in ("Purity", "Ploidy") else ""
        return f"<th{cls}{title}>{esc(label)}</th>"

    head = [th("Case")]
    if have_cnv:
        head += [th("CNV QC"), th("Purity"), th("Ploidy")]
    head += [
        th("TMB", "tmb", "mut/Mb (DRAGEN TMB)"),
        th("MSI %", "msi", "% unstable MS sites (DRAGEN PercentageUnstableSites)"),
        th("HRD", "hrd", "HRD sum = LOH + TAI + LST (DRAGEN)"),
    ]
    body = []
    for r in rows:
        tds = [f"<td>{esc(r['case'])}</td>"]
        if have_cnv:
            state = {"PASS": "exs-ok", "FAIL": "exs-bad"}.get(r["cnv_qc"], "")
            tds.append(cell(r["cnv_qc"], state, r["cnv_reasons"]))
            tds.append(cell(r["purity"], "num"))
            tds.append(cell(r["ploidy"], "num"))
        for key in ("tmb", "msi", "hrd"):
            cls = "num exs-hi" if r[f"{key}_high"] else "num"
            tds.append(cell(r[key], cls))
        body.append(f"<tr>{''.join(tds)}</tr>")
    parts.append(
        '<div class="exs-scroll"><table class="table table-sm table-condensed '
        f'table-hover"><thead><tr>{"".join(head)}</tr></thead>'
        f'<tbody>{"".join(body)}</tbody></table></div>'
    )

    # --- hidden details ----------------------------------------------------
    def details(title: str, inner: str) -> str:
        return f"<details><summary>{esc(title)}</summary>{inner}</details>"

    def table(header: list[str], data: list[list[str]]) -> str:
        h = "".join(f"<th>{esc(x)}</th>" for x in header)
        b = "".join("<tr>" + "".join(cell(x) for x in r) + "</tr>" for r in data)
        return (
            '<div class="exs-scroll"><table class="table table-sm '
            f'table-condensed"><thead><tr>{h}</tr></thead><tbody>{b}'
            "</tbody></table></div>"
        )

    if biomarker:
        data = [
            [case] + [biomarker[case].get(c) or NA for c in columns]
            for case in sorted(biomarker)
        ]
        parts.append(details("Biomarkers: all values", table(["case"] + columns, data)))
    failed = [r for r in rows if r["cnv_qc"] == "FAIL"]
    if failed:
        parts.append(
            details(
                f"CNV QC failures ({len(failed)})",
                table(
                    ["case", "reasons"],
                    [[r["case"], r["cnv_reasons"].replace(";", "; ")] for r in failed],
                ),
            )
        )
    if have_fp:
        odd = [p for p in pairs if p["result"] not in FP_EXPECTED]
        if odd:
            parts.append(
                details(
                    f"Fingerprint: not as expected ({len(odd)} of {len(pairs)} pairs)",
                    table(
                        ["file 1", "file 2", "result", "LOD"],
                        [[p["left"], p["right"], p["result"], p["lod"]] for p in odd],
                    ),
                )
            )
    if missing:
        parts.append(
            details(
                f"Incomplete cases ({len(missing)})",
                table(
                    ["case", "missing"],
                    [[case, ", ".join(miss)] for case, miss in missing],
                ),
            )
        )
    parts.append("</div>")
    return "\n".join(parts) + "\n"


def write_tsv(path: Path, header: list[str], data: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(data)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[sequencingqc_summary] {msg}", file=sys.stderr)


def positive(value: str) -> float:
    num = float(value)
    if num < 0:
        raise argparse.ArgumentTypeError(f"threshold must be >= 0, got {value}")
    return num


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--dragen",
        nargs="+",
        type=Path,
        required=True,
        help="DRAGEN output directories (searched recursively)",
    )
    p.add_argument(
        "--cnv-metrics",
        nargs="*",
        type=Path,
        default=[],
        help="cnv-metrics.csv from analyze_cnv_metrics.R",
    )
    p.add_argument(
        "--crosscheck",
        nargs="*",
        type=Path,
        default=[],
        help="Picard CrosscheckFingerprints metrics",
    )
    p.add_argument(
        "--tmb-high",
        type=positive,
        default=None,
        help="Flag TMB >= value (mut/Mb). Default: no flag",
    )
    p.add_argument(
        "--msi-high",
        type=positive,
        default=None,
        help="Flag MSI PercentageUnstableSites >= value. Default: no flag",
    )
    p.add_argument(
        "--hrd-high",
        type=positive,
        default=None,
        help="Flag HRD_sum >= value. Default: no flag",
    )
    p.add_argument("--prefix", default="", help="Output file name prefix")
    p.add_argument("--version", action="version", version=__version__)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    thresholds = {"tmb": args.tmb_high, "msi": args.msi_high, "hrd": args.hrd_high}
    try:
        tmb = {c: read_tmb(p) for c, p in find_files(args.dragen, TMB_SUFFIX).items()}
        msi = {c: read_msi(p) for c, p in find_files(args.dragen, MSI_SUFFIX).items()}
        hrd = {c: read_hrd(p) for c, p in find_files(args.dragen, HRD_SUFFIX).items()}
        cnv = read_cnv(args.cnv_metrics)
        pairs = read_crosscheck(args.crosscheck)
    except InputError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1
    log(
        f"found TMB={len(tmb)} MSI={len(msi)} HRD={len(hrd)} CNV={len(cnv)} "
        f"fingerprint pairs={len(pairs)}"
    )

    columns, biomarker = combine(tmb, msi, hrd)
    rows = summary_rows(biomarker, cnv, thresholds)
    sources = {"TMB": tmb, "MSI": msi, "HRD": hrd}
    if args.cnv_metrics:
        sources["CNV"] = cnv
    missing = missing_inputs((r["case"] for r in rows), sources)
    for case, miss in missing:
        log(f"case '{case}' has no {', '.join(miss)} input; shown as {NA}")
    if not rows:
        log("no TMB/MSI/HRD/CNV inputs found; summary lists run-level checks only")

    prefix = args.prefix
    fragment = render_html(
        rows,
        columns,
        biomarker,
        pairs,
        missing,
        thresholds,
        have_cnv=bool(args.cnv_metrics),
        have_fp=bool(args.crosscheck),
    )
    Path(f"{prefix}sequencingqc_summary.html").write_text(fragment)
    # MultiQC renders intro_text right below the title, above General Statistics;
    # JSON is valid YAML and keeps the HTML string safely quoted.
    Path(f"{prefix}{CONFIG_NAME}").write_text(
        json.dumps({"intro_text": fragment}, ensure_ascii=False) + "\n"
    )
    keys = ["case", "cnv_qc", "purity", "ploidy", "tmb", "msi", "hrd"]
    write_tsv(
        Path(f"{prefix}sequencingqc_summary.tsv"),
        ["sample"] + keys[1:],
        [[r[k] for k in keys] for r in rows],
    )
    write_tsv(
        Path(f"{prefix}biomarker.tsv"),
        ["sample"] + columns,
        [
            [case] + [biomarker[case].get(c) or NA for c in columns]
            for case in sorted(biomarker)
        ],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
