#!/usr/bin/env python3
"""Validate the Stage 1 paper/sample metadata without downloading or modifying data.

Usage: python3 stages/01-paper-and-samples/validate.py [path/to/manifest.tsv]
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

DEFAULT_TSV = Path(__file__).resolve().parent / "data/GSE89223_samples.tsv"
REQUIRED = {
    "gsm", "title", "patient_id", "diagnosis", "analysis_group", "pair_status",
    "paper_final_set", "paper_group", "paper_exclusion_reason", "srx", "run",
    "biosample", "library_layout", "platform", "model", "spots", "review_flag",
}
EXPECTED_BIOLOGY = {
    "pca_tumor": 15,
    "pca_adjacent_normal": 14,
    "bph_adenoma": 2,
    "bph_normal": 1,
}
EXPECTED_FINAL_TUMOR = {
    "BP3", "CP1", "CP2", "CP3", "CP8", "CP11", "CP12", "CP13",
    "CP14", "CP15",
}
EXPECTED_FINAL_CONTROL = {
    "BN3", "CN1", "CN3", "CN8", "CN11", "CN12", "CN13",
    "CN14", "CN15", "BN1", "BP1", "BP2",
}
EXPECTED_EXCLUDED = {
    "CP5", "CN5", "CP6", "CN6", "CP7", "CN7",
    "CP9", "CN9", "CP10", "CN10",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate(path: Path) -> dict[str, int]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        names = reader.fieldnames or []
        require(len(names) == len(set(names)), "Duplicate TSV columns")
        require(REQUIRED <= set(names), f"Missing columns: {sorted(REQUIRED - set(names))}")
        require("paired_pca_sensitivity" not in names,
                "Project-specific downstream design must not enter Stage 1")
        rows = list(reader)
    require(len(rows) == 32, f"Expected 32 sequenced Runs, got {len(rows)}")

    for column, pattern in [
        ("gsm", r"GSM[0-9]+"),
        ("srx", r"SRX[0-9]+"),
        ("run", r"SRR[0-9]+"),
        ("biosample", r"SAMN[0-9]+"),
    ]:
        values = []
        for line, row in enumerate(rows, start=2):
            require(None not in row and None not in row.values(),
                    f"Malformed TSV at line {line}")
            require(bool(re.fullmatch(pattern, row[column])),
                    f"Invalid {column} in row {line}: {row[column]!r}")
            values.append(row[column])
        require(len(set(values)) == len(values), f"Duplicate {column} identifier")

    titles = [row["title"] for row in rows]
    require(len(set(titles)) == 32, "Duplicate sample titles")
    require(Counter(row["analysis_group"] for row in rows) == EXPECTED_BIOLOGY,
            "Biological sample composition differs from the source study")
    require(Counter(row["paper_group"] for row in rows)
            == {"tumor": 10, "control": 12, "excluded": 10},
            "Published 10/12/10 study grouping changed")
    require({row["title"] for row in rows if row["paper_group"] == "tumor"}
            == EXPECTED_FINAL_TUMOR, "Published tumor membership changed")
    require({row["title"] for row in rows if row["paper_group"] == "control"}
            == EXPECTED_FINAL_CONTROL, "Published control membership changed")
    require({row["title"] for row in rows if row["paper_group"] == "excluded"}
            == EXPECTED_EXCLUDED, "Published excluded sample membership changed")

    patient_counts = Counter(row["patient_id"] for row in rows)
    require(len(patient_counts) == 17, "Expected exactly 17 sequenced patients")
    require(Counter(patient_counts.values()) == {2: 15, 1: 2},
            "Expected 15 patient pairs and 2 unpaired samples")
    for row in rows:
        title = row["title"]
        included = row["paper_group"] != "excluded"
        require(row["paper_final_set"] == ("yes" if included else "no"),
                f"Final-paper inclusion inconsistent: {title}")
        require(bool(row["paper_exclusion_reason"]) == (not included),
                f"Exclusion provenance inconsistent: {title}")
        expected_pair = "paired" if patient_counts[row["patient_id"]] == 2 else "unpaired"
        require(row["pair_status"] == expected_pair, f"Pair status inconsistent: {title}")
        require(row["library_layout"] == "SINGLE", f"Not a single-end study Run: {title}")
        require(row["platform"] == "ION_TORRENT" and row["model"] == "Ion Torrent Proton",
                f"Unexpected sequencing platform: {title}")
        require(row["spots"].isdecimal() and int(row["spots"]) > 0,
                f"Missing or invalid SRA spots: {title}")

    by_title = {row["title"]: row for row in rows}
    for title, group in [("BP3", "pca_tumor"), ("BN3", "pca_adjacent_normal")]:
        row = by_title[title]
        require(row["analysis_group"] == group and row["diagnosis"] == "prostate cancer",
                f"B3 corrected diagnosis misrepresented: {title}")
        require("B-prefix" in row["review_flag"],
                f"Missing diagnosis-correction provenance: {title}")

    require(by_title["CP2"]["pair_status"] == "unpaired"
            and by_title["BP2"]["pair_status"] == "unpaired",
            "Expected unpaired study samples changed")
    require(by_title["CP2"]["psa"] == "13.5"
            and "8.6" in by_title["CP2"]["review_flag"],
            "CP2 GEO-vs-paper PSA discrepancy was erased")

    return {
        "sequenced_runs": len(rows),
        "patients": len(patient_counts),
        "final_paper_samples": sum(x["paper_final_set"] == "yes" for x in rows),
        "excluded_from_paper_comparison": len(EXPECTED_EXCLUDED),
    }


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) > 1 or (args and args[0] in {"-h", "--help"}):
        print(__doc__)
        return 0 if args and args[0] in {"-h", "--help"} else 2
    path = Path(args[0]) if args else DEFAULT_TSV
    try:
        summary = validate(path)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        print(f"STAGE 1 FAIL: {exc}", file=sys.stderr)
        return 1
    print("STAGE 1 PASS: " + ", ".join(f"{key}={value}" for key, value in summary.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
