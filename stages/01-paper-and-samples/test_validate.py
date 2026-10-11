#!/usr/bin/env python3
"""Offline regression checks for the Stage 1 identity and paper-cohort contract."""

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("stage1_validator", HERE / "validate.py")
assert _spec and _spec.loader
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)


class Stage1ManifestTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.file = Path(self.tmp.name) / "samples.tsv"
        self.file.write_bytes(mod.DEFAULT_TSV.read_bytes())

    def rewrite(self, transform):
        with self.file.open(encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        transform(rows)
        with self.file.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_canonical_source_contract(self):
        self.assertEqual(mod.validate(self.file)["sequenced_runs"], 32)

    def test_duplicate_run_is_rejected(self):
        self.rewrite(lambda rows: rows[1].update(run=rows[0]["run"]))
        with self.assertRaisesRegex(ValueError, "Duplicate run"):
            mod.validate(self.file)

    def test_paper_group_relabel_is_rejected(self):
        self.rewrite(lambda rows: rows[0].update(paper_group="control"))
        with self.assertRaises(ValueError):
            mod.validate(self.file)

    def test_b3_diagnosis_correction_cannot_be_erased(self):
        self.rewrite(lambda rows: next(x for x in rows if x["title"] == "BP3").update(diagnosis="benign prostatic hyperplasia"))
        with self.assertRaisesRegex(ValueError, "B3 corrected diagnosis"):
            mod.validate(self.file)

    def test_cp2_psa_discrepancy_must_remain_visible(self):
        self.rewrite(lambda rows: next(x for x in rows if x["title"] == "CP2").update(review_flag=""))
        with self.assertRaisesRegex(ValueError, "PSA discrepancy"):
            mod.validate(self.file)

    def test_missing_run_is_rejected(self):
        self.rewrite(lambda rows: rows.pop())
        with self.assertRaisesRegex(ValueError, "Expected 32"):
            mod.validate(self.file)


if __name__ == "__main__":
    unittest.main()
