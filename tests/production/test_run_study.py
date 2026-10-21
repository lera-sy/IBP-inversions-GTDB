"""Synthetic study-runner tests: isolation, missingness and aggregation."""

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from ibp.module5.results import COMPARISON_FIELDS
from ibp.module6.catalog import SUMMARY_FIELDS
from ibp.module7.aggregate import aggregate, read_table
from ibp.module8.analyze import SPECIES_METRICS
from ibp.production import run_study as study
from tests.module7.test_aggregate import fixture


class ProductionStudyTests(unittest.TestCase):
    def test_aggregate_can_normalize_isolated_source_cohort_ids(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            a = fixture(folder, "species_a", cohort="study__species_a")
            b = fixture(folder, "species_b", cohort="study__species_b")
            output = folder / "combined"
            aggregate([a, b], output, "study", source_cohorts={
                "species_a": "study__species_a", "species_b": "study__species_b"
            })
            summary = read_table(output / "species_summary.tsv", SUMMARY_FIELDS)
            self.assertEqual({r["cohort_id"] for r in summary}, {"study"})
            # Neither input was modified. Source provenance remains auditable.
            self.assertEqual(read_table(a, SUMMARY_FIELDS)[0]["cohort_id"], "study__species_a")

    def test_missing_source_cohort_mapping_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            a = fixture(folder, "species_a", cohort="study__species_a")
            with self.assertRaises(ValueError):
                aggregate([a], folder / "out", "study", source_cohorts={})

    def test_failure_isolated_and_partial_study_aggregates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            repo = root / "repo"
            (repo / "data/processed").mkdir(parents=True)
            (repo / "config").mkdir()
            pd.DataFrame([
                {"Organism Taxonomic ID": "100", "Assembly Accession": "GCF_A"},
                {"Organism Taxonomic ID": "200", "Assembly Accession": "GCF_B"},
            ]).to_csv(repo / "data/processed/selected_assemblies.tsv", sep="\t", index=False)
            pd.DataFrame([
                {"taxid": "100", "species_id": "species_a"},
                {"taxid": "200", "species_id": "species_b"},
            ]).to_csv(repo / "data/processed/species.tsv", sep="\t", index=False)
            config = repo / "config/study.yaml"
            config.write_text(yaml.safe_dump({"cohort": "study", "species": ["species_a", "species_b"],
                                              "run_assembly_limit": 1}), encoding="utf-8")
            scratch = root / "scratch"

            class Result:
                def __init__(self, code):
                    self.returncode = code

            def fake_snakemake(command, cwd, env, stdout, stderr, check, timeout):
                arguments = yaml.safe_load(Path(command[command.index("--configfile") + 1]).read_text())
                name = arguments["species"][0]
                if name == "species_b":
                    return Result(3)
                _, summary, comparisons = study.paths_for(scratch, "study", name)
                fixture(summary.parent.parent.parent, name, cohort="study__species_a")
                comparisons.parent.mkdir(parents=True, exist_ok=True)
                with comparisons.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=COMPARISON_FIELDS, delimiter="\t")
                    writer.writeheader()
                    writer.writerow({"status": "success"})
                return Result(0)

            with patch.object(study.subprocess, "run", side_effect=fake_snakemake):
                rows = study.run_study(repo, scratch, config, 2, species_workers=2)
            self.assertEqual({r["species_id"]: r["status"] for r in rows}, {
                "species_a": "completed", "species_b": "failed"
            })
            study_dir = scratch / "IBP-studies/study"
            statuses = study.read_table(study_dir / "species_status.tsv", study.STATUS_FIELDS)
            self.assertEqual(len(statuses), 2)
            self.assertIn("code 3", next(r for r in statuses if r["species_id"] == "species_b")["error"])
            stats = study.read_table(study_dir / "module8/study/species_metrics.tsv",
                                     SPECIES_METRICS)
            self.assertEqual({r["species_id"] for r in stats}, {"species_a"})

    def test_unknown_comparisons_are_marked_partial(self):
        with tempfile.TemporaryDirectory() as temp:
            from unittest.mock import Mock

            root = Path(temp)
            fake_summaries = [
                {"caller": "leen", "candidate_families": "1"},
                {"caller": "syri", "candidate_families": "0"},
            ]
            counters = {"success": 1, "unresolved": 1, "unsupported": 0, "failed": 0}
            with (
                patch.object(study.subprocess, "run", return_value=Mock(returncode=0)),
                patch.object(study, "validate_catalogue", return_value=(fake_summaries, counters)),
            ):
                row, path = study.run_one(
                    root, root, root / "study", "study", {},
                    "species_a", "100", cores=1
                )
            self.assertEqual(row["status"], "partial")
            self.assertEqual(row["comparisons_unresolved"], "1")
            self.assertEqual(row["leen_candidate_families"], "1")
            self.assertEqual(path.name, "summary.tsv")

    def test_input_names_reject_path_traversal(self):
        with self.assertRaises(ValueError):
            study.validate_name("../../bad", "study")
