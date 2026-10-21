"""Topology lookup must match exact sequence versions and preserve QC decisions."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import requests

from ibp.module3 import fetch_chromosome_topology as topology


class TopologyTest(unittest.TestCase):
    def test_matching_and_review_cases(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            qc = root / "qc.tsv"
            output = root / "topology.tsv"
            rows = [{"Assembly Accession": f"GCF_{i}.1", "Organism Taxonomic ID": "1",
                     "QC Status": "PASS" if i < 5 else "REVIEW",
                     "Chromosome Accession": f"CHR_{i}.1", "Chromosome Length": "100"}
                    for i in range(6)]
            pd.DataFrame(rows).to_csv(qc, sep="\t", index=False)
            records = [
                {"accessionversion": "CHR_0.1", "slen": 100, "topology": "circular"},
                {"accessionversion": "CHR_1.1", "slen": 100, "topology": "linear"},
                {"accessionversion": "CHR_2.1", "slen": 100, "topology": None},
                {"accessionversion": "CHR_3.1", "slen": 101, "topology": "circular"},
                # A newer version must not substitute for the selected chromosome.
                {"accessionversion": "CHR_4.2", "slen": 100, "topology": "circular"},
            ]
            response = {"uids": [str(i) for i in range(5)]}
            response.update({str(i): record for i, record in enumerate(records)})
            with (
                patch.object(topology, "fetch_summaries", return_value=response) as fetch,
                patch.object(sys, "argv", ["topology", "--qc", str(qc), "--output", str(output)]),
            ):
                topology.main()
            self.assertEqual(fetch.call_args.args[0], [f"CHR_{i}.1" for i in range(5)])
            result = pd.read_csv(output, sep="\t").fillna("")
            self.assertEqual(result["Topology Status"].tolist(), ["PASS"] * 2 + ["REVIEW"] * 3)
            self.assertEqual(result["Topology Notes"].tolist(), ["", "", "unknown_topology",
                                                               "chromosome_length_mismatch", "missing_exact_sequence_summary"])
            self.assertEqual(len(result), 5)

    def test_http_error_retries_with_pauses(self):
        failed = Mock()
        failed.raise_for_status.side_effect = requests.HTTPError("429")
        success = Mock()
        success.json.return_value = {"result": {"uids": []}}
        with (
            patch.object(topology.requests, "get", side_effect=[failed, success]) as get,
            patch.object(topology.time, "sleep") as sleep,
        ):
            self.assertEqual(topology.fetch_summaries(["CHR_0.1"]), {"uids": []})
        self.assertEqual(get.call_count, 2)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [1, 30])

    def test_api_error_inside_successful_http_response_is_not_accepted(self):
        response = Mock()
        response.json.return_value = {"error": "API rate limit exceeded"}
        with (
            patch.object(topology.requests, "get", return_value=response) as get,
            patch.object(topology.time, "sleep") as sleep,
            self.assertRaisesRegex(ValueError, "invalid sequence-summary"),
        ):
            topology.fetch_summaries(["CHR_0.1"])
        self.assertEqual(get.call_count, 3)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [1, 30, 60])


if __name__ == "__main__":
    unittest.main()
