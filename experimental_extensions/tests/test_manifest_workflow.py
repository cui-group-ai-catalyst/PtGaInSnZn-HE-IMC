from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from experimental_extensions.contracts import load_manifest, validate_manifest
from experimental_extensions.run_manifold import run as run_manifold
from experimental_extensions.run_validation import run_validation


class GenericManifoldTest(unittest.TestCase):
    def test_three_elements_and_shuffled_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            rows = []
            site_count = 4
            for a_count in range(site_count + 1):
                for b_count in range(site_count - a_count + 1):
                    c_count = site_count - a_count - b_count
                    a, b, c = (
                        a_count / site_count,
                        b_count / site_count,
                        c_count / site_count,
                    )
                    energy = (
                        1.0 * a + 2.0 * b + 3.0 * c
                        + 0.8 * a * b - 0.4 * a * c + 0.6 * b * c
                    )
                    rows.append(
                        {
                            "A_count": a_count,
                            "B_count": b_count,
                            "C_count": c_count,
                            "energy": energy,
                        }
                    )
            data = pd.DataFrame(rows).sample(frac=1.0, random_state=17)
            data.to_csv(temp / "synthetic.csv", index=False)
            config = {
                "schema_version": 1,
                "system_id": "synthetic_three_element_shuffled",
                "status": "test-only",
                "host": "X",
                "prototype": "test",
                "mixing_elements": ["A", "B", "C"],
                "mixing_site_count": site_count,
                "count_columns": {
                    "A": "A_count", "B": "B_count", "C": "C_count"
                },
                "energy_column": "energy",
                "energy_unit": "arbitrary",
                "input_csv": "synthetic.csv",
                "output_dir": "unused",
                "require_complete_integer_manifold": True,
                "group_holdout_element": "A",
                "scientific_scope": "Synthetic interface test only.",
            }
            (temp / "config.json").write_text(
                json.dumps(config), encoding="utf-8"
            )
            summary = run_manifold(temp / "config.json", temp / "outputs")
        self.assertEqual(summary["n_compositions"], 15)
        self.assertEqual(summary["n_pair_parameters"], 3)
        self.assertEqual(summary["design_rank"], 3)
        self.assertAlmostEqual(summary["training_R2"], 1.0, places=12)
        self.assertLess(summary["nonendmember_LOOCV_metrics"]["RMSE"], 1e-12)


class ManifestWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(__file__).resolve().parents[1]
        self.manifest_path = self.root / "system_manifest.json"

    def test_manifest_rejects_nonmetal_transferability_overclaim(self):
        manifest = load_manifest(self.manifest_path)
        unsafe = copy.deepcopy(manifest)
        unsafe["claim_boundary"]["nonmetal_compound_transferability"] = "demonstrated"
        with self.assertRaises(ValueError):
            validate_manifest(unsafe)

    def test_unified_validation_outputs_bounded_evidence(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            results = run_validation(
                self.manifest_path, output_dir, create_pdf=False
            )
            saved = json.loads(
                (output_dir / "validation_results.json").read_text(encoding="utf-8")
            )
            html = (output_dir / "evidence_report.html").read_text(encoding="utf-8")
            reversals = pd.read_csv(
                output_dir
                / "uma_mp_dft_binary_rank"
                / "ranking_reversals.csv"
            )
            self.assertTrue((output_dir / "validation_evidence.png").exists())
            self.assertTrue((output_dir / "validation_evidence.pdf").exists())
        self.assertEqual(results["overall_status"], "passed")
        self.assertEqual(saved["scientific_transferability_claim"], "none")
        pairwise = saved["modules"]["uma_mp_dft_binary_rank"]["result"][
            "pairwise_results"
        ]
        all_pairs = {
            row["pair_id"]: row for row in pairwise if row["subset_id"] == "all_hosts"
        }
        self.assertAlmostEqual(
            all_pairs["mp_dft__uma_s_1p1"]["spearman_rho"], 0.8679, places=4
        )
        self.assertAlmostEqual(
            all_pairs["mp_dft__chgnet"]["spearman_rho"], 0.9500, places=4
        )
        self.assertGreater(len(reversals), 0)
        self.assertIn("does not predict synthesizability", html)
        self.assertIn("not_evaluated", html)

    def test_unified_entry_accepts_renamed_three_element_and_backend_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            rows = []
            site_count = 4
            for a_count in range(site_count + 1):
                for b_count in range(site_count - a_count + 1):
                    c_count = site_count - a_count - b_count
                    a, b, c = (
                        a_count / site_count,
                        b_count / site_count,
                        c_count / site_count,
                    )
                    rows.append(
                        {
                            "A_sites": a_count,
                            "B_sites": b_count,
                            "C_sites": c_count,
                            "custom_energy": (
                                a + 2.0 * b + 3.0 * c
                                + 0.8 * a * b - 0.4 * a * c + 0.6 * b * c
                            ),
                        }
                    )
            pd.DataFrame(rows).sample(frac=1.0, random_state=23).to_csv(
                temp / "custom_manifold.csv", index=False
            )
            pd.DataFrame(
                {
                    "Material": ["m1", "m2", "m3", "m4", "m5"],
                    "reference_energy": [5.0, 4.0, 3.0, 2.0, 1.0],
                    "candidate_energy": [5.1, 3.8, 3.2, 1.9, 1.2],
                }
            ).to_csv(temp / "custom_backends.csv", index=False)
            config = {
                "schema_version": 1,
                "system_id": "custom_three_element",
                "status": "test-only",
                "host": "X",
                "prototype": "custom_A3B",
                "mixing_elements": ["A", "B", "C"],
                "mixing_site_count": site_count,
                "count_columns": {
                    "A": "A_sites", "B": "B_sites", "C": "C_sites"
                },
                "energy_column": "custom_energy",
                "energy_unit": "arbitrary energy unit",
                "input_csv": "custom_manifold.csv",
                "output_dir": "unused",
                "require_complete_integer_manifold": True,
                "group_holdout_element": "B",
                "scientific_scope": "Synthetic interface substitution test only.",
            }
            (temp / "custom_config.json").write_text(
                json.dumps(config), encoding="utf-8"
            )
            manifest = load_manifest(self.manifest_path)
            manifest["system_id"] = "custom_manifest_ids_and_columns"
            manifest["scientific_scope"] = "Synthetic end-to-end adaptability test only."
            manifest["contract_catalog"] = str(self.root / "module_contracts.json")
            manifest["modules"] = [
                {
                    "id": "renamed_composition_module",
                    "kind": "manifold_regression",
                    "enabled": True,
                    "config_path": "custom_config.json",
                },
                {
                    "id": "renamed_backend_module",
                    "kind": "energy_backend_comparison",
                    "enabled": True,
                    "config": {
                        "input_csv": "custom_backends.csv",
                        "key_column": "Material",
                        "backends": [
                            {
                                "id": "reference_backend",
                                "label": "Reference backend",
                                "column": "reference_energy",
                            },
                            {
                                "id": "candidate_backend",
                                "label": "Candidate backend",
                                "column": "candidate_energy",
                            },
                        ],
                        "unit": "arbitrary energy unit",
                        "primary_metric": "spearman_rank_correlation",
                        "ranking_direction": "descending",
                        "top_k": 2,
                        "subsets": [{"id": "all_materials", "expected_n": 5}],
                        "interpretation": "Synthetic matched-column interface check.",
                        "limitation": "No scientific material claim.",
                    },
                },
            ]
            custom_manifest = temp / "system_manifest.json"
            custom_manifest.write_text(json.dumps(manifest), encoding="utf-8")
            output_dir = temp / "outputs"
            results = run_validation(custom_manifest, output_dir, create_pdf=False)
            html = (output_dir / "evidence_report.html").read_text(encoding="utf-8")

        self.assertEqual(results["overall_status"], "passed")
        custom_manifold = results["modules"]["renamed_composition_module"]["result"]
        custom_backends = results["modules"]["renamed_backend_module"]["result"]
        self.assertEqual(custom_manifold["mixing_elements"], ["A", "B", "C"])
        self.assertEqual(custom_manifold["n_compositions"], 15)
        self.assertEqual(custom_backends["key_column"], "Material")
        self.assertEqual(
            [item["id"] for item in custom_backends["backends"]],
            ["reference_backend", "candidate_backend"],
        )
        self.assertIn("renamed_composition_module", html)
        self.assertIn("Candidate backend", html)


if __name__ == "__main__":
    unittest.main()
