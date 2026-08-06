from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from demo.run_reviewer_demo import extract_bounded_evidence, run_demo


SAFE_BOUNDARY = {
    "software_configurability": "demonstrated_for_bundled_inputs",
    "within_manifold_interpolation": "validated_internal",
    "binary_energy_rank_consistency": "observed_reference_set",
    "new_host_or_prototype_transferability": "not_evaluated",
    "nonmetal_compound_transferability": "not_evaluated",
    "synthesizability": "not_predicted",
}


def sample_results() -> dict:
    return {
        "overall_status": "passed",
        "claim_boundary": SAFE_BOUNDARY,
        "modules": {
            "renamed_manifold": {
                "kind": "manifold_regression",
                "result": {
                    "system_id": "synthetic",
                    "host": "X",
                    "prototype": "A3B",
                    "mixing_elements": ["A", "B", "C"],
                    "mixing_site_count": 4,
                    "n_compositions": 15,
                    "design_rank": 3,
                    "training_R2": 1.0,
                    "training_metrics": {"RMSE": 0.0},
                    "nonendmember_LOOCV_metrics": {"RMSE": 0.0},
                    "group_holdout": {"RMSE_min": 0.0, "RMSE_max": 0.0},
                    "transferability_claim": "none",
                },
            },
            "renamed_backends": {
                "kind": "energy_backend_comparison",
                "result": {
                    "unit": "arbitrary",
                    "backends": [{"id": "a"}, {"id": "b"}],
                    "pairwise_results": [
                        {
                            "subset_id": "all",
                            "pair_id": "a__b",
                            "n": 5,
                            "spearman_rho": 0.9,
                            "RMSE": 1.0,
                            "MAE": 0.8,
                            "top_k": 2,
                            "top_k_overlap_count": 2,
                            "top_k_jaccard": 1.0,
                            "ranking_reversals": 1,
                            "ranking_reversal_fraction": 0.1,
                        }
                    ],
                    "transferability_claim": "none",
                },
            },
        },
        "artifacts": {
            "machine_readable": "validation_results.json",
            "reviewer_report_html": "evidence_report.html",
        },
    }


class EvidenceExtractionTest(unittest.TestCase):
    def test_module_ids_are_not_hard_coded(self):
        evidence = extract_bounded_evidence(sample_results())
        self.assertEqual(
            evidence["manifold_validations"][0]["module_id"], "renamed_manifold"
        )
        self.assertEqual(
            evidence["energy_backend_comparisons"][0]["module_id"],
            "renamed_backends",
        )
        self.assertEqual(
            evidence["energy_backend_comparisons"][0]["pairwise_results"][0][
                "spearman_rho"
            ],
            0.9,
        )

    def test_quick_orchestration_writes_auditable_summary(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = {
                "schema_version": 1,
                "system_id": "custom_test",
                "status": "bounded-validation-demo",
                "scientific_scope": "Test-only supplied data.",
                "output_dir": "unused",
                "contract_catalog": "contracts.json",
                "claim_boundary": SAFE_BOUNDARY,
                "modules": [
                    {
                        "id": "renamed_manifold",
                        "kind": "manifold_regression",
                        "enabled": True,
                        "config_path": "unused.json",
                    }
                ],
            }
            manifest_path = temp / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            output_dir = temp / "outputs"

            def fake_step(label, command, root, logs_dir, env=None):
                logs_dir.mkdir(parents=True, exist_ok=True)
                (logs_dir / f"{label}.stdout.log").write_text("ok", encoding="utf-8")
                (logs_dir / f"{label}.stderr.log").write_text("", encoding="utf-8")
                if label == "bounded_validation":
                    bounded_dir = output_dir / "bounded_validation"
                    bounded_dir.mkdir(parents=True, exist_ok=True)
                    (bounded_dir / "validation_results.json").write_text(
                        json.dumps(sample_results()), encoding="utf-8"
                    )
                return {
                    "id": label,
                    "status": "passed",
                    "return_code": 0,
                    "duration_seconds": 0.01,
                    "command": command,
                    "stdout_log": str(logs_dir / f"{label}.stdout.log"),
                    "stderr_log": str(logs_dir / f"{label}.stderr.log"),
                }

            summary = run_demo(
                "quick",
                manifest_path=manifest_path,
                output_dir=output_dir,
                release_root=temp,
                step_runner=fake_step,
            )
            saved = json.loads(
                (output_dir / "demo_summary.json").read_text(encoding="utf-8")
            )
            markdown = (output_dir / "demo_summary.md").read_text(encoding="utf-8")

        self.assertEqual(summary["status"], "passed")
        self.assertEqual(len(saved["steps"]), 2)
        self.assertFalse(saved["input_contract"]["automatic_arbitrary_material_scoring"])
        self.assertIn("Not established", markdown)

    def test_full_uma_summary_uses_fresh_isolated_reports(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = {
                "schema_version": 1,
                "system_id": "full_mode_test",
                "status": "bounded-validation-demo",
                "scientific_scope": "Test-only supplied data.",
                "output_dir": "unused",
                "contract_catalog": "contracts.json",
                "claim_boundary": SAFE_BOUNDARY,
                "modules": [
                    {
                        "id": "renamed_manifold",
                        "kind": "manifold_regression",
                        "enabled": True,
                        "config_path": "unused.json",
                    }
                ],
            }
            manifest_path = temp / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            output_dir = temp / "outputs"

            def fake_step(label, command, root, logs_dir, env=None):
                logs_dir.mkdir(parents=True, exist_ok=True)
                (logs_dir / f"{label}.stdout.log").write_text("ok", encoding="utf-8")
                (logs_dir / f"{label}.stderr.log").write_text("", encoding="utf-8")
                if label == "bounded_validation":
                    bounded_dir = output_dir / "bounded_validation"
                    bounded_dir.mkdir(parents=True, exist_ok=True)
                    (bounded_dir / "validation_results.json").write_text(
                        json.dumps(sample_results()), encoding="utf-8"
                    )
                elif label == "panel_c_uma_rerun":
                    report_path = Path(command[command.index("--report") + 1])
                    report_path.parent.mkdir(parents=True, exist_ok=True)
                    report_path.write_text(
                        json.dumps({"status": "passed", "marker": "fresh"}),
                        encoding="utf-8",
                    )
                elif label == "si_fig3_uma_rerun":
                    source = root / "SI_Figures" / "SI_Fig03_TripleConsensus"
                    source.mkdir(parents=True, exist_ok=True)
                    for name in (
                        "data_FigF_TripleConsensus_Data_uma_regen.csv",
                        "data_FigF_TripleConsensus_Origin_uma_regen.csv",
                        "data_FigF_TripleConsensus_Summary_uma_regen.csv",
                    ):
                        (source / name).write_text("fixture\n", encoding="utf-8")
                    (source / "preview_FigF_TripleConsensus_uma_regen.png").write_bytes(
                        b"fixture"
                    )
                    (source / "validation_FigF_UMA_Rerun.json").write_text(
                        json.dumps({"status": "passed"}), encoding="utf-8"
                    )
                return {
                    "id": label,
                    "status": "passed",
                    "return_code": 0,
                    "duration_seconds": 0.01,
                    "command": command,
                    "stdout_log": str(logs_dir / f"{label}.stdout.log"),
                    "stderr_log": str(logs_dir / f"{label}.stderr.log"),
                }

            summary = run_demo(
                "full-uma",
                manifest_path=manifest_path,
                output_dir=output_dir,
                release_root=temp,
                step_runner=fake_step,
            )

        self.assertEqual(summary["full_uma_validation"]["panel_c"]["marker"], "fresh")
        self.assertIn(
            "uma_recomputation/panel_c",
            summary["artifacts"]["full_uma"]["panel_c_report"].replace("\\", "/"),
        )


if __name__ == "__main__":
    unittest.main()
