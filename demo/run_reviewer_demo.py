"""Run the reviewer-facing adaptability and bounded-validation demo.

The demo orchestrates existing scientific modules. It does not generate new
structures for arbitrary chemistries and does not claim scientific transfer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable


RELEASE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = RELEASE_ROOT / "experimental_extensions" / "system_manifest.json"
DEFAULT_OUTPUT_ROOT = RELEASE_ROOT / "demo" / "outputs"

SUPPORTED_CONCLUSIONS = [
    "The workflow accepts a schema-checked system manifest instead of relying on a fixed row order.",
    "Supplied composition-energy tables can be checked for uniqueness and manifold completeness.",
    "Pairwise CEF interpolation can be evaluated by training fit, non-endmember LOOCV, and composition-family holdout.",
    "Matched supplied energy backends can be compared by error, rank correlation, top-k overlap, and ranking reversals.",
]

EXCLUDED_CONCLUSIONS = [
    "Accuracy or transferability for a new host, ordered prototype, or bonding class.",
    "Applicability to N-, O-, S-, or P-containing compounds without new models and independent validation.",
    "Universal accuracy of UMA, CHGNet, DFT, or any other energy backend.",
    "Global phase stability, kinetic accessibility, or synthesizability.",
]


class DemoError(RuntimeError):
    """Raised after a failed step has been recorded in the demo summary."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _display_path(path: Path, root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def _git_sha(root: Path) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() or None


def _portable_command(command: list[str], root: Path) -> list[str]:
    portable = []
    for token in command:
        if token == sys.executable:
            portable.append("<python>")
            continue
        path = Path(token)
        if path.is_absolute():
            try:
                portable.append(path.resolve().relative_to(root.resolve()).as_posix())
            except ValueError:
                portable.append(path.name)
        else:
            portable.append(token.replace("\\", "/"))
    return portable


def _run_step(
    label: str,
    command: list[str],
    root: Path,
    logs_dir: Path,
    env: dict[str, str] | None = None,
) -> dict[str, Any]:
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / f"{label}.stdout.log"
    stderr_path = logs_dir / f"{label}.stderr.log"
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    duration = time.perf_counter() - started
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    return {
        "id": label,
        "status": "passed" if completed.returncode == 0 else "failed",
        "return_code": completed.returncode,
        "duration_seconds": round(duration, 3),
        "command": _portable_command(command, root),
        "stdout_log": _display_path(stdout_path, root),
        "stderr_log": _display_path(stderr_path, root),
    }


def extract_bounded_evidence(results: dict[str, Any]) -> dict[str, Any]:
    """Extract reviewer-readable metrics without relying on fixed module IDs."""
    manifolds: list[dict[str, Any]] = []
    backend_comparisons: list[dict[str, Any]] = []
    for module_id, module in results.get("modules", {}).items():
        result = module.get("result", {})
        if module.get("kind") == "manifold_regression":
            manifolds.append(
                {
                    "module_id": module_id,
                    "system_id": result.get("system_id"),
                    "host": result.get("host"),
                    "prototype": result.get("prototype"),
                    "mixing_elements": result.get("mixing_elements"),
                    "mixing_site_count": result.get("mixing_site_count"),
                    "n_compositions": result.get("n_compositions"),
                    "design_rank": result.get("design_rank"),
                    "training_R2": result.get("training_R2"),
                    "training_metrics": result.get("training_metrics"),
                    "nonendmember_LOOCV_metrics": result.get(
                        "nonendmember_LOOCV_metrics"
                    ),
                    "group_holdout": result.get("group_holdout"),
                    "transferability_claim": result.get("transferability_claim"),
                }
            )
        elif module.get("kind") == "energy_backend_comparison":
            pairs = []
            for row in result.get("pairwise_results", []):
                pairs.append(
                    {
                        key: row.get(key)
                        for key in (
                            "subset_id",
                            "pair_id",
                            "n",
                            "spearman_rho",
                            "RMSE",
                            "MAE",
                            "top_k",
                            "top_k_overlap_count",
                            "top_k_jaccard",
                            "ranking_reversals",
                            "ranking_reversal_fraction",
                        )
                    }
                )
            backend_comparisons.append(
                {
                    "module_id": module_id,
                    "unit": result.get("unit"),
                    "backends": result.get("backends"),
                    "pairwise_results": pairs,
                    "interpretation": result.get("interpretation"),
                    "limitation": result.get("limitation"),
                    "transferability_claim": result.get("transferability_claim"),
                }
            )
    return {
        "overall_status": results.get("overall_status"),
        "manifold_validations": manifolds,
        "energy_backend_comparisons": backend_comparisons,
    }


def _write_markdown_summary(summary: dict[str, Any], path: Path) -> None:
    evidence = summary.get("bounded_validation", {})
    lines = [
        "# Reviewer demo summary",
        "",
        f"- Status: **{summary['status']}**",
        f"- Mode: `{summary['mode']}`",
        f"- System: `{summary['input_contract']['system_id']}`",
        f"- Manifest: `{summary['input_contract']['manifest']}`",
        f"- Generated: `{summary['provenance']['generated_at_utc']}`",
        "",
        "## What was executed",
        "",
    ]
    for step in summary.get("steps", []):
        lines.append(
            f"- `{step['id']}`: {step['status']} ({step['duration_seconds']:.3f} s)"
        )
    lines.extend(["", "## Bounded numerical evidence", ""])
    for item in evidence.get("manifold_validations", []):
        loo = item.get("nonendmember_LOOCV_metrics") or {}
        holdout = item.get("group_holdout") or {}
        lines.extend(
            [
                f"### Manifold: `{item.get('module_id')}`",
                "",
                f"- Composition rows: {item.get('n_compositions')}",
                f"- Mixing elements: {', '.join(item.get('mixing_elements') or [])}",
                f"- Training R2: {item.get('training_R2')}",
                f"- Non-endmember LOOCV RMSE: {loo.get('RMSE')}",
                f"- Group-holdout RMSE range: {holdout.get('RMSE_min')} to {holdout.get('RMSE_max')}",
                "",
            ]
        )
    for item in evidence.get("energy_backend_comparisons", []):
        lines.extend([f"### Energy backends: `{item.get('module_id')}`", ""])
        for pair in item.get("pairwise_results", []):
            lines.append(
                f"- {pair.get('subset_id')} / {pair.get('pair_id')}: "
                f"Spearman rho={pair.get('spearman_rho')}, "
                f"top-k overlap={pair.get('top_k_overlap_count')}/{pair.get('top_k')}, "
                f"rank reversals={pair.get('ranking_reversals')}"
            )
        lines.append("")
    lines.extend(["## Supported conclusions", ""])
    lines.extend(f"- {item}" for item in summary["supported_conclusions"])
    lines.extend(["", "## Not established", ""])
    lines.extend(f"- {item}" for item in summary["excluded_conclusions"])
    lines.extend(
        [
            "",
            "The HTML/PDF report under `bounded_validation/` is a static view of the",
            "computed JSON/CSV evidence; it is not the computational engine.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_summary(summary: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "demo_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    _write_markdown_summary(summary, output_dir / "demo_summary.md")


def _default_output_dir(mode: str, manifest_path: Path) -> Path:
    if manifest_path.resolve() == DEFAULT_MANIFEST.resolve():
        return DEFAULT_OUTPUT_ROOT / mode.replace("-", "_")
    manifest = _read_json(manifest_path)
    safe_id = "".join(
        char if char.isalnum() or char in "-_" else "_"
        for char in str(manifest.get("system_id", "custom"))
    )
    return DEFAULT_OUTPUT_ROOT / f"custom_{safe_id}"


def _copy_full_uma_artifacts(root: Path, output_dir: Path) -> dict[str, Any]:
    si_source = root / "SI_Figures" / "SI_Fig03_TripleConsensus"
    si_target = output_dir / "si_fig3"
    si_target.mkdir(parents=True, exist_ok=True)
    names = [
        "data_FigF_TripleConsensus_Data_uma_regen.csv",
        "data_FigF_TripleConsensus_Origin_uma_regen.csv",
        "data_FigF_TripleConsensus_Summary_uma_regen.csv",
        "preview_FigF_TripleConsensus_uma_regen.png",
        "validation_FigF_UMA_Rerun.json",
    ]
    for name in names:
        shutil.copy2(si_source / name, si_target / name)
    return _read_json(si_target / "validation_FigF_UMA_Rerun.json")


def run_demo(
    mode: str,
    manifest_path: Path = DEFAULT_MANIFEST,
    output_dir: Path | None = None,
    checkpoint: Path | None = None,
    device: str = "cpu",
    create_pdf: bool = True,
    release_root: Path = RELEASE_ROOT,
    step_runner: Callable[..., dict[str, Any]] = _run_step,
) -> dict[str, Any]:
    root = release_root.resolve()
    manifest_path = manifest_path.resolve()
    manifest = _read_json(manifest_path)
    output_dir = (
        output_dir.resolve()
        if output_dir is not None
        else _default_output_dir(mode, manifest_path).resolve()
    )
    logs_dir = output_dir / "logs"
    bounded_dir = output_dir / "bounded_validation"
    summary: dict[str, Any] = {
        "schema_version": 1,
        "demo_id": "reviewer_adaptability_and_bounded_validation",
        "mode": mode,
        "status": "running",
        "purpose": (
            "Demonstrate manifest-driven software adaptability and bounded numerical "
            "validation on supplied data; not scientific transferability."
        ),
        "provenance": {
            "execution": "computed_from_manifest_selected_supplied_tables",
            "generated_at_utc": _utc_now(),
            "python_executable": Path(sys.executable).name,
            "python_version": platform.python_version(),
            "git_sha": _git_sha(root),
        },
        "input_contract": {
            "manifest": _display_path(manifest_path, root),
            "manifest_sha256": _sha256(manifest_path),
            "system_id": manifest.get("system_id"),
            "scientific_scope": manifest.get("scientific_scope"),
            "configured_modules": [
                {"id": item.get("id"), "kind": item.get("kind")}
                for item in manifest.get("modules", [])
                if item.get("enabled")
            ],
            "energy_data_origin": "precomputed tables referenced by the manifest",
            "automatic_new_structure_generation": False,
            "automatic_arbitrary_material_scoring": False,
        },
        "claim_boundary": manifest.get("claim_boundary"),
        "supported_conclusions": SUPPORTED_CONCLUSIONS + [
            (
                "Bundled calculations and validation artifacts can be reproduced "
                "deterministically within the stated scope."
                if manifest_path == DEFAULT_MANIFEST.resolve()
                else "The custom supplied tables completed the same contract, validation, "
                "visualization, and reporting path."
            )
        ],
        "excluded_conclusions": EXCLUDED_CONCLUSIONS,
        "steps": [],
        "artifacts": {},
    }

    commands: list[tuple[str, list[str], dict[str, str] | None]] = [
        (
            "release_verifier",
            [sys.executable, "scripts/verify_release.py"],
            None,
        ),
        (
            "bounded_validation",
            [
                sys.executable,
                "experimental_extensions/run_validation.py",
                str(manifest_path),
                "--output-dir",
                str(bounded_dir),
            ]
            + ([] if create_pdf else ["--no-report-pdf"]),
            None,
        ),
    ]

    if mode == "full-uma":
        full_dir = output_dir / "uma_recomputation"
        panel_dir = full_dir / "panel_c"
        uma_env = os.environ.copy()
        uma_env["UMA_DEVICE"] = device
        if checkpoint is not None:
            uma_env["UMA_CHECKPOINT"] = str(checkpoint.resolve())
        panel_command = [
            sys.executable,
            "Panel_c_OrderedVsDisordered/script_FigC_Rerun_UMA.py",
            "--device",
            device,
            "--output",
            str(panel_dir / "data_FigC_Raw_UMA_Energies_regen.csv"),
            "--report",
            str(panel_dir / "validation_FigC_UMA_Rerun.json"),
        ]
        if checkpoint is not None:
            panel_command.extend(["--checkpoint", str(checkpoint.resolve())])
        commands.extend(
            [
                ("panel_c_uma_rerun", panel_command, uma_env),
                (
                    "si_fig3_uma_rerun",
                    [
                        sys.executable,
                        "SI_Figures/SI_Fig03_TripleConsensus/script_FigF_TripleConsensus.py",
                        "--rerun-uma",
                    ],
                    uma_env,
                ),
            ]
        )

    try:
        for label, command, env in commands:
            step = step_runner(label, command, root, logs_dir, env=env)
            summary["steps"].append(step)
            _write_summary(summary, output_dir)
            if step["status"] != "passed":
                raise DemoError(
                    f"Step {label!r} failed; see {step['stderr_log']} and {step['stdout_log']}"
                )

        validation_path = bounded_dir / "validation_results.json"
        validation = _read_json(validation_path)
        summary["bounded_validation"] = extract_bounded_evidence(validation)
        summary["artifacts"]["bounded_validation"] = {
            name: _display_path(bounded_dir / value, root)
            for name, value in validation.get("artifacts", {}).items()
            if value is not None
        }
        if mode == "full-uma":
            panel_report_path = (
                output_dir
                / "uma_recomputation"
                / "panel_c"
                / "validation_FigC_UMA_Rerun.json"
            )
            panel_report = _read_json(panel_report_path)
            si_report = _copy_full_uma_artifacts(root, output_dir / "uma_recomputation")
            summary["full_uma_validation"] = {
                "panel_c": panel_report,
                "si_fig3": si_report,
                "checkpoint_requirement": (
                    "External gated UMA-s-1p1 checkpoint; checkpoint weights are not redistributed."
                ),
                "scientific_transferability_claim": "none",
            }
            summary["artifacts"]["full_uma"] = {
                "panel_c_report": _display_path(panel_report_path, root),
                "si_fig3_report": _display_path(
                    output_dir
                    / "uma_recomputation"
                    / "si_fig3"
                    / "validation_FigF_UMA_Rerun.json",
                    root,
                ),
            }
        summary["status"] = "passed"
    except Exception as exc:
        summary["status"] = "failed"
        summary["error"] = str(exc)
        _write_summary(summary, output_dir)
        raise

    _write_summary(summary, output_dir)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("quick", "full-uma"), default="quick")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--checkpoint", type=Path, help="Optional path to uma-s-1p1.pt")
    parser.add_argument("--device", default=os.environ.get("UMA_DEVICE", "cpu"))
    parser.add_argument("--no-report-pdf", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary = run_demo(
            mode=args.mode,
            manifest_path=args.manifest,
            output_dir=args.output_dir,
            checkpoint=args.checkpoint,
            device=args.device,
            create_pdf=not args.no_report_pdf,
        )
    except Exception as exc:
        print(f"Reviewer demo failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    actual_output = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else _default_output_dir(args.mode, args.manifest)
    )
    print(f"Reviewer demo outputs: {actual_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
