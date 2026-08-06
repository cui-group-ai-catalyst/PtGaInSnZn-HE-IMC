"""Run all bounded P1 checks and create static reviewer-facing evidence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from experimental_extensions.compare_energy_backends import run_comparison
    from experimental_extensions.contracts import (
        load_contract_catalog,
        load_manifest,
        resolve_relative,
    )
    from experimental_extensions.reporting import write_html, write_pdf
    from experimental_extensions.run_manifold import run as run_manifold
    from experimental_extensions.visualization import write_validation_figure
except ModuleNotFoundError:
    from compare_energy_backends import run_comparison
    from contracts import load_contract_catalog, load_manifest, resolve_relative
    from reporting import write_html, write_pdf
    from run_manifold import run as run_manifold
    from visualization import write_validation_figure


def run_validation(
    manifest_path: Path,
    output_dir_override: Path | None = None,
    create_pdf: bool = True,
) -> dict:
    manifest_path = manifest_path.resolve()
    manifest = load_manifest(manifest_path)
    contract_catalog = load_contract_catalog(manifest_path, manifest)
    output_dir = (
        output_dir_override.resolve()
        if output_dir_override is not None
        else resolve_relative(manifest_path, manifest["output_dir"])
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    module_results: dict[str, dict] = {}
    for module in manifest["modules"]:
        if not module["enabled"]:
            continue
        module_dir = output_dir / module["id"]
        if module["kind"] == "manifold_regression":
            result = run_manifold(
                resolve_relative(manifest_path, module["config_path"]), module_dir
            )
        elif module["kind"] == "energy_backend_comparison":
            result = run_comparison(module, manifest_path, module_dir)
        else:
            raise AssertionError(f"Unhandled module kind: {module['kind']}")
        module_results[module["id"]] = {
            "kind": module["kind"],
            "status": "passed",
            "result": result,
        }

    results = {
        "schema_version": 1,
        "system_id": manifest["system_id"],
        "scientific_scope": manifest["scientific_scope"],
        "overall_status": "passed",
        "claim_boundary": manifest["claim_boundary"],
        "modules": module_results,
        "scientific_transferability_claim": "none",
        "artifacts": {
            "machine_readable": "validation_results.json",
            "si_figure_png": "validation_evidence.png",
            "si_figure_pdf": "validation_evidence.pdf",
            "reviewer_report_html": "evidence_report.html",
            "reviewer_report_pdf": "evidence_report.pdf" if create_pdf else None,
        },
    }
    (output_dir / "validation_results.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    write_validation_figure(results, output_dir / "validation_evidence")
    write_html(manifest, results, contract_catalog, output_dir / "evidence_report.html")
    if create_pdf:
        write_pdf(manifest, results, output_dir / "evidence_report.pdf")
    print(json.dumps(results, indent=2, ensure_ascii=True))
    print(f"Evidence outputs: {output_dir}")
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "manifest",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().with_name("system_manifest.json"),
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-report-pdf", action="store_true")
    args = parser.parse_args()
    run_validation(
        args.manifest,
        output_dir_override=args.output_dir,
        create_pdf=not args.no_report_pdf,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
