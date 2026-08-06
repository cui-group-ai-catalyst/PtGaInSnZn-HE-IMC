"""Create deterministic reports from manifest-selected validation modules."""
from __future__ import annotations

from html import escape
from pathlib import Path
from textwrap import wrap


CLAIM_EXPLANATIONS = {
    "software_configurability": "Bundled modules are driven by external JSON inputs.",
    "within_manifold_interpolation": "Validation is internal to the supplied fixed composition manifold.",
    "binary_energy_rank_consistency": "Agreement is observed only on the supplied matched reference set.",
    "new_host_or_prototype_transferability": "No cross-host or cross-prototype validation is supplied.",
    "nonmetal_compound_transferability": "N/O/S/P systems require replacement or revalidation of scientific modules.",
    "synthesizability": "Continuous calculations prioritize candidates; they do not classify synthesis success.",
}


def _fmt(value: object, digits: int = 6) -> str:
    return f"{value:.{digits}f}" if isinstance(value, float) else str(value)


def _manifold_html(module_id: str, result: dict) -> str:
    rows: list[tuple[str, object]] = [
        ("Composition rows", result["n_compositions"]),
        ("Mixing elements", ", ".join(result["mixing_elements"])),
        ("Pair parameters", result["n_pair_parameters"]),
        ("Training R2", result["training_R2"]),
        ("Training RMSE", result["training_metrics"]["RMSE"]),
        ("Endmember-only RMSE", result["endmember_only_nonendmember_metrics"]["RMSE"]),
        ("Non-endmember LOOCV RMSE", result["nonendmember_LOOCV_metrics"]["RMSE"]),
        ("Energy unit", result.get("energy_unit", "not reported")),
    ]
    group = result.get("group_holdout")
    if group:
        rows.append(
            (
                "Group-holdout RMSE range",
                f"{group['RMSE_min']:.6f}-{group['RMSE_max']:.6f}",
            )
        )
    table = "".join(
        f"<tr><th>{escape(label)}</th><td>{escape(_fmt(value))}</td></tr>"
        for label, value in rows
    )
    return (
        f"<section><h3>{escape(module_id)} | Manifold regression</h3>"
        f"<table>{table}</table><p>{escape(result['scientific_scope'])}</p>"
        "<p>This is an internal interpolation and ablation check on the supplied "
        "energy table, not an external validation of its energy backend.</p></section>"
    )


def _comparison_html(module_id: str, result: dict) -> str:
    rows = "".join(
        "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            escape(row["subset_id"]), escape(row["pair_id"]), row["n"],
            _fmt(row["spearman_rho"]), _fmt(row["RMSE"]),
            row["top_k_overlap_count"], row["ranking_reversals"],
            _fmt(row["ranking_reversal_fraction"]),
        )
        for row in result["pairwise_results"]
    )
    labels = ", ".join(
        backend.get("label", backend["id"]) for backend in result["backends"]
    )
    return f"""<section><h3>{escape(module_id)} | Energy-backend comparison</h3>
<p><strong>Backends:</strong> {escape(labels)}; <strong>unit:</strong> {escape(result['unit'])}</p>
<table><tr><th>Subset</th><th>Pair</th><th>N</th><th>Spearman rho</th><th>RMSE</th><th>Top-k overlap</th><th>Reversals</th><th>Reversal fraction</th></tr>{rows}</table>
<p>{escape(result['interpretation'])}</p><p><strong>Limitation:</strong> {escape(result['limitation'])}</p></section>"""


def write_html(
    manifest: dict, results: dict, contract_catalog: dict, output_path: Path
) -> None:
    claims = "".join(
        "<tr><th>{}</th><td><code>{}</code></td><td>{}</td></tr>".format(
            escape(key), escape(value), escape(CLAIM_EXPLANATIONS[key])
        )
        for key, value in manifest["claim_boundary"].items()
    )
    module_sections = []
    for module_id, module in results["modules"].items():
        if module["kind"] == "manifold_regression":
            module_sections.append(_manifold_html(module_id, module["result"]))
        elif module["kind"] == "energy_backend_comparison":
            module_sections.append(_comparison_html(module_id, module["result"]))
    contracts = "".join(
        "<section><h3>{}</h3><p><strong>Evidence:</strong> {}<br><strong>Excluded:</strong> {}</p></section>".format(
            escape(item["kind"]), escape(item["evidence_level"]),
            escape("; ".join(item["excluded_inferences"])),
        )
        for item in contract_catalog["contracts"]
    )
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Bounded validation report</title><style>
body{{font-family:Arial,sans-serif;margin:0;color:#202124;background:#f3f5f6;line-height:1.45}}main{{max-width:1000px;margin:auto;padding:30px 24px 50px;background:white}}h1{{font-size:28px;margin:0 0 8px}}h2{{font-size:20px;margin-top:30px;border-bottom:2px solid #263238;padding-bottom:6px}}h3{{font-size:16px;margin:18px 0 4px}}.scope,.warning{{padding:12px 16px;margin:12px 0}}.scope{{border-left:5px solid #00897b;background:#eef8f6}}.warning{{border-left:5px solid #c62828;background:#fff4f3}}table{{width:100%;border-collapse:collapse;margin:12px 0;font-size:14px}}th,td{{border:1px solid #c9ced3;padding:8px;text-align:left;vertical-align:top}}th{{background:#eef1f3}}code{{font-family:Consolas,monospace}}img{{display:block;max-width:100%;height:auto;margin:14px auto;border:1px solid #c9ced3}}footer{{margin-top:34px;color:#5f6368;font-size:13px}}@media(max-width:650px){{main{{padding:18px 10px 36px}}table{{font-size:12px}}th,td{{padding:5px}}}}
</style></head><body><main>
<h1>Bounded validation report</h1><p><strong>System:</strong> {escape(manifest['system_id'])} &nbsp; <strong>Status:</strong> {escape(results['overall_status'])}</p>
<div class="scope"><strong>Validated scope.</strong> {escape(manifest['scientific_scope'])}</div>
<div class="warning"><strong>Non-claim.</strong> This report does not establish transferability to a new host, prototype, or bonding class and does not predict synthesizability.</div>
<h2>Claim boundary</h2><table>{claims}</table>
<h2>Visual evidence</h2><img src="validation_evidence.png" alt="Manifest-selected bounded validation evidence">
<h2>Module results</h2>{''.join(module_sections)}
<h2>Module contracts</h2>{contracts}<footer>Generated by experimental_extensions/run_validation.py from manifest-selected supplied inputs.</footer>
</main></body></html>"""
    output_path.write_text(html, encoding="utf-8")


def _manifold_pdf_lines(module_id: str, result: dict) -> list[str]:
    group = result.get("group_holdout")
    lines = [
        f"Module: {module_id}",
        f"Composition rows: {result['n_compositions']}",
        f"Mixing elements: {', '.join(result['mixing_elements'])}",
        f"Pair parameters: {result['n_pair_parameters']}",
        f"Training R2: {result['training_R2']:.9f}",
        f"Training RMSE: {result['training_metrics']['RMSE']:.6f} {result.get('energy_unit', '')}",
        f"Endmember-only RMSE: {result['endmember_only_nonendmember_metrics']['RMSE']:.6f}",
        f"Non-endmember LOOCV RMSE: {result['nonendmember_LOOCV_metrics']['RMSE']:.6f}",
    ]
    if group:
        lines.append(
            f"Group-holdout RMSE range: {group['RMSE_min']:.6f}-{group['RMSE_max']:.6f}"
        )
    lines.extend(["", "Interpretation: internal interpolation and ablation only."])
    return lines


def _comparison_pdf_lines(module_id: str, result: dict) -> list[str]:
    lines = [
        f"Module: {module_id}",
        "Backends: " + ", ".join(
            backend.get("label", backend["id"]) for backend in result["backends"]
        ),
        f"Unit: {result['unit']}",
        "",
    ]
    lines.extend(
        f"{row['subset_id']} {row['pair_id']}: N={row['n']}, rho={row['spearman_rho']:.4f}, RMSE={row['RMSE']:.3f}, top-{row['top_k']} overlap={row['top_k_overlap_count']}, reversals={row['ranking_reversals']}/{row['comparable_item_pairs']}"
        for row in result["pairwise_results"]
    )
    lines.extend(["", f"Limitation: {result['limitation']}"])
    return lines


def write_pdf(manifest: dict, results: dict, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    pages: list[tuple[str, list[str]]] = [
        (
            "Bounded validation report",
            [
                f"System: {manifest['system_id']}",
                f"Status: {results['overall_status']}", "",
                f"Validated scope: {manifest['scientific_scope']}", "",
                *[f"{key}: {value}" for key, value in manifest["claim_boundary"].items()],
                "", "This report does not predict synthesizability.",
            ],
        )
    ]
    for module_id, module in results["modules"].items():
        if module["kind"] == "manifold_regression":
            pages.append(
                ("Fixed-manifold validation", _manifold_pdf_lines(module_id, module["result"]))
            )
        elif module["kind"] == "energy_backend_comparison":
            pages.append(
                ("Energy-backend reference comparison", _comparison_pdf_lines(module_id, module["result"]))
            )
    with PdfPages(output_path) as pdf:
        for title, lines in pages:
            fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
            fig.text(0.08, 0.94, title, fontsize=18, weight="bold", va="top")
            y = 0.89
            for line in lines:
                pieces = wrap(line, width=98) if line else [""]
                for piece in pieces:
                    fig.text(0.08, y, piece, fontsize=10, va="top")
                    y -= 0.027
                if not line:
                    y -= 0.008
            fig.text(0.08, 0.035, "Bounded evidence artifact", fontsize=8, color="#666666")
            plt.axis("off")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
