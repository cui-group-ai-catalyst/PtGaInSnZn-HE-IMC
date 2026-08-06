"""Compare matched energy backends without making transferability claims."""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from experimental_extensions.contracts import resolve_relative
except ModuleNotFoundError:
    from contracts import resolve_relative


def _pair_id(left: dict, right: dict) -> str:
    return f"{left['id']}__{right['id']}"


def _metrics(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    delta = right - left
    left_rank = pd.Series(left).rank(method="average").to_numpy(dtype=float)
    right_rank = pd.Series(right).rank(method="average").to_numpy(dtype=float)
    return {
        "spearman_rho": float(np.corrcoef(left_rank, right_rank)[0, 1]),
        "pearson_r": float(np.corrcoef(left, right)[0, 1]),
        "RMSE": float(np.sqrt(np.mean(delta**2))),
        "MAE": float(np.mean(np.abs(delta))),
        "bias_right_minus_left": float(np.mean(delta)),
    }


def _subset_mask(data: pd.DataFrame, subset: dict) -> np.ndarray:
    rule = subset.get("filter")
    if rule is None:
        return np.ones(len(data), dtype=bool)
    column = rule["column"]
    if column not in data.columns:
        raise ValueError(f"Missing subset-filter column: {column}")
    return (data[column] == rule["equals"]).to_numpy(dtype=bool)


def _ranking_reversals(
    keys: list[str], left: np.ndarray, right: np.ndarray
) -> tuple[int, int, int, list[dict]]:
    reversals: list[dict] = []
    comparable = 0
    ties = 0
    for first, second in combinations(range(len(keys)), 2):
        left_delta = float(left[first] - left[second])
        right_delta = float(right[first] - right[second])
        if left_delta == 0.0 or right_delta == 0.0:
            ties += 1
            continue
        comparable += 1
        if left_delta * right_delta < 0.0:
            reversals.append(
                {
                    "item_a": keys[first],
                    "item_b": keys[second],
                    "left_preferred": keys[first] if left_delta > 0 else keys[second],
                    "right_preferred": keys[first] if right_delta > 0 else keys[second],
                    "left_difference": left_delta,
                    "right_difference": right_delta,
                }
            )
    return len(reversals), comparable, ties, reversals


def run_comparison(module: dict, manifest_path: Path, output_dir: Path) -> dict:
    config = module["config"]
    input_csv = resolve_relative(manifest_path, config["input_csv"])
    data = pd.read_csv(input_csv)
    key_column = config["key_column"]
    backends = config["backends"]
    required = [key_column, *[backend["column"] for backend in backends]]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing comparison columns: {missing}")
    if data[key_column].duplicated().any():
        raise ValueError(f"Comparison key column {key_column!r} contains duplicates")
    if config.get("ranking_direction") != "descending":
        raise ValueError("Schema v1 supports descending energy-magnitude ranks only")

    top_k = int(config["top_k"])
    pairwise_results: list[dict] = []
    reversal_rows: list[dict] = []
    top_k_rows: list[dict] = []
    ranking_table = data[[key_column]].copy()
    for backend in backends:
        ranking_table[f"{backend['id']}_value"] = data[backend["column"]].astype(float)
        ranking_table[f"{backend['id']}_rank"] = data[backend["column"]].rank(
            method="average", ascending=False
        )

    for subset in config["subsets"]:
        mask = _subset_mask(data, subset)
        subset_data = data.loc[mask].reset_index(drop=True)
        n_rows = len(subset_data)
        if n_rows != int(subset["expected_n"]):
            raise AssertionError(
                f"Subset {subset['id']!r}: observed {n_rows}, expected {subset['expected_n']}"
            )
        if n_rows < max(3, top_k):
            raise ValueError(f"Subset {subset['id']!r} is too small for top-{top_k}")
        keys = subset_data[key_column].astype(str).tolist()
        expected_rho = subset.get("expected_spearman_rho", {})

        for left_backend, right_backend in combinations(backends, 2):
            pair_id = _pair_id(left_backend, right_backend)
            left = subset_data[left_backend["column"]].to_numpy(dtype=float)
            right = subset_data[right_backend["column"]].to_numpy(dtype=float)
            metrics = _metrics(left, right)
            if pair_id in expected_rho:
                tolerance = float(subset["tolerance"])
                if abs(metrics["spearman_rho"] - float(expected_rho[pair_id])) > tolerance:
                    raise AssertionError(
                        f"{subset['id']} {pair_id}: Spearman rho "
                        f"{metrics['spearman_rho']:.6f} does not match "
                        f"{float(expected_rho[pair_id]):.6f} +/- {tolerance:.6f}"
                    )
            left_order = np.argsort(-left, kind="stable")[:top_k]
            right_order = np.argsort(-right, kind="stable")[:top_k]
            left_top = [keys[index] for index in left_order]
            right_top = [keys[index] for index in right_order]
            overlap = sorted(set(left_top) & set(right_top))
            union = set(left_top) | set(right_top)
            n_reversals, comparable, ties, reversals = _ranking_reversals(
                keys, left, right
            )
            result = {
                "subset_id": subset["id"],
                "pair_id": pair_id,
                "left_backend": left_backend["id"],
                "right_backend": right_backend["id"],
                "n": n_rows,
                **metrics,
                "top_k": top_k,
                "top_k_overlap_count": len(overlap),
                "top_k_jaccard": float(len(overlap) / len(union)),
                "top_k_overlap": overlap,
                "ranking_reversals": n_reversals,
                "comparable_item_pairs": comparable,
                "tied_item_pairs": ties,
                "ranking_reversal_fraction": float(n_reversals / comparable),
            }
            pairwise_results.append(result)
            for backend_id, entries in (
                (left_backend["id"], left_top),
                (right_backend["id"], right_top),
            ):
                for rank, key in enumerate(entries, start=1):
                    top_k_rows.append(
                        {
                            "subset_id": subset["id"],
                            "pair_id": pair_id,
                            "backend": backend_id,
                            "rank": rank,
                            "item": key,
                        }
                    )
            for reversal in reversals:
                reversal_rows.append(
                    {"subset_id": subset["id"], "pair_id": pair_id, **reversal}
                )

    summary = {
        "module_id": module["id"],
        "status": "passed",
        "backends": backends,
        "key_column": key_column,
        "unit": config["unit"],
        "ranking_direction": config["ranking_direction"],
        "primary_metric": config["primary_metric"],
        "pairwise_results": pairwise_results,
        "ranking_rows": ranking_table.to_dict(orient="records"),
        "interpretation": config["interpretation"],
        "limitation": config["limitation"],
        "transferability_claim": "none",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    flat_metrics = [
        {key: value for key, value in row.items() if not isinstance(value, list)}
        for row in pairwise_results
    ]
    pd.DataFrame(flat_metrics).to_csv(
        output_dir / "comparison_metrics.csv", index=False
    )
    pd.DataFrame(reversal_rows).to_csv(
        output_dir / "ranking_reversals.csv", index=False
    )
    pd.DataFrame(top_k_rows).to_csv(output_dir / "top_k_members.csv", index=False)
    ranking_table.to_csv(output_dir / "ranking_table.csv", index=False)
    return summary
