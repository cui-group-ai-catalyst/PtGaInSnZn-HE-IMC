"""Configuration-driven validation of an integer composition manifold.

This is an experimental software-extension utility, not a manuscript
prediction engine. It validates the tabular composition manifold and fits a
pairwise CEF representation to supplied energies. Scientific validity of the
input structures and energies remains the user's responsibility.
"""
from __future__ import annotations

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


def _metrics(actual: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    error = predicted - actual
    return {
        "RMSE": float(np.sqrt(np.mean(error**2))),
        "MAE": float(np.mean(np.abs(error))),
        "max_abs_error": float(np.max(np.abs(error))),
    }


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def load_config(config_path: Path) -> dict:
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("schema_version") != 1:
        raise ValueError("Only schema_version=1 is supported")
    return config


def prepare(config: dict, config_path: Path):
    elements = tuple(config["mixing_elements"])
    site_count = int(config["mixing_site_count"])
    count_columns = config["count_columns"]
    if set(count_columns) != set(elements):
        raise ValueError("count_columns must map every mixing element exactly once")

    input_csv = _resolve(config_path, config["input_csv"])
    data = pd.read_csv(input_csv)
    required = [count_columns[element] for element in elements]
    required.append(config["energy_column"])
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing input columns: {missing}")

    count_frame = data[[count_columns[element] for element in elements]].copy()
    numeric = count_frame.to_numpy(dtype=float)
    if np.any(numeric < 0) or not np.allclose(numeric, np.round(numeric)):
        raise ValueError("Composition counts must be non-negative integers")
    counts = numeric.astype(int)
    if not np.all(counts.sum(axis=1) == site_count):
        raise ValueError(f"Every composition must sum to {site_count} mixing sites")
    if len({tuple(row) for row in counts}) != len(counts):
        raise ValueError("Duplicate integer composition vectors found")

    expected_rows = math.comb(site_count + len(elements) - 1, len(elements) - 1)
    if config.get("require_complete_integer_manifold", False) and len(data) != expected_rows:
        raise ValueError(
            f"Expected complete integer manifold with {expected_rows} rows, got {len(data)}"
        )

    fractions = counts / float(site_count)
    endmember = {}
    for index, element in enumerate(elements):
        mask = counts[:, index] == site_count
        if int(mask.sum()) != 1:
            raise ValueError(f"Expected exactly one {element} endmember")
        endmember[element] = float(data.loc[mask, config["energy_column"]].iloc[0])

    baseline = sum(
        fractions[:, index] * endmember[element]
        for index, element in enumerate(elements)
    )
    pairs = tuple(combinations(range(len(elements)), 2))
    design = np.column_stack(
        [fractions[:, left] * fractions[:, right] for left, right in pairs]
    )
    target = data[config["energy_column"]].to_numpy(dtype=float)
    is_endmember = counts.max(axis=1) == site_count
    return data, elements, counts, fractions, pairs, endmember, baseline, design, target, is_endmember


def fit(design: np.ndarray, residual: np.ndarray, mask: np.ndarray | None = None):
    if mask is None:
        mask = np.ones(len(residual), dtype=bool)
    omega, _, rank, singular = np.linalg.lstsq(design[mask], residual[mask], rcond=None)
    return omega, int(rank), singular


def run(config_path: Path, output_dir_override: Path | None = None) -> dict:
    config_path = config_path.resolve()
    config = load_config(config_path)
    (
        data,
        elements,
        counts,
        fractions,
        pairs,
        endmember,
        baseline,
        design,
        target,
        is_endmember,
    ) = prepare(config, config_path)

    residual = target - baseline
    omega, rank, singular = fit(design, residual)
    fitted = baseline + design @ omega
    train_metrics = _metrics(target, fitted)
    ss_res = float(np.sum((target - fitted) ** 2))
    ss_tot = float(np.sum((target - target.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot

    nonendmember = ~is_endmember
    baseline_metrics = _metrics(target[nonendmember], baseline[nonendmember])
    loocv = np.full(len(data), np.nan)
    for index in np.flatnonzero(nonendmember):
        train_mask = np.ones(len(data), dtype=bool)
        train_mask[index] = False
        fold_omega, _, _ = fit(design, residual, train_mask)
        loocv[index] = baseline[index] + design[index] @ fold_omega
    loocv_metrics = _metrics(target[nonendmember], loocv[nonendmember])

    group_rows = []
    group_element = config.get("group_holdout_element")
    if group_element:
        group_index = elements.index(group_element)
        for count in sorted(set(counts[:, group_index])):
            test_mask = (counts[:, group_index] == count) & nonendmember
            train_mask = counts[:, group_index] != count
            if not np.any(test_mask):
                continue
            fold_omega, fold_rank, _ = fit(design, residual, train_mask)
            predicted = baseline[test_mask] + design[test_mask] @ fold_omega
            group_rows.append(
                {
                    "held_out_element": group_element,
                    "held_out_count": int(count),
                    "n_train": int(train_mask.sum()),
                    "n_test_nonendmember": int(test_mask.sum()),
                    "design_rank": fold_rank,
                    **_metrics(target[test_mask], predicted),
                }
            )

    expected_rows = math.comb(
        int(config["mixing_site_count"]) + len(elements) - 1,
        len(elements) - 1,
    )
    summary = {
        "system_id": config["system_id"],
        "status": config["status"],
        "host": config.get("host"),
        "prototype": config.get("prototype"),
        "energy_unit": config["energy_unit"],
        "mixing_elements": list(elements),
        "mixing_site_count": int(config["mixing_site_count"]),
        "n_compositions": int(len(data)),
        "expected_complete_manifold_rows": int(expected_rows),
        "n_pair_parameters": int(len(pairs)),
        "design_rank": rank,
        "design_condition_number": float(singular[0] / singular[-1]),
        "training_R2": r2,
        "training_metrics": train_metrics,
        "endmember_only_nonendmember_metrics": baseline_metrics,
        "nonendmember_LOOCV_metrics": loocv_metrics,
        "group_holdout": (
            {
                "element": group_element,
                "n_groups": len(group_rows),
                "RMSE_min": float(min(row["RMSE"] for row in group_rows)),
                "RMSE_max": float(max(row["RMSE"] for row in group_rows)),
                "RMSE_mean": float(np.mean([row["RMSE"] for row in group_rows])),
            }
            if group_rows
            else None
        ),
        "scientific_scope": config["scientific_scope"],
        "transferability_claim": "none",
    }

    output_dir = (
        output_dir_override.resolve()
        if output_dir_override is not None
        else _resolve(config_path, config["output_dir"])
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    pd.DataFrame(
        [
            {
                "pair": f"{elements[left]}-{elements[right]}",
                "omega_per_input_energy_basis": float(omega[index]),
            }
            for index, (left, right) in enumerate(pairs)
        ]
    ).to_csv(output_dir / "pair_parameters.csv", index=False)
    if group_rows:
        pd.DataFrame(group_rows).to_csv(output_dir / "group_holdout.csv", index=False)

    print(json.dumps(summary, indent=2, ensure_ascii=True))
    print(f"Outputs: {output_dir}")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    run(args.config, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
