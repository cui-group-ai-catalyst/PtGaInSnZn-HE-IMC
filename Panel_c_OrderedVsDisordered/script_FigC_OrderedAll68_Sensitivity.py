"""Check whether the 30-class ordered sample represents all 68 symmetry classes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DEFAULT_SAMPLE = HERE / "data_FigC_OrderedEnsemble_Raw_UMA_Energies_regen.csv"
DEFAULT_ALL = HERE / "data_FigC_OrderedAll68_Raw_UMA_Energies_sensitivity.csv"
DEFAULT_REPORT = HERE / "validation_FigC_OrderedAll68_Sensitivity.json"
EV_TO_KJMOL = 96.485


def stats(values: np.ndarray) -> dict[str, float | int]:
    return {
        "n": len(values),
        "mean_eV_atom": float(np.mean(values)),
        "sample_sd_eV_atom": float(np.std(values, ddof=1)),
        "range_eV_atom": float(np.max(values) - np.min(values)),
        "sample_sd_kJ_mol_atom": float(np.std(values, ddof=1) * EV_TO_KJMOL),
        "range_kJ_mol_atom": float((np.max(values) - np.min(values)) * EV_TO_KJMOL),
        "min_eV_atom": float(np.min(values)),
        "max_eV_atom": float(np.max(values)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--all-classes", type=Path, default=DEFAULT_ALL)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    sample = pd.read_csv(args.sample)
    all_classes = pd.read_csv(args.all_classes)
    if len(sample) != 30 or sample["Canonical_Class"].nunique() != 30:
        raise ValueError("Expected 30 distinct sampled symmetry classes")
    if len(all_classes) != 68 or all_classes["Canonical_Class"].nunique() != 68:
        raise ValueError("Expected all 68 distinct symmetry classes")
    if int(all_classes["Class_Degeneracy"].sum()) != 2520:
        raise ValueError("Class degeneracies do not recover all 2520 labelled assignments")

    sample_values = sample["Energy_eV_atom"].to_numpy(float)
    all_values = all_classes["Energy_eV_atom"].to_numpy(float)
    weights = all_classes["Class_Degeneracy"].to_numpy(float)
    weighted_mean = float(np.average(all_values, weights=weights))
    weighted_variance = float(np.average((all_values - weighted_mean) ** 2, weights=weights))
    sample_mean = float(np.mean(sample_values))
    mean_difference = (weighted_mean - sample_mean) * EV_TO_KJMOL
    anchor_difference = abs(float(sample_values[0]) - float(all_values[0]))

    report = {
        "status": "passed",
        "provenance_mode": "measured_computational",
        "energy_protocol": "fixed-cell fixed-coordinate UMA-s-1p1 single-point",
        "sample_30": stats(sample_values),
        "all_68_equal_class_weight": stats(all_values),
        "all_68_degeneracy_weighted": {
            "raw_assignment_count": int(weights.sum()),
            "mean_eV_atom": weighted_mean,
            "population_sd_eV_atom": math.sqrt(weighted_variance),
            "population_sd_kJ_mol_atom": math.sqrt(weighted_variance) * EV_TO_KJMOL,
        },
        "weighted_all68_minus_sample30_mean_kJ_mol_atom": mean_difference,
        "historical_anchor_abs_difference_eV_atom": anchor_difference,
        "interpretation": (
            "The 30-class sample reproduces the degeneracy-weighted all-class mean; "
            "the narrow spread is not caused by duplicate or symmetry-equivalent structures."
        ),
    }
    if abs(mean_difference) > 0.01:
        report["status"] = "failed"
        raise AssertionError("The 30-class sample mean differs materially from the all-class mean")
    if anchor_difference > 1.0e-6:
        report["status"] = "failed"
        raise AssertionError("Historical anchor is inconsistent between the two runs")

    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
