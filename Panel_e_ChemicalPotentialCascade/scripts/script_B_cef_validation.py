"""Validate the six-parameter CEF used for HEI chemical potentials.

The tests cover interpolation inside the fixed L1_2-Pt3(Ga,In,Sn,Zn)
composition manifold. They do not test UMA accuracy or transferability to new
elements, hosts, structure prototypes, or bonding chemistries.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ELEMENTS = ("Ga", "In", "Sn", "Zn")
PAIRS = tuple(
    (left, right)
    for index, left in enumerate(ELEMENTS)
    for right in ELEMENTS[index + 1 :]
)
SCRIPT_DIR = Path(__file__).resolve().parent
RELEASE_ROOT = SCRIPT_DIR.parents[1]
INPUT_CSV = (
    RELEASE_ROOT
    / "SI_Figures"
    / "SI_Fig04_165CompositionLandscape"
    / "data_FigG_165_ElementReferenced_Hf.csv"
)
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs"


def prepare_design(df: pd.DataFrame):
    data = df.copy()
    for element in ELEMENTS:
        data[f"y_{element}"] = data[f"{element}_count"] / 8.0

    endmember = {
        element: float(
            data.loc[
                data[f"{element}_count"] == 8,
                "ElementRef_Hf_kJ_mol",
            ].iloc[0]
        )
        for element in ELEMENTS
    }
    baseline = sum(
        data[f"y_{element}"] * endmember[element] for element in ELEMENTS
    ).to_numpy()
    design = np.column_stack(
        [
            (data[f"y_{left}"] * data[f"y_{right}"]).to_numpy()
            for left, right in PAIRS
        ]
    )
    target = data["ElementRef_Hf_kJ_mol"].to_numpy()
    is_endmember = np.array(
        [
            max(int(row[f"{element}_count"]) for element in ELEMENTS) == 8
            for _, row in data.iterrows()
        ]
    )
    return data, baseline, design, target, is_endmember


def fit(design: np.ndarray, residual: np.ndarray, mask=None):
    if mask is None:
        mask = np.ones(len(residual), dtype=bool)
    omega, _, rank, singular = np.linalg.lstsq(
        design[mask], residual[mask], rcond=None
    )
    return omega, int(rank), singular


def metrics(actual: np.ndarray, predicted: np.ndarray):
    error = predicted - actual
    return {
        "RMSE": float(np.sqrt(np.mean(error**2))),
        "MAE": float(np.mean(np.abs(error))),
        "max_abs_error": float(np.max(np.abs(error))),
    }


def main():
    df = pd.read_csv(INPUT_CSV)
    if len(df) != 165:
        raise ValueError(f"Expected 165 rows, got {len(df)}")

    data, baseline, design, target, is_endmember = prepare_design(df)
    residual = target - baseline
    nonendmember = ~is_endmember

    omega, rank, singular = fit(design, residual)
    fitted = baseline + design @ omega
    train = metrics(target, fitted)
    ss_res = float(np.sum((target - fitted) ** 2))
    ss_tot = float(np.sum((target - target.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot
    baseline_metrics = metrics(
        target[nonendmember], baseline[nonendmember]
    )

    loocv = np.full(len(data), np.nan)
    for index in np.flatnonzero(nonendmember):
        train_mask = np.ones(len(data), dtype=bool)
        train_mask[index] = False
        fold_omega, _, _ = fit(design, residual, train_mask)
        loocv[index] = baseline[index] + design[index] @ fold_omega
    loocv_metrics = metrics(
        target[nonendmember], loocv[nonendmember]
    )

    group_rows = []
    for ga_count in sorted(data["Ga_count"].unique()):
        test_mask = (data["Ga_count"].to_numpy() == ga_count) & nonendmember
        train_mask = data["Ga_count"].to_numpy() != ga_count
        if not np.any(test_mask):
            continue
        fold_omega, fold_rank, _ = fit(design, residual, train_mask)
        predicted = baseline[test_mask] + design[test_mask] @ fold_omega
        fold_metrics = metrics(target[test_mask], predicted)
        group_rows.append(
            {
                "held_out_Ga_count": int(ga_count),
                "n_train": int(np.sum(train_mask)),
                "n_test_nonendmember": int(np.sum(test_mask)),
                "design_rank": fold_rank,
                **fold_metrics,
            }
        )

    summary_rows = [
        ("n_compositions", len(data), "count", "full dataset"),
        ("n_pair_parameters", len(PAIRS), "count", "pairwise CEF"),
        ("design_rank", rank, "count", "pairwise CEF"),
        (
            "design_condition_number",
            float(singular[0] / singular[-1]),
            "dimensionless",
            "pairwise CEF",
        ),
        ("training_R2", r2, "dimensionless", "all 165 compositions"),
        (
            "training_RMSE",
            train["RMSE"],
            "kJ mol^-1 atom^-1",
            "all 165 compositions",
        ),
        (
            "training_max_abs_error",
            train["max_abs_error"],
            "kJ mol^-1 atom^-1",
            "all 165 compositions",
        ),
        (
            "endmember_only_RMSE",
            baseline_metrics["RMSE"],
            "kJ mol^-1 atom^-1",
            "161 non-endmember compositions",
        ),
        (
            "nonendmember_LOOCV_RMSE",
            loocv_metrics["RMSE"],
            "kJ mol^-1 atom^-1",
            "161 non-endmember compositions",
        ),
        (
            "nonendmember_LOOCV_MAE",
            loocv_metrics["MAE"],
            "kJ mol^-1 atom^-1",
            "161 non-endmember compositions",
        ),
        (
            "nonendmember_LOOCV_max_abs_error",
            loocv_metrics["max_abs_error"],
            "kJ mol^-1 atom^-1",
            "161 non-endmember compositions",
        ),
    ]
    summary = pd.DataFrame(
        summary_rows, columns=["metric", "value", "unit", "scope"]
    )

    prediction_columns = [
        "Composition",
        "Ga_count",
        "In_count",
        "Sn_count",
        "Zn_count",
    ]
    predictions = data[prediction_columns].copy()
    predictions["is_endmember"] = is_endmember
    predictions["actual_kJmol_atom"] = target
    predictions["endmember_only_prediction_kJmol_atom"] = baseline
    predictions["cef_fitted_prediction_kJmol_atom"] = fitted
    predictions["cef_fitted_error_kJmol_atom"] = fitted - target
    predictions["cef_loocv_prediction_kJmol_atom"] = loocv
    predictions["cef_loocv_error_kJmol_atom"] = loocv - target

    parameter_rows = [
        {
            "pair": f"{left}-{right}",
            "omega_kJmol_atom": float(omega[index]),
            "Omega_kJmol_beta_site": 4.0 * float(omega[index]),
            "unit_relation": "omega_atom = Omega_beta_site / 4",
        }
        for index, (left, right) in enumerate(PAIRS)
    ]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_DIR / "cef_validation_summary.csv", index=False)
    predictions.to_csv(
        OUTPUT_DIR / "cef_validation_predictions.csv", index=False
    )
    pd.DataFrame(group_rows).to_csv(
        OUTPUT_DIR / "cef_validation_ga_group_cv.csv", index=False
    )
    pd.DataFrame(parameter_rows).to_csv(
        OUTPUT_DIR / "cef_validation_parameters.csv", index=False
    )

    print(summary.to_string(index=False))
    print("\nGa-count group holdout validation")
    print(pd.DataFrame(group_rows).to_string(index=False))


if __name__ == "__main__":
    main()
