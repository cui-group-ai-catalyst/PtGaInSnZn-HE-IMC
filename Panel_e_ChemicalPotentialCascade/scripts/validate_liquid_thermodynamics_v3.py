"""Validate the corrected multicomponent liquid regular-solution implementation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import script_A_mu_liquid_v3 as liquid


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs"
REPORT = OUTPUT_DIR / "validation_liquid_thermodynamics_v3.json"
ATOL = 1.0e-10


def main() -> int:
    elements = list(liquid.LIQUID_COMPOSITION)
    x = np.array([liquid.LIQUID_COMPOSITION[element] for element in elements])
    matrix, _ = liquid.compute_omega_matrix(elements)

    if not np.isclose(x.sum(), 1.0, atol=ATOL):
        raise AssertionError("Liquid composition does not sum to one")
    if not np.allclose(matrix, matrix.T, atol=ATOL):
        raise AssertionError("Omega matrix is not symmetric")
    if not np.allclose(np.diag(matrix), 0.0, atol=ATOL):
        raise AssertionError("Omega matrix diagonal is not zero")

    ga_in_expected = (
        -liquid.P_MIEDEMA * (4.10 - 3.90) ** 2
        + liquid.Q_MIEDEMA * (1.34 - 1.17) ** 2
    )
    ga_index = elements.index("Ga")
    in_index = elements.index("In")
    if not np.isclose(matrix[ga_index, in_index], ga_in_expected, atol=ATOL):
        raise AssertionError("Ga-In Omega does not use the tabulated n_WS^(1/3) values")

    g_excess = sum(
        matrix[i, j] * x[i] * x[j]
        for i in range(len(elements))
        for j in range(i + 1, len(elements))
    )
    partial_excess = np.array(
        [
            sum(matrix[i, j] * x[j] for j in range(len(elements)) if j != i)
            - g_excess
            for i in range(len(elements))
        ]
    )
    gibbs_duhem_residual = float(np.dot(x, partial_excess) - g_excess)
    if abs(gibbs_duhem_residual) > ATOL:
        raise AssertionError("Partial molar excess terms violate the Euler identity")

    calculated = liquid.mu_liquid_SER(x, elements, matrix, T=0.0)
    expected_rows = {row["element"]: row for row in calculated}
    csv_rows = pd.read_csv(OUTPUT_DIR / "mu_liquid_v3_0K.csv").set_index("element")
    max_csv_difference = max(
        abs(
            float(csv_rows.loc[element, "mu_total_SER_kJmol"])
            - float(expected_rows[element]["mu_total_SER_kJmol"])
        )
        for element in expected_rows
    )
    if max_csv_difference > 5.0e-4:
        raise AssertionError("Liquid CSV is inconsistent with the corrected implementation")

    hei = pd.read_csv(OUTPUT_DIR / "mu_HEI_v3_0K.csv").set_index("element")
    reaction = pd.read_csv(OUTPUT_DIR / "delta_G_rxn_v3_summary.csv").iloc[0]
    g_atom_hei = float(hei.loc["Pt", "mu_HEI_per_atom_kJmol"])
    sum_y_mu_liquid = sum(
        liquid.LIQUID_COMPOSITION[element]
        * float(csv_rows.loc[element, "mu_total_SER_kJmol"])
        for element in elements
    )
    reconstructed_delta_g = 4.0 * g_atom_hei - sum_y_mu_liquid
    reported_delta_g = float(reaction["delta_G_rxn_per_fu_kJmol"])
    delta_g_residual = reconstructed_delta_g - reported_delta_g
    if abs(delta_g_residual) > 1.0e-3:
        raise AssertionError("Downstream reaction energy is inconsistent with A and B outputs")

    report = {
        "status": "passed",
        "liquid_model": "symmetric multicomponent regular solution",
        "n_ws_convention": "stored values are n_WS^(1/3)",
        "g_excess_formula": "sum_(i<j) Omega_ij*x_i*x_j",
        "partial_molar_excess_formula": "sum_(j!=i) Omega_ij*x_j - G_excess",
        "ga_in_omega_kJmol": float(matrix[ga_index, in_index]),
        "g_excess_kJmol": float(g_excess),
        "gibbs_duhem_euler_residual_kJmol": gibbs_duhem_residual,
        "max_csv_difference_kJmol": float(max_csv_difference),
        "reconstructed_delta_G_rxn_per_fu_kJmol": float(reconstructed_delta_g),
        "reported_delta_G_rxn_per_fu_kJmol": reported_delta_g,
        "delta_G_identity_residual_kJmol": float(delta_g_residual),
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
