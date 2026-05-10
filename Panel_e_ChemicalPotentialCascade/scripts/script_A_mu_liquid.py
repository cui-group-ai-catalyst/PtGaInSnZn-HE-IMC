"""
script_A_mu_liquid.py
=====================
Compute liquid-phase partial molar chemical potentials mu_i^(0) for the
representative liquid cocktail (Ga 0.65 / In 0.20 / Sn 0.10 / Zn 0.05)
plus Pt pure solid reference.

Method
------
Miedema-derived binary interaction parameter Omega_ij (simplified form,
matching Panel d / b convention: no volume prefactor) applied within the
Hildebrand multicomponent regular-solution framework to obtain per-element
partial molar chemical potentials.

Convention
----------
- 0 K enthalpic approximation (matches Panel b gamma_SL Miedema and
  Panel g DFT single-point energies).
- Pure-element liquid reference state: mu_i^pure_L = 0 by construction.
- Pt at t=0 is pure solid substrate (not dissolved in liquid), so
  mu_Pt^(0) = 0 by pure-Pt reference convention.

Outputs
-------
- outputs/mu_liquid_0K.csv              Main Panel a numbers
- outputs/mu_liquid_T_table.csv         SI temperature sensitivity
- outputs/omega_liquid_binary.csv       Audit trail of Miedema Omegas

References
----------
- Miedema, de Chatel, de Boer 1980, Physica B 100, 1.
- de Boer et al. 1988, Cohesion in Metals, North-Holland.
- Hildebrand 1929, J. Am. Chem. Soc. 51, 66.
- Guggenheim 1935, Proc. Roy. Soc. A 148, 304.

Notes
-----
The simplified Miedema form -P(dPhi)^2 + Q(dn_WS^(1/3))^2 omits the
element-specific volume prefactor used in the full de Boer 1988 expression.
This is the same simplification used by Panel d (`script_FigD_compute_and_plot.py`);
consistency across Panels a / b / d / g is therefore preserved. A volume-
corrected variant can be toggled via VOLUME_CORRECTED = True for cross-check.

Author: 2026-04-23 (PtGaInSnZn submission prep)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# INPUTS (hard-coded to isolate Script A from upstream data dependencies;
#         values taken from 01_scripts/data_periodic_table.py and
#         04_Panel_d_Wetting/script_FigD_compute_and_plot.py)
# --------------------------------------------------------------------------
PERIODIC: dict[str, dict[str, float]] = {
    # Phi in V, n_ws in density units (Miedema 1980 convention), V in cm^3/mol
    "Ga": {"Phi": 4.10, "n_ws": 1.34, "V": 11.8},
    "In": {"Phi": 3.90, "n_ws": 1.17, "V": 15.7},
    "Sn": {"Phi": 4.15, "n_ws": 1.25, "V": 16.3},
    "Zn": {"Phi": 4.10, "n_ws": 1.32, "V": 9.17},
    "Pt": {"Phi": 5.65, "n_ws": 1.78, "V": 9.10},
}

# Miedema constants for transition / non-transition liquid metals
# (Panel d convention: no *10 scaling)
P_MIEDEMA: float = 14.1  # kJ/mol per (dPhi)^2
Q_MIEDEMA: float = 9.4   # kJ/mol per (dn_WS^(1/3))^2

# Liquid cocktail composition -- LOCKED to Panel b baseline
LIQUID_COMPOSITION: dict[str, float] = {
    "Ga": 0.65,
    "In": 0.20,
    "Sn": 0.10,
    "Zn": 0.05,
}

# Temperatures (K) for SI sensitivity table
SI_TEMPERATURES: list[int] = [300, 800, 1000, 1200]

# Gas constant in kJ/(mol*K)
R_KJ_MOL_K: float = 8.31446261815324e-3

# Optional: use volume-prefactor-corrected Miedema for cross-check
VOLUME_CORRECTED: bool = False


# --------------------------------------------------------------------------
# CORE FUNCTIONS
# --------------------------------------------------------------------------
def miedema_omega_simple(elem_i: str, elem_j: str) -> float:
    """Simplified Miedema binary interaction parameter.

    Omega_ij = -P * (Phi_i - Phi_j)^2 + Q * (n_ws_i^(1/3) - n_ws_j^(1/3))^2   [kJ/mol]

    This is the form used in Panel d (`script_FigD_compute_and_plot.py`).
    """
    pi = PERIODIC[elem_i]
    pj = PERIODIC[elem_j]
    d_phi = pi["Phi"] - pj["Phi"]
    d_n13 = pi["n_ws"] ** (1.0 / 3.0) - pj["n_ws"] ** (1.0 / 3.0)
    return -P_MIEDEMA * d_phi**2 + Q_MIEDEMA * d_n13**2


def miedema_omega_volume_corrected(elem_i: str, elem_j: str) -> float:
    """Volume-prefactor-corrected Miedema (de Boer 1988, simplified).

    Omega_ij = 2 * (V_i^(2/3) + V_j^(2/3)) / (n_ws_i^(-1/3) + n_ws_j^(-1/3))
              * [-P * (dPhi)^2 + Q * (dn_WS^(1/3))^2]   [kJ/mol]

    This converges on the full de Boer expression at equiatomic and is
    offered here only as a sanity-check alternative.
    """
    pi = PERIODIC[elem_i]
    pj = PERIODIC[elem_j]
    d_phi = pi["Phi"] - pj["Phi"]
    d_n13 = pi["n_ws"] ** (1.0 / 3.0) - pj["n_ws"] ** (1.0 / 3.0)
    vf = (pi["V"] ** (2.0 / 3.0) + pj["V"] ** (2.0 / 3.0))
    nf = (pi["n_ws"] ** (-1.0 / 3.0) + pj["n_ws"] ** (-1.0 / 3.0))
    prefactor = 2.0 * vf / nf
    return prefactor * (-P_MIEDEMA * d_phi**2 + Q_MIEDEMA * d_n13**2)


def compute_omega_matrix(elements: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    """Return symmetric matrix of Omega_ij and an audit DataFrame."""
    omega_fn = miedema_omega_volume_corrected if VOLUME_CORRECTED else miedema_omega_simple
    n = len(elements)
    M = np.zeros((n, n))
    rows: list[dict] = []
    for i, ei in enumerate(elements):
        for j, ej in enumerate(elements):
            if i == j:
                continue
            M[i, j] = omega_fn(ei, ej)
            if i < j:
                rows.append(
                    {
                        "i": ei,
                        "j": ej,
                        "dPhi_V": round(PERIODIC[ei]["Phi"] - PERIODIC[ej]["Phi"], 3),
                        "dn_WS_13": round(
                            PERIODIC[ei]["n_ws"] ** (1.0 / 3.0)
                            - PERIODIC[ej]["n_ws"] ** (1.0 / 3.0),
                            4,
                        ),
                        "Omega_ij_kJmol": round(M[i, j], 3),
                    }
                )
    return M, pd.DataFrame(rows)


def chemical_potential_excess(
    x: np.ndarray, elements: list[str], M: np.ndarray
) -> dict[str, float]:
    """Hildebrand multicomponent excess chemical potential (regular solution).

    mu_i^E = sum_{j != i} Omega_ij * x_j^2
    """
    mu_e: dict[str, float] = {}
    for i, ei in enumerate(elements):
        s = 0.0
        for j, _ej in enumerate(elements):
            if j == i:
                continue
            s += M[i, j] * x[j] ** 2
        mu_e[ei] = s
    return mu_e


def chemical_potential_total(
    x: np.ndarray, elements: list[str], M: np.ndarray, T: float
) -> tuple[dict[str, float], dict[str, float]]:
    """Total partial molar chemical potential.

    mu_i = mu_i^E + RT ln x_i        (pure-liquid reference, G_i^0 = 0)

    Returns (mu_total, breakdown) where breakdown has 'mu_E', 'RT_ln_x'.
    """
    mu_e = chemical_potential_excess(x, elements, M)
    mu_total: dict[str, float] = {}
    breakdown: dict[str, dict[str, float]] = {}
    for i, ei in enumerate(elements):
        if x[i] <= 0:
            mu_total[ei] = float("nan")
            breakdown[ei] = {"mu_E": float("nan"), "RT_ln_x": float("nan")}
            continue
        rt_ln = R_KJ_MOL_K * T * np.log(x[i]) if T > 0 else 0.0
        mu_total[ei] = mu_e[ei] + rt_ln
        breakdown[ei] = {"mu_E": mu_e[ei], "RT_ln_x": rt_ln}
    return mu_total, breakdown


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main() -> None:
    liquid_elements = list(LIQUID_COMPOSITION.keys())
    x_liquid = np.array([LIQUID_COMPOSITION[e] for e in liquid_elements], dtype=float)
    assert abs(x_liquid.sum() - 1.0) < 1e-9, "Composition must sum to 1"

    # Step 1: Omega matrix for 4 liquid elements
    M_liq, omega_df = compute_omega_matrix(liquid_elements)

    # Step 2: chemical potentials @ 0 K (main Panel a output)
    mu_0K, bd_0K = chemical_potential_total(x_liquid, liquid_elements, M_liq, T=0.0)
    main_rows: list[dict] = []
    for e in liquid_elements:
        main_rows.append(
            {
                "element": e,
                "phase": "liquid_cocktail",
                "x_i_or_y_i": LIQUID_COMPOSITION[e],
                "mu_E_kJmol": round(bd_0K[e]["mu_E"], 3),
                "mu_ideal_kJmol": round(bd_0K[e]["RT_ln_x"], 3),
                "mu_total_kJmol": round(mu_0K[e], 3),
                "note": "0K enthalpic (excess only)",
            }
        )
    main_rows.append(
        {
            "element": "Pt",
            "phase": "solid_fcc_pure",
            "x_i_or_y_i": 1.0,
            "mu_E_kJmol": 0.0,
            "mu_ideal_kJmol": 0.0,
            "mu_total_kJmol": 0.0,
            "note": "Pure Pt substrate reference at t=0",
        }
    )
    df_main = pd.DataFrame(main_rows)

    # Step 3: SI temperature-sensitivity table
    si_rows: list[dict] = []
    for T in SI_TEMPERATURES:
        mu_T, bd_T = chemical_potential_total(x_liquid, liquid_elements, M_liq, T=float(T))
        for e in liquid_elements:
            si_rows.append(
                {
                    "T_K": T,
                    "element": e,
                    "phase": "liquid_cocktail",
                    "x_i": LIQUID_COMPOSITION[e],
                    "mu_E_kJmol": round(bd_T[e]["mu_E"], 3),
                    "RT_ln_x_kJmol": round(bd_T[e]["RT_ln_x"], 3),
                    "mu_total_kJmol": round(mu_T[e], 3),
                }
            )
        si_rows.append(
            {
                "T_K": T,
                "element": "Pt",
                "phase": "solid_fcc_pure",
                "x_i": 1.0,
                "mu_E_kJmol": 0.0,
                "RT_ln_x_kJmol": 0.0,
                "mu_total_kJmol": 0.0,
            }
        )
    df_T = pd.DataFrame(si_rows)

    # Step 4: write outputs
    out_dir = Path(__file__).resolve().parent.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out_dir / "mu_liquid_0K.csv", index=False, encoding="utf-8")
    df_T.to_csv(out_dir / "mu_liquid_T_table.csv", index=False, encoding="utf-8")
    omega_df.to_csv(out_dir / "omega_liquid_binary.csv", index=False, encoding="utf-8")

    meta = {
        "script": "script_A_mu_liquid.py",
        "composition": LIQUID_COMPOSITION,
        "miedema": {"P": P_MIEDEMA, "Q": Q_MIEDEMA, "volume_corrected": VOLUME_CORRECTED},
        "temperatures_K": SI_TEMPERATURES,
        "reference_state": "pure_liquid_for_cocktail_elements; pure_solid_fcc_for_Pt",
        "convention": "0K enthalpic main output; RT ln x included only in T table",
    }
    (out_dir / "script_A_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # Step 5: print summary
    print("=" * 72)
    print("Script A: Liquid-phase chemical potentials mu_i^(0)")
    print("=" * 72)
    print(f"Composition      : {LIQUID_COMPOSITION}")
    print(f"Miedema P, Q     : {P_MIEDEMA}, {Q_MIEDEMA}")
    print(f"Volume-corrected : {VOLUME_CORRECTED}")
    print()
    print("-- Binary Omega_ij^L (kJ/mol) --")
    print(omega_df.to_string(index=False))
    print()
    print("-- Main output: mu_i^(0) @ 0 K enthalpic (for Panel a) --")
    print(
        df_main[
            ["element", "phase", "x_i_or_y_i", "mu_E_kJmol", "mu_total_kJmol", "note"]
        ].to_string(index=False)
    )
    print()
    print("-- SI temperature sensitivity (sample head) --")
    print(df_T.head(10).to_string(index=False))
    print()
    print(f"Outputs written to: {out_dir}")
    print("Files:")
    for f in sorted(out_dir.glob("*")):
        print(f"  - {f.name}")


if __name__ == "__main__":
    sys.exit(main() or 0)
