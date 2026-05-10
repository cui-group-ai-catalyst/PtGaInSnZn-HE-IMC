"""
script_A_mu_liquid_v2.py
========================
Compute liquid-phase chemical potentials for the representative cocktail
(Ga 0.65 / In 0.20 / Sn 0.10 / Zn 0.05) and Pt solid substrate,
referenced to **SER (Stable Element Reference)** = 0 K stable solid phase.

v2 change vs. v1
-----------------
v1 used pure-liquid reference (mu_i^pure_L = 0), making all liquid mu ~ 0.
v2 uses SER: mu_i^L = +DH_fusion(i) + Miedema excess, so liquid metals sit
ABOVE the solid reference, correctly reflecting their higher energy state.

Physical meaning
----------------
Liquid metals are metastable relative to their ground-state solids at 0 K.
The latent heat of fusion DH_fus represents the enthalpic cost of the liquid
state. In the SER convention this appears as a POSITIVE offset, placing the
liquid cocktail above the Pt solid substrate on the chemical-potential axis.

Data sources
------------
- DH_fusion: CRC Handbook of Chemistry and Physics, 97th ed.
- Miedema Omega_ij: de Boer et al. 1988, Cohesion in Metals (simplified,
  no V-prefactor, matching Panel d convention).

Outputs (all in ../outputs/)
-----------------------------
- mu_liquid_v2_0K.csv          Main Panel a starting points (SER ref)
- mu_liquid_v2_T_table.csv     SI temperature sensitivity
- omega_liquid_binary_v2.csv   Omega audit (unchanged from v1)
- script_A_v2_meta.json        Provenance

Author: 2026-04-23 v2
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ?? INPUTS ??????????????????????????????????????????????????????????????
PERIODIC: dict[str, dict[str, float]] = {
    "Ga": {"Phi": 4.10, "n_ws": 1.34, "V": 11.8},
    "In": {"Phi": 3.90, "n_ws": 1.17, "V": 15.7},
    "Sn": {"Phi": 4.15, "n_ws": 1.25, "V": 16.3},
    "Zn": {"Phi": 4.10, "n_ws": 1.32, "V": 9.17},
    "Pt": {"Phi": 5.65, "n_ws": 1.78, "V": 9.10},
}

# Latent heat of fusion (kJ/mol), CRC Handbook 97th ed.
DH_FUSION: dict[str, float] = {
    "Ga": 5.59,
    "In": 3.28,
    "Sn": 7.03,
    "Zn": 7.32,
    "Pt": 0.0,   # Pt stays solid; SER reference = 0
}

P_MIEDEMA: float = 14.1
Q_MIEDEMA: float = 9.4
R_KJ_MOL_K: float = 8.31446261815324e-3

LIQUID_COMPOSITION: dict[str, float] = {
    "Ga": 0.65, "In": 0.20, "Sn": 0.10, "Zn": 0.05,
}

SI_TEMPERATURES: list[int] = [300, 800, 1000, 1200]

BETA_ELEMENTS: list[str] = ["Ga", "In", "Sn", "Zn"]

# ?? CORE FUNCTIONS ??????????????????????????????????????????????????????
def miedema_omega(ei: str, ej: str) -> float:
    pi, pj = PERIODIC[ei], PERIODIC[ej]
    d_phi = pi["Phi"] - pj["Phi"]
    d_n13 = pi["n_ws"] ** (1.0 / 3.0) - pj["n_ws"] ** (1.0 / 3.0)
    return -P_MIEDEMA * d_phi**2 + Q_MIEDEMA * d_n13**2


def compute_omega_matrix(elements: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    n = len(elements)
    M = np.zeros((n, n))
    rows: list[dict] = []
    for i, ei in enumerate(elements):
        for j, ej in enumerate(elements):
            if i == j:
                continue
            M[i, j] = miedema_omega(ei, ej)
            if i < j:
                rows.append({
                    "i": ei, "j": ej,
                    "dPhi_V": round(PERIODIC[ei]["Phi"] - PERIODIC[ej]["Phi"], 3),
                    "dn_WS_13": round(
                        PERIODIC[ei]["n_ws"] ** (1/3) - PERIODIC[ej]["n_ws"] ** (1/3), 4),
                    "Omega_ij_kJmol": round(M[i, j], 3),
                })
    return M, pd.DataFrame(rows)


def mu_liquid_SER(
    x: np.ndarray, elements: list[str], M: np.ndarray, T: float,
) -> list[dict]:
    """Chemical potentials in SER reference.

    mu_i^L(SER) = DH_fus(i) + RT ln x_i + sum_j Omega_ij x_j^2
    """
    rows: list[dict] = []
    for i, ei in enumerate(elements):
        dh_fus = DH_FUSION[ei]
        excess = sum(M[i, j] * x[j] ** 2 for j in range(len(elements)) if j != i)
        rt_ln = R_KJ_MOL_K * T * np.log(x[i]) if (T > 0 and x[i] > 0) else 0.0
        total = dh_fus + excess + rt_ln
        rows.append({
            "element": ei,
            "phase": "liquid_cocktail",
            "x_i": x[i],
            "DH_fus_kJmol": round(dh_fus, 3),
            "mu_E_kJmol": round(excess, 3),
            "RT_ln_x_kJmol": round(rt_ln, 3),
            "mu_total_SER_kJmol": round(total, 3),
        })
    # Pt solid substrate
    rows.append({
        "element": "Pt",
        "phase": "solid_fcc_pure",
        "x_i": 1.0,
        "DH_fus_kJmol": 0.0,
        "mu_E_kJmol": 0.0,
        "RT_ln_x_kJmol": 0.0,
        "mu_total_SER_kJmol": 0.0,
    })
    return rows


# ?? MAIN ????????????????????????????????????????????????????????????????
def main() -> None:
    elements = list(LIQUID_COMPOSITION.keys())
    x = np.array([LIQUID_COMPOSITION[e] for e in elements], dtype=float)
    assert abs(x.sum() - 1.0) < 1e-9

    M, omega_df = compute_omega_matrix(elements)

    # 0 K main output
    rows_0K = mu_liquid_SER(x, elements, M, T=0.0)
    df_0K = pd.DataFrame(rows_0K)

    # SI temperature table
    si_rows: list[dict] = []
    for T in SI_TEMPERATURES:
        for r in mu_liquid_SER(x, elements, M, float(T)):
            r["T_K"] = T
            si_rows.append(r)
    df_T = pd.DataFrame(si_rows)

    # Write
    out = Path(__file__).resolve().parent.parent / "outputs"
    out.mkdir(parents=True, exist_ok=True)
    df_0K.to_csv(out / "mu_liquid_v2_0K.csv", index=False, encoding="utf-8")
    df_T.to_csv(out / "mu_liquid_v2_T_table.csv", index=False, encoding="utf-8")
    omega_df.to_csv(out / "omega_liquid_binary_v2.csv", index=False, encoding="utf-8")
    meta = {
        "script": "script_A_mu_liquid_v2.py",
        "reference": "SER (0 K stable solid = 0)",
        "DH_fusion_source": "CRC Handbook 97th ed.",
        "composition": LIQUID_COMPOSITION,
        "miedema": {"P": P_MIEDEMA, "Q": Q_MIEDEMA},
        "v2_change": "Added DH_fusion offset; SER reference instead of pure-liquid=0",
    }
    (out / "script_A_v2_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    # Print summary
    print("=" * 70)
    print("Script A v2: liquid mu_i (SER reference)")
    print("=" * 70)
    print(df_0K[["element", "phase", "x_i", "DH_fus_kJmol",
                  "mu_E_kJmol", "mu_total_SER_kJmol"]].to_string(index=False))
    print(f"\nOutputs -> {out}")


if __name__ == "__main__":
    sys.exit(main() or 0)
