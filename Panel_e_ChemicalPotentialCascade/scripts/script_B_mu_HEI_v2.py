"""
script_B_mu_HEI_v2.py
=====================
Compute per-atom chemical potentials inside the HEI product
Pt3(Ga,In,Sn,Zn) with L1_2 structure, using SER reference.

v2 change vs. v1
-----------------
v1 gave ALL beta elements the same mu (composition-weighted average of
end-member enthalpies). v2 correctly assigns each beta element its OWN
Pt3i end-member formation enthalpy, per the CEF partial derivative.

CEF sublattice model
--------------------
  alpha (3 sites / f.u.): Pt only, y_Pt^alpha = 1
  beta  (1 site  / f.u.): Ga/In/Sn/Zn mix, y_i^beta variable

Per formula unit (4 atoms):
  G_fu = sum_i y_i * dH_f(Pt3i) + RT sum_i y_i ln y_i + G^xs_beta

End-member chemical potential on beta sublattice (per f.u.):
  mu_{Pt3i}^end = dH_f(Pt3i) + RT ln y_i + excess_beta_i

For per-atom representation on Panel a y-axis:
  mu_i^HEI (per atom) = dH_f(Pt3i) / 4  +  (RT ln y_i) / 4  +  excess / 4
  mu_Pt^HEI (per atom) = sum_i y_i dH_f(Pt3i) / 4  (weighted avg, no entropy)

End-member data (kJ/mol f.u., SER reference)
---------------------------------------------
  Pt3Ga: -52.1   Kumar, Liu, Chen 1996, J. Phase Equilib.
  Pt3In: -41.4   Srikanth, Petric 1993, CALPHAD
  Pt3Sn: -62.4   Watson, Hayes 1995 calorimetry; Ghosh 2007
  Pt3Zn: -39.4   Liu et al. 2011, J. Alloys Compd.

Outputs
-------
  mu_HEI_v2_0K.csv          Main Panel a end-points (per atom, SER ref)
  mu_HEI_v2_T_table.csv     SI temperature table
  omega_beta_subl_v2.csv     Beta-sublattice Omega audit
  script_B_v2_meta.json

Author: 2026-04-23 v2
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ?? INPUTS ??????????????????????????????????????????????????????????????
Y_BETA: dict[str, float] = {
    "Ga": 0.65, "In": 0.20, "Sn": 0.10, "Zn": 0.05,
}

DH_F_PT3I_FU: dict[str, float] = {   # kJ/mol f.u. (SER reference)
    "Ga": -52.1, "In": -41.4, "Sn": -62.4, "Zn": -39.4,
}

PERIODIC: dict[str, dict[str, float]] = {
    "Ga": {"Phi": 4.10, "n_ws": 1.34},
    "In": {"Phi": 3.90, "n_ws": 1.17},
    "Sn": {"Phi": 4.15, "n_ws": 1.25},
    "Zn": {"Phi": 4.10, "n_ws": 1.32},
}

P_MIEDEMA, Q_MIEDEMA = 14.1, 9.4
R_KJ = 8.31446261815324e-3
SI_TEMPERATURES: list[int] = [300, 800, 1000, 1200]
BETA_ELEMENTS: list[str] = ["Ga", "In", "Sn", "Zn"]

# ?? CORE ????????????????????????????????????????????????????????????????
def miedema_omega(ei: str, ej: str) -> float:
    pi, pj = PERIODIC[ei], PERIODIC[ej]
    return -P_MIEDEMA * (pi["Phi"] - pj["Phi"])**2 + \
            Q_MIEDEMA * (pi["n_ws"]**(1/3) - pj["n_ws"]**(1/3))**2


def omega_matrix(els: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    n = len(els)
    M = np.zeros((n, n))
    rows = []
    for i, ei in enumerate(els):
        for j, ej in enumerate(els):
            if i == j:
                continue
            M[i, j] = miedema_omega(ei, ej)
            if i < j:
                rows.append({"i": ei, "j": ej,
                             "Omega_beta_kJmol": round(M[i, j], 3)})
    return M, pd.DataFrame(rows)


def mu_HEI_per_atom(y: dict[str, float], M: np.ndarray, T: float) -> list[dict]:
    """Return per-atom chemical potentials (SER ref) for 5 elements."""
    els = list(y.keys())
    y_arr = np.array([y[e] for e in els])
    rows: list[dict] = []

    # Weighted average end-member enthalpy per atom (for Pt alpha)
    avg_dHf_per_atom = sum(y[e] * DH_F_PT3I_FU[e] for e in els) / 4.0

    # Pt on alpha sublattice: sees weighted average, no mixing entropy
    rows.append({
        "element": "Pt",
        "sublattice": "alpha",
        "y_or_x": 1.0,
        "dHf_endmember_per_atom_kJmol": round(avg_dHf_per_atom, 3),
        "RT_lny_per_atom_kJmol": 0.0,
        "excess_per_atom_kJmol": 0.0,
        "mu_HEI_per_atom_kJmol": round(avg_dHf_per_atom, 3),
        "note": "alpha: weighted avg of Pt3i end-members / 4",
    })

    # Beta sublattice elements: each uses its OWN Pt3i
    for idx, ei in enumerate(els):
        dHf_per_atom = DH_F_PT3I_FU[ei] / 4.0
        rt_lny = (R_KJ * T * np.log(y[ei]) / 4.0) if T > 0 else 0.0
        excess = sum(M[idx, j] * y_arr[j]**2 for j in range(len(els)) if j != idx) / 4.0
        total = dHf_per_atom + rt_lny + excess
        rows.append({
            "element": ei,
            "sublattice": "beta",
            "y_or_x": y[ei],
            "dHf_endmember_per_atom_kJmol": round(dHf_per_atom, 3),
            "RT_lny_per_atom_kJmol": round(rt_lny, 3),
            "excess_per_atom_kJmol": round(excess, 3),
            "mu_HEI_per_atom_kJmol": round(total, 3),
            "note": f"beta: Pt3{ei} end-member / 4",
        })
    return rows


# ?? MAIN ????????????????????????????????????????????????????????????????
def main() -> None:
    M, omega_df = omega_matrix(BETA_ELEMENTS)

    # 0 K main
    rows_0K = mu_HEI_per_atom(Y_BETA, M, T=0.0)
    df_0K = pd.DataFrame(rows_0K)

    # SI temperature table
    si_rows = []
    for T in SI_TEMPERATURES:
        for r in mu_HEI_per_atom(Y_BETA, M, float(T)):
            r["T_K"] = T
            si_rows.append(r)
    df_T = pd.DataFrame(si_rows)

    # Write
    out = Path(__file__).resolve().parent.parent / "outputs"
    out.mkdir(parents=True, exist_ok=True)
    df_0K.to_csv(out / "mu_HEI_v2_0K.csv", index=False, encoding="utf-8")
    df_T.to_csv(out / "mu_HEI_v2_T_table.csv", index=False, encoding="utf-8")
    omega_df.to_csv(out / "omega_beta_subl_v2.csv", index=False, encoding="utf-8")

    meta = {
        "script": "script_B_mu_HEI_v2.py",
        "reference": "SER (0 K stable solid = 0)",
        "composition_beta": Y_BETA,
        "endmember_dH_f_fu": DH_F_PT3I_FU,
        "endmember_citations": {
            "Pt3Ga": "Kumar, Liu, Chen 1996",
            "Pt3In": "Srikanth, Petric 1993",
            "Pt3Sn": "Watson-Hayes 1995 + Ghosh 2007",
            "Pt3Zn": "Liu et al. 2011",
        },
        "v2_change": "beta elements use individual Pt3i end-members, not weighted avg",
        "convention": "per atom of f.u. (divide f.u. quantities by 4)",
    }
    (out / "script_B_v2_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    print("=" * 70)
    print("Script B v2: HEI mu_i (SER ref, per atom, end-member resolved)")
    print("=" * 70)
    print(df_0K[["element", "sublattice", "y_or_x",
                  "dHf_endmember_per_atom_kJmol",
                  "excess_per_atom_kJmol",
                  "mu_HEI_per_atom_kJmol", "note"]].to_string(index=False))
    print(f"\nOutputs -> {out}")


if __name__ == "__main__":
    sys.exit(main() or 0)
