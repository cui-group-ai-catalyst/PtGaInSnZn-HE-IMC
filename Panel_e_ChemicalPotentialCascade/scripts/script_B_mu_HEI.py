"""
script_B_mu_HEI.py
==================
Compute partial molar chemical potentials of each element inside the final
HEI product Pt3(Ga,In,Sn,Zn) with L1_2 structure.

Structural model (Compound Energy Formalism, CEF; Hillert-Sundman-Agren 1997)
-----------------------------------------------------------------------------
- Sublattice alpha (3 sites / formula unit): fully occupied by Pt.
- Sublattice beta  (1 site  / formula unit): mixed (Ga 0.65, In 0.20,
                                              Sn 0.10, Zn 0.05).
- Total atoms per formula unit: 4.

Gibbs energy per mole formula unit (4 atoms), with pure-element reference:
    G_fu = sum_i y_i * dH_f^{Pt3i, f.u.}
         + R*T * sum_i y_i * ln(y_i)       # ideal mixing on beta
         + sum_{i<j} Omega^beta_{ij} * y_i * y_j    # beta excess

Per-atom chemical potentials (Hillert CEF identities, regular-solution
approximation on beta sublattice; pure-element 0 K reference):

mu_i^beta(HEI) per mol atom of i     =  dH_f^{Pt3i, f.u.} / 4 (per atom share)
                                      + R T ln y_i
                                      + sum_{j != i, j in beta} Omega^beta_ij * y_j^2
                                      [plus smaller ternary terms; see note]

mu_Pt^alpha(HEI) per mol atom Pt     =  sum_i y_i * dH_f^{Pt3i, f.u.} / 4
                                      (the "alpha-host" sees a weighted mix
                                       of end-member compounds; y_Pt^alpha=1
                                       so no ideal-mixing entropy term)

Euler consistency check (per mole formula unit, 4 atoms):
    G_fu  =?=  3 * mu_Pt^alpha + 1 * mu_i_weighted^beta        (to verify)

Data sources
------------
- End-member formation enthalpies dH_f^{Pt3i, f.u.} (L1_2):
    Pt3Ga : -52.1 kJ/mol f.u.   (Kumar, Liu, Chen 1996 J. Phase Equil.)
    Pt3In : -41.4 kJ/mol f.u.   (Srikanth, Petric 1993 CALPHAD)
    Pt3Sn : -62.4 kJ/mol f.u.   (Ghosh 2007 CALPHAD; Watson-Hayes 1995 calorimetry)
    Pt3Zn : -39.4 kJ/mol f.u.   (Liu et al. 2011 J. Alloys Compd.; Kumar 2004)
  (Estimated uncertainty ~ +/- 5 kJ/mol f.u.; will be propagated in Script E.)

- Beta-sublattice interaction Omega^beta_ij (Ga/In/Sn/Zn pairs):
    Inherited from Script A's Miedema Omega^L values as first-order estimate
    (Neumann-Kopp approximation: vibrational dS_mix ~ 0; solid/liquid Omega
    differ only by order 20% for isostructural p-block metals).

Convention
----------
- 0 K enthalpic approximation for main output (RT ln y term dropped).
- SI temperature table adds RT ln y with T in {300, 800, 1000, 1200 K}.
- Pure-element reference state (all mu_i^pure = 0).

Outputs
-------
- outputs/mu_HEI_0K.csv                 Main 5 numbers for Panel a
- outputs/mu_HEI_T_table.csv            SI temperature sensitivity
- outputs/omega_beta_subl.csv           Beta-sublattice Omega audit
- outputs/script_B_meta.json            Provenance + uncertainties

Author: 2026-04-23 (PtGaInSnZn submission prep)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# INPUTS
# --------------------------------------------------------------------------
# Composition on beta sublattice (mirrors Panel b liquid cocktail)
Y_BETA: dict[str, float] = {
    "Ga": 0.65,
    "In": 0.20,
    "Sn": 0.10,
    "Zn": 0.05,
}

# End-member L1_2 formation enthalpies per formula unit (4 atoms)
# Source citations in module docstring.
DH_F_PT3I_FU: dict[str, float] = {      # kJ/mol f.u.
    "Ga": -52.1,
    "In": -41.4,
    "Sn": -62.4,
    "Zn": -39.4,
}

# Estimated standard uncertainty per end-member (for SI error propagation)
DH_F_PT3I_SIGMA: dict[str, float] = {   # kJ/mol f.u.
    "Ga": 5.0,
    "In": 5.0,
    "Sn": 5.0,
    "Zn": 5.0,
}

# Element Miedema parameters (only for beta-subl. pair interactions;
# Pt parameters kept for consistency checks)
PERIODIC: dict[str, dict[str, float]] = {
    "Ga": {"Phi": 4.10, "n_ws": 1.34, "V": 11.8},
    "In": {"Phi": 3.90, "n_ws": 1.17, "V": 15.7},
    "Sn": {"Phi": 4.15, "n_ws": 1.25, "V": 16.3},
    "Zn": {"Phi": 4.10, "n_ws": 1.32, "V": 9.17},
    "Pt": {"Phi": 5.65, "n_ws": 1.78, "V": 9.10},
}

P_MIEDEMA: float = 14.1
Q_MIEDEMA: float = 9.4
R_KJ_MOL_K: float = 8.31446261815324e-3   # kJ/(mol*K)

# SI temperature sensitivity list
SI_TEMPERATURES: list[int] = [300, 800, 1000, 1200]

# Neumann-Kopp scaling for solid/liquid Omega  (1.0 = identical;
# typical p-block metals 0.8 - 1.0; we keep 1.0 and flag sensitivity)
NEUMANN_KOPP_FACTOR: float = 1.0

BETA_ELEMENTS: list[str] = ["Ga", "In", "Sn", "Zn"]


# --------------------------------------------------------------------------
# CORE FUNCTIONS
# --------------------------------------------------------------------------
def miedema_omega_simple(elem_i: str, elem_j: str) -> float:
    """Simplified Miedema Omega_ij in kJ/mol (no V prefactor, matches Panel d)."""
    pi = PERIODIC[elem_i]
    pj = PERIODIC[elem_j]
    d_phi = pi["Phi"] - pj["Phi"]
    d_n13 = pi["n_ws"] ** (1.0 / 3.0) - pj["n_ws"] ** (1.0 / 3.0)
    return -P_MIEDEMA * d_phi**2 + Q_MIEDEMA * d_n13**2


def omega_beta_matrix(elements: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    """Beta-sublattice Omega matrix for Ga/In/Sn/Zn pairs, plus audit DataFrame."""
    n = len(elements)
    M = np.zeros((n, n))
    rows: list[dict] = []
    for i, ei in enumerate(elements):
        for j, ej in enumerate(elements):
            if i == j:
                continue
            M[i, j] = NEUMANN_KOPP_FACTOR * miedema_omega_simple(ei, ej)
            if i < j:
                rows.append(
                    {
                        "i": ei,
                        "j": ej,
                        "Omega_beta_kJmol": round(M[i, j], 3),
                        "source": "Miedema simple x Neumann-Kopp factor",
                    }
                )
    return M, pd.DataFrame(rows)


def avg_endmember_per_atom(y: dict[str, float]) -> float:
    """Composition-weighted average end-member formation enthalpy, per mol atom.

    For each Pt3i (4 atoms), a "per-atom share" is dH_f_Pt3i / 4.
    The HEI has y_i moles of "Pt3i" flavour per mole beta; each atom (alpha
    or beta) carries this same per-atom share once we average over flavours.
    """
    s = 0.0
    for el, yi in y.items():
        s += yi * DH_F_PT3I_FU[el] / 4.0
    return s


def mu_HEI_0K(y: dict[str, float], M_beta: np.ndarray) -> dict[str, dict]:
    """Compute 0 K enthalpic chemical potentials (RT ln y -> 0).

    Returns dict keyed by element with breakdown fields:
        endmember_share_kJmol, excess_beta_kJmol, total_kJmol
    """
    avg_end = avg_endmember_per_atom(y)
    elements = list(y.keys())
    y_arr = np.array([y[e] for e in elements])

    result: dict[str, dict] = {}

    # Pt on alpha sublattice: endmember share only, no beta excess
    result["Pt"] = {
        "sublattice": "alpha",
        "y_or_x": 1.0,
        "endmember_share_kJmol": round(avg_end, 3),
        "RT_ln_y_kJmol": 0.0,
        "excess_beta_kJmol": 0.0,
        "total_kJmol": round(avg_end, 3),
        "note": "alpha site, y_Pt=1, no ideal entropy, no beta excess",
    }

    # Beta-sublattice elements (Ga, In, Sn, Zn)
    for i, ei in enumerate(elements):
        excess = 0.0
        for j, _ej in enumerate(elements):
            if j == i:
                continue
            excess += M_beta[i, j] * y_arr[j] ** 2
        total = avg_end + excess  # 0 K enthalpic; RT ln y = 0
        result[ei] = {
            "sublattice": "beta",
            "y_or_x": y[ei],
            "endmember_share_kJmol": round(avg_end, 3),
            "RT_ln_y_kJmol": 0.0,
            "excess_beta_kJmol": round(excess, 3),
            "total_kJmol": round(total, 3),
            "note": "beta site, 0K enthalpic",
        }
    return result


def mu_HEI_at_T(y: dict[str, float], M_beta: np.ndarray, T: float) -> dict[str, dict]:
    """Full chemical potentials at finite T, including ideal mixing."""
    avg_end = avg_endmember_per_atom(y)
    elements = list(y.keys())
    y_arr = np.array([y[e] for e in elements])
    result: dict[str, dict] = {}

    result["Pt"] = {
        "sublattice": "alpha",
        "y_or_x": 1.0,
        "endmember_share_kJmol": avg_end,
        "RT_ln_y_kJmol": 0.0,
        "excess_beta_kJmol": 0.0,
        "total_kJmol": avg_end,
    }

    for i, ei in enumerate(elements):
        excess = 0.0
        for j, _ej in enumerate(elements):
            if j == i:
                continue
            excess += M_beta[i, j] * y_arr[j] ** 2
        rt_ln_y = R_KJ_MOL_K * T * np.log(y[ei]) if y[ei] > 0 else float("nan")
        total = avg_end + rt_ln_y + excess
        result[ei] = {
            "sublattice": "beta",
            "y_or_x": y[ei],
            "endmember_share_kJmol": avg_end,
            "RT_ln_y_kJmol": rt_ln_y,
            "excess_beta_kJmol": excess,
            "total_kJmol": total,
        }
    return result


def euler_check(y: dict[str, float], mu_0K: dict[str, dict]) -> dict[str, float]:
    """Verify Euler relation (within the chosen CEF convention):

    G_fu/4 = 3/4 * mu_Pt^alpha/3 + 1/4 * sum_i y_i * mu_i^beta
           = 1/4 * mu_Pt^alpha + 1/4 * sum_i y_i * mu_i^beta ? (per atom)

    A cleaner check: per-atom average of (3*mu_Pt + 1*mu_beta_avg) / 4 should
    equal the composition-weighted per-atom enthalpy avg_endmember_per_atom
    (since Omega^beta excess averages to 0 when weighted by y_i for the
    regular-solution Hildebrand form -- a built-in consistency).
    """
    avg_end = avg_endmember_per_atom(y)
    mu_Pt = mu_0K["Pt"]["total_kJmol"]
    mu_beta_avg = sum(y[e] * mu_0K[e]["total_kJmol"] for e in y)
    per_atom_avg = (3.0 * mu_Pt + 1.0 * mu_beta_avg) / 4.0

    # Independent: sum of y_i * excess_i (should equal 2 * dH_mix^beta)
    excess_weighted = sum(y[e] * mu_0K[e]["excess_beta_kJmol"] for e in y)
    # Half of that is dH_mix^beta_per_atom; stored for info

    return {
        "avg_endmember_per_atom": avg_end,
        "mu_weighted_per_atom": per_atom_avg,
        "residual_kJmol": per_atom_avg - avg_end - 0.5 * excess_weighted,
        "excess_weighted": excess_weighted,
    }


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main() -> None:
    # Step 1: beta-sublattice Omega matrix
    M_beta, omega_df = omega_beta_matrix(BETA_ELEMENTS)

    # Step 2: 0 K enthalpic chemical potentials (main Panel a output)
    mu_0K = mu_HEI_0K(Y_BETA, M_beta)

    main_rows: list[dict] = []
    for el in ["Pt"] + BETA_ELEMENTS:
        r = mu_0K[el]
        main_rows.append(
            {
                "element": el,
                "sublattice": r["sublattice"],
                "y_or_x": r["y_or_x"],
                "endmember_share_kJmol": r["endmember_share_kJmol"],
                "RT_ln_y_kJmol": r["RT_ln_y_kJmol"],
                "excess_beta_kJmol": r["excess_beta_kJmol"],
                "mu_HEI_total_kJmol": r["total_kJmol"],
                "note": r["note"],
            }
        )
    df_main = pd.DataFrame(main_rows)

    # Step 3: SI temperature sensitivity
    si_rows: list[dict] = []
    for T in SI_TEMPERATURES:
        mu_T = mu_HEI_at_T(Y_BETA, M_beta, float(T))
        for el in ["Pt"] + BETA_ELEMENTS:
            r = mu_T[el]
            si_rows.append(
                {
                    "T_K": T,
                    "element": el,
                    "sublattice": r["sublattice"],
                    "y_or_x": r["y_or_x"],
                    "endmember_share_kJmol": round(r["endmember_share_kJmol"], 3),
                    "RT_ln_y_kJmol": round(r["RT_ln_y_kJmol"], 3),
                    "excess_beta_kJmol": round(r["excess_beta_kJmol"], 3),
                    "mu_HEI_total_kJmol": round(r["total_kJmol"], 3),
                }
            )
    df_T = pd.DataFrame(si_rows)

    # Step 4: Euler consistency check
    check = euler_check(Y_BETA, mu_0K)

    # Step 5: write outputs
    out_dir = Path(__file__).resolve().parent.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out_dir / "mu_HEI_0K.csv", index=False, encoding="utf-8")
    df_T.to_csv(out_dir / "mu_HEI_T_table.csv", index=False, encoding="utf-8")
    omega_df.to_csv(out_dir / "omega_beta_subl.csv", index=False, encoding="utf-8")

    meta = {
        "script": "script_B_mu_HEI.py",
        "composition_beta": Y_BETA,
        "endmember_dH_f_fu": DH_F_PT3I_FU,
        "endmember_dH_f_sigma": DH_F_PT3I_SIGMA,
        "endmember_citations": {
            "Pt3Ga": "Kumar, Liu, Chen 1996, J. Phase Equilib. 17, 482",
            "Pt3In": "Srikanth, Petric 1993, CALPHAD 17, 39",
            "Pt3Sn": "Ghosh 2007 CALPHAD; Watson, Hayes 1995 calorimetry",
            "Pt3Zn": "Liu et al. 2011, J. Alloys Compd.; Kumar 2004",
        },
        "miedema": {"P": P_MIEDEMA, "Q": Q_MIEDEMA, "neumann_kopp_factor": NEUMANN_KOPP_FACTOR},
        "temperatures_K": SI_TEMPERATURES,
        "convention": "per mole atom of species; 0 K enthalpic main; Hillert CEF on L1_2",
        "euler_check": check,
    }
    (out_dir / "script_B_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # Step 6: print summary
    print("=" * 72)
    print("Script B: HEI chemical potentials mu_i^HEI (0 K enthalpic)")
    print("=" * 72)
    print(f"beta composition: {Y_BETA}")
    print(f"End-member dH_f (kJ/mol f.u.): {DH_F_PT3I_FU}")
    print(f"Weighted avg endmember share (per atom): "
          f"{avg_endmember_per_atom(Y_BETA):.3f} kJ/mol")
    print()
    print("-- Omega^beta (Ga/In/Sn/Zn pairs, kJ/mol) --")
    print(omega_df.to_string(index=False))
    print()
    print("-- Main output: mu_i^HEI @ 0 K (for Panel a) --")
    print(
        df_main[
            [
                "element", "sublattice", "y_or_x",
                "endmember_share_kJmol", "excess_beta_kJmol",
                "mu_HEI_total_kJmol",
            ]
        ].to_string(index=False)
    )
    print()
    print("-- Euler consistency check --")
    for k, v in check.items():
        print(f"  {k:35s} : {v:+.4f}" if isinstance(v, float) else f"  {k}: {v}")
    print()
    print(f"Outputs written to: {out_dir}")


if __name__ == "__main__":
    sys.exit(main() or 0)
