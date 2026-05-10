"""
script_C_drive_force.py
=======================
Combine liquid-phase (Script A) and HEI-phase (Script B) chemical potentials
into:

  (1) Per-element driving forces     Delta_mu_i = mu_i^(0) - mu_i^HEI
  (2) Reaction total driving force   Delta_G_rxn (per mole atom of HEI)
  (3) Panel a visualization-ready summary table (three tiers)

Reaction (normalized to 1 mole f.u. HEI, 4 atoms):
    3 Pt(s)  +  Ga_0.65 In_0.20 Sn_0.10 Zn_0.05 (liquid)  ->  Pt3M(L1_2)

Per-atom total driving force (exhaustion-type termination, y_i^beta = x_i^L):
    Delta_G_rxn  =  (3 * mu_Pt^solid + sum_i y_i mu_i^L) / 4
                  -  (3 * mu_Pt^alpha + sum_i y_i mu_i^beta) / 4
                  =  0  -  avg_endmember_per_atom  -  <excess_beta>/2    (since RT ln y drops at 0 K)

For Panel a "height drop" visualization, we need per-element Delta_mu_i
                                           and the overall Delta_G_rxn.

Inputs (read from Script A & B outputs):
-----------------------------------------
- outputs/mu_liquid_0K.csv    (Script A main output)
- outputs/mu_HEI_0K.csv       (Script B main output)

Outputs:
--------
- outputs/delta_mu_0K.csv           Per-element Delta_mu (Panel a heights)
- outputs/delta_G_rxn_summary.csv   One-row summary with total driving force
- outputs/panel_a_tier_summary.csv  Three-tier table (liquid / Pt-solid / HEI)
- outputs/script_C_meta.json

Author: 2026-04-23 (PtGaInSnZn submission prep)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


# Composition (identical to Scripts A & B)
Y_COMPOSITION: dict[str, float] = {
    "Ga": 0.65,
    "In": 0.20,
    "Sn": 0.10,
    "Zn": 0.05,
}
BETA_ELEMENTS: list[str] = ["Ga", "In", "Sn", "Zn"]

# Reaction stoichiometric counts per formula unit HEI (4 atoms)
N_PT_PER_FU: int = 3
N_BETA_PER_FU: int = 1


def load_script_outputs(out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_L = pd.read_csv(out_dir / "mu_liquid_0K.csv")
    df_H = pd.read_csv(out_dir / "mu_HEI_0K.csv")
    return df_L, df_H


def build_mu_maps(df_L: pd.DataFrame, df_H: pd.DataFrame) -> tuple[dict, dict]:
    mu_start = {row["element"]: row["mu_total_kJmol"] for _, row in df_L.iterrows()}
    mu_hei = {row["element"]: row["mu_HEI_total_kJmol"] for _, row in df_H.iterrows()}
    return mu_start, mu_hei


def compute_delta_mu(mu_start: dict[str, float], mu_hei: dict[str, float]) -> pd.DataFrame:
    rows: list[dict] = []
    for el in ["Pt"] + BETA_ELEMENTS:
        start = mu_start[el]
        end = mu_hei[el]
        rows.append(
            {
                "element": el,
                "mu_start_kJmol": round(start, 3),
                "mu_HEI_kJmol": round(end, 3),
                "delta_mu_kJmol": round(start - end, 3),
                "start_phase": "solid_fcc_pure" if el == "Pt" else "liquid_cocktail",
                "end_phase": "HEI_L1_2",
            }
        )
    return pd.DataFrame(rows)


def compute_delta_G_rxn(delta_mu_df: pd.DataFrame, y: dict[str, float]) -> dict:
    """Per-atom driving force AND thermodynamic Delta_G_rxn.

    Driving force (positive = spontaneous to proceed):
        F^{f.u.}   =  3 * Delta_mu_Pt  +  1 * sum_i y_i Delta_mu_i   (kJ/mol f.u.)
        F^{atom}   =  F^{f.u.} / 4                                   (kJ/mol atom)

    Thermodynamic Delta_G_rxn (negative = spontaneous):
        Delta_G_rxn = -F                                              (kJ/mol)
    """
    row_by_el = {r["element"]: r for _, r in delta_mu_df.iterrows()}
    dmu_Pt = row_by_el["Pt"]["delta_mu_kJmol"]
    dmu_beta_weighted = sum(y[e] * row_by_el[e]["delta_mu_kJmol"] for e in y)

    F_per_fu = N_PT_PER_FU * dmu_Pt + N_BETA_PER_FU * dmu_beta_weighted
    F_per_atom = F_per_fu / (N_PT_PER_FU + N_BETA_PER_FU)

    return {
        "delta_mu_Pt_kJmol": round(dmu_Pt, 3),
        "delta_mu_beta_weighted_kJmol": round(dmu_beta_weighted, 3),
        "driving_force_per_fu_kJmol": round(F_per_fu, 3),
        "driving_force_per_atom_kJmol": round(F_per_atom, 3),
        "delta_G_rxn_per_fu_kJmol": round(-F_per_fu, 3),
        "delta_G_rxn_per_atom_kJmol": round(-F_per_atom, 3),
        "notes": (
            "Driving force F = sum(stoich * Delta_mu) > 0 means reaction proceeds "
            "spontaneously; Delta_G_rxn = -F is the thermodynamic free-energy "
            "change of the reaction (should be < 0 for spontaneity)."
        ),
    }


def build_panel_a_tier_table(
    mu_start: dict[str, float], mu_hei: dict[str, float], y: dict[str, float]
) -> pd.DataFrame:
    """Three-tier table: liquid / Pt-solid starting / HEI end.

    For Panel a, we need three heights:
      * Liquid tier (top-left)       : elemental mu_i^L for Ga/In/Sn/Zn only.
      * Pt-substrate tier (middle)   : mu_Pt^solid = 0.
      * HEI tier (bottom-right)      : mu_i^HEI for all 5 elements.
    """
    rows: list[dict] = []

    for el in BETA_ELEMENTS:
        rows.append(
            {
                "tier": "liquid_top",
                "element": el,
                "y_or_x": y[el],
                "mu_kJmol": round(mu_start[el], 3),
                "phase_label": "Ga/In/Sn/Zn (l, cocktail)",
            }
        )
    rows.append(
        {
            "tier": "Pt_middle",
            "element": "Pt",
            "y_or_x": 1.0,
            "mu_kJmol": round(mu_start["Pt"], 3),
            "phase_label": "Pt (s, fcc)",
        }
    )
    for el in ["Pt"] + BETA_ELEMENTS:
        rows.append(
            {
                "tier": "HEI_bottom",
                "element": el,
                "y_or_x": 1.0 if el == "Pt" else y[el],
                "mu_kJmol": round(mu_hei[el], 3),
                "phase_label": "Pt3(Ga,In,Sn,Zn) L1_2",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "outputs"
    df_L, df_H = load_script_outputs(out_dir)
    mu_start, mu_hei = build_mu_maps(df_L, df_H)

    # (1) Delta_mu table per element
    dmu_df = compute_delta_mu(mu_start, mu_hei)

    # (2) total reaction driving force
    dG = compute_delta_G_rxn(dmu_df, Y_COMPOSITION)

    # (3) tier table for Panel a
    tier_df = build_panel_a_tier_table(mu_start, mu_hei, Y_COMPOSITION)

    # write
    dmu_df.to_csv(out_dir / "delta_mu_0K.csv", index=False, encoding="utf-8")
    pd.DataFrame([dG]).to_csv(out_dir / "delta_G_rxn_summary.csv", index=False, encoding="utf-8")
    tier_df.to_csv(out_dir / "panel_a_tier_summary.csv", index=False, encoding="utf-8")

    meta = {
        "script": "script_C_drive_force.py",
        "composition": Y_COMPOSITION,
        "reaction_per_fu": f"3 Pt(s) + GaInSnZn (liq) -> Pt3(Ga{Y_COMPOSITION['Ga']}In{Y_COMPOSITION['In']}Sn{Y_COMPOSITION['Sn']}Zn{Y_COMPOSITION['Zn']}) (L1_2)",
        "termination_type": "exhaustion (liquid limiting reagent, y_i^beta = x_i^L)",
        "convention": "per mole atom of species; 0 K enthalpic",
        "driving_force_summary": dG,
    }
    (out_dir / "script_C_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # print
    print("=" * 72)
    print("Script C: Driving forces Delta_mu_i and Delta_G_rxn")
    print("=" * 72)
    print(f"Composition: {Y_COMPOSITION}")
    print()
    print("-- Delta_mu_i per element (kJ/mol atom of species) --")
    print(dmu_df.to_string(index=False))
    print()
    print("-- Total reaction driving force --")
    for k, v in dG.items():
        if isinstance(v, str):
            print(f"  {k:35s} : {v}")
        else:
            print(f"  {k:38s} : {v:+.3f} kJ/mol")
    print()
    print("-- Panel a three-tier summary --")
    print(tier_df.to_string(index=False))
    print()
    print(f"Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
