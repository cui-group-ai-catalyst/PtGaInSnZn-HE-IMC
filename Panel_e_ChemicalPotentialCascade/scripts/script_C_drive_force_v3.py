# -*- coding: utf-8 -*-
"""
script_C_drive_force_v3.py
==========================
Combine SER-referenced liquid (A_v3) and UMA-CEF HEI (B_v3) chemical
potentials into per-element driving forces and total reaction energy.

v3 changes vs v2
-----------------
- HEI chemical potentials now come from CEF fitted to Panel g 165-pt UMA
  data (script_B_mu_HEI_v3), replacing literature values.
- Liquid side unchanged (script_A_mu_liquid_v3, same as v2).

Legacy element-resolved decomposition (manuscript eq. 3, SI Table 1):
    F_i = mu_i_source - mu_i_HEI       (kJ/mol per atom of species i)
        mu_i_source : Pt = 0 (SER); X = liquid cocktail chemical potential
        mu_i_HEI    : chemical potential of i inside the HEI (CEF fit, per atom)

The fixed-composition Pt3X energy manifold does not uniquely identify every
absolute elemental mu_i. The F_i values are retained as an Euler-consistent
reference-convention decomposition for manuscript regression. The total
reaction energy below is the primary gauge-invariant thermodynamic result;
beta-sublattice diffusion potentials are written by script_B_mu_HEI_v3.py.

Total reaction Gibbs energy per formula unit (primary result):
    dG_rxn = 4*G_atom_HEI - Sum y_i * mu_i^L            (negative = spontaneous)
           = -Sum nu_i * F_i,    nu_Pt = 3,  nu_X = y_X (manuscript eq 3)

Outputs
-------
  delta_mu_v3_0K.csv            Per-element driving forces
  delta_G_rxn_v3_summary.csv    Total reaction summary
  panel_a_v3_tier_summary.csv   Three-tier data for plotting
  script_C_v3_meta.json

Author: 2026-04-23 v3
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

Y_COMP: dict[str, float] = {"Ga": 0.65, "In": 0.20, "Sn": 0.10, "Zn": 0.05}
BETA_ELEMENTS = ["Ga", "In", "Sn", "Zn"]


def main() -> None:
    out = Path(__file__).resolve().parent.parent / "outputs"
    df_L = pd.read_csv(out / "mu_liquid_v3_0K.csv")
    df_H = pd.read_csv(out / "mu_HEI_v3_0K.csv")

    mu_start = {r["element"]: r["mu_total_SER_kJmol"] for _, r in df_L.iterrows()}
    mu_end = {r["element"]: r["mu_HEI_per_atom_kJmol"] for _, r in df_H.iterrows()}

    dmu_rows: list[dict] = []
    for el in BETA_ELEMENTS:
        # F_i per manuscript eq (3) and SI Table 1:
        #   F_i = mu_i_source - mu_i_HEI   (per atom of species i)
        # (previous v3 used (3*mu_Pt + mu_i)/4 as the "start", which mixed in
        #  Pt and therefore did not reproduce SI Table 1's F_i column.)
        start_pa = mu_start[el]
        end_pa = mu_end[el]
        F = start_pa - end_pa
        dmu_rows.append({
            "element": el,
            "sublattice": "beta",
            "y_i": Y_COMP[el],
            "mu_start_liquid_kJmol": round(mu_start[el], 4),
            "start_per_atom_kJmol": round(start_pa, 4),
            "mu_HEI_per_atom_kJmol": round(end_pa, 4),
            "driving_force_per_atom_kJmol": round(F, 4),
            "driving_force_per_fu_kJmol": round(F * 4, 4),
            "interpretation": (
                "legacy Euler-consistent decomposition; not a unique absolute "
                "elemental chemical potential"
            ),
        })

    pt_start = mu_start["Pt"]
    pt_end = mu_end["Pt"]
    F_Pt = pt_start - pt_end
    dmu_rows.append({
        "element": "Pt",
        "sublattice": "alpha",
        "y_i": 1.0,
        "mu_start_liquid_kJmol": round(pt_start, 4),
        "start_per_atom_kJmol": round(pt_start, 4),
        "mu_HEI_per_atom_kJmol": round(pt_end, 4),
        "driving_force_per_atom_kJmol": round(F_Pt, 4),
        "driving_force_per_fu_kJmol": round(F_Pt * 4, 4),
        "interpretation": (
            "legacy Euler-consistent decomposition; not a unique absolute "
            "elemental chemical potential"
        ),
    })
    dmu_df = pd.DataFrame(dmu_rows)

    G_atom_HEI = mu_end["Pt"]
    sum_y_mu_L = sum(Y_COMP[el] * mu_start[el] for el in BETA_ELEMENTS)
    dG_rxn_fu = 4.0 * G_atom_HEI - sum_y_mu_L
    dG_rxn_atom = dG_rxn_fu / 4.0

    dG_summary = {
        "delta_G_rxn_per_fu_kJmol": round(dG_rxn_fu, 4),
        "delta_G_rxn_per_atom_kJmol": round(dG_rxn_atom, 4),
        "G_atom_HEI_kJmol": round(G_atom_HEI, 4),
        "sum_y_mu_L_kJmol": round(sum_y_mu_L, 4),
        "notes": (
            "delta_G_rxn<0 is thermodynamically favourable within the stated "
            "model and reference states; kinetics and competing phases are not implied"
        ),
    }

    tier_rows: list[dict] = []
    for el in BETA_ELEMENTS:
        start_pa = (3 * mu_start["Pt"] + mu_start[el]) / 4.0
        tier_rows.append({
            "tier": "liquid_top",
            "element": el,
            "y_or_x": Y_COMP[el],
            "mu_per_atom_kJmol": round(start_pa, 4),
            "phase_label": f"{el} (liquid cocktail)",
        })
    tier_rows.append({
        "tier": "Pt_middle",
        "element": "Pt",
        "y_or_x": 1.0,
        "mu_per_atom_kJmol": 0.0,
        "phase_label": "Pt (s, fcc, SER=0)",
    })
    for el in ["Pt"] + BETA_ELEMENTS:
        tier_rows.append({
            "tier": "HEI_bottom",
            "element": el,
            "y_or_x": 1.0 if el == "Pt" else Y_COMP[el],
            "mu_per_atom_kJmol": round(mu_end[el], 4),
            "phase_label": "Pt3(Ga,In,Sn,Zn) L1_2",
        })
    tier_df = pd.DataFrame(tier_rows)

    dmu_df.to_csv(out / "delta_mu_v3_0K.csv", index=False, encoding="utf-8")
    pd.DataFrame([dG_summary]).to_csv(
        out / "delta_G_rxn_v3_summary.csv", index=False, encoding="utf-8")
    tier_df.to_csv(out / "panel_a_v3_tier_summary.csv", index=False, encoding="utf-8")

    meta = {
        "script": "script_C_drive_force_v3.py",
        "convention": "per atom of f.u.; SER reference; 0K enthalpic",
        "interpretation_priority": (
            "delta_G_rxn is primary and gauge invariant within the model; "
            "elementwise F_i is a reference-convention decomposition"
        ),
        "beta_diffusion_potentials": "beta_diffusion_potentials_v3_0K.csv",
        "data_sources": {
            "liquid": "script_A_mu_liquid_v3 (Miedema + CRC DH_fus)",
            "HEI": "script_B_mu_HEI_v3 (CEF fitted to Panel g 165 UMA points)",
        },
        "driving_force_summary": dG_summary,
        "v3_change": "HEI from UMA-CEF instead of literature values",
    }
    (out / "script_C_v3_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    print("=" * 70)
    print("Script C v3: element-resolved driving forces (per atom, SER)")
    print("=" * 70)
    print(dmu_df.to_string(index=False))
    print()
    print("-- Total reaction --")
    for k, v in dG_summary.items():
        if isinstance(v, str):
            print(f"  {k:40s}: {v}")
        else:
            print(f"  {k:40s}: {v:+.4f} kJ/mol")
    print()
    print("-- Tier summary for Panel a --")
    print(tier_df.to_string(index=False))
    print(f"\nOutputs -> {out}")


if __name__ == "__main__":
    main()
