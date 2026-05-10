"""
script_C_drive_force_v2.py
==========================
Combine SER-referenced liquid (A_v2) and HEI (B_v2) chemical potentials
into per-element driving forces and total reaction Gibbs energy.

v2 changes
----------
- Liquid starting mu now includes DH_fusion (SER reference).
- HEI end mu uses individual Pt3i end-members for beta elements.
- All quantities are per atom of formula unit (kJ/mol).

Transfer event for beta element i:
    3 Pt(s, SER=0) + 1 i(l, SER=+DH_fus) -> Pt3i(HEI, SER=dHf/4)

Per-atom driving force:
    F_i = [start_per_atom] - [end_per_atom]
        = (3*mu_Pt + 1*mu_i^L) / 4 - dH_f(Pt3i) / 4
    (positive = spontaneous)

Total reaction driving force (weighted by beta composition):
    delta_G_rxn = -sum_i y_i * F_i * 4   (per f.u., negative = spontaneous)

Outputs
-------
  delta_mu_v2_0K.csv              Per-element driving forces
  delta_G_rxn_v2_summary.csv      Total reaction summary
  panel_a_v2_tier_summary.csv     Three-tier data for plotting
  script_C_v2_meta.json

Author: 2026-04-23 v2
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

Y_COMP: dict[str, float] = {"Ga": 0.65, "In": 0.20, "Sn": 0.10, "Zn": 0.05}
BETA_ELEMENTS = ["Ga", "In", "Sn", "Zn"]


def main() -> None:
    out = Path(__file__).resolve().parent.parent / "outputs"
    df_L = pd.read_csv(out / "mu_liquid_v2_0K.csv")
    df_H = pd.read_csv(out / "mu_HEI_v2_0K.csv")

    mu_start = {r["element"]: r["mu_total_SER_kJmol"] for _, r in df_L.iterrows()}
    mu_end = {r["element"]: r["mu_HEI_per_atom_kJmol"] for _, r in df_H.iterrows()}

    # Per-element driving forces
    dmu_rows: list[dict] = []
    for el in BETA_ELEMENTS:
        start_per_atom = (3 * mu_start["Pt"] + mu_start[el]) / 4.0
        end_per_atom = mu_end[el]
        F = start_per_atom - end_per_atom
        dmu_rows.append({
            "element": el,
            "mu_start_liquid_kJmol": round(mu_start[el], 3),
            "start_per_atom_kJmol": round(start_per_atom, 3),
            "mu_HEI_per_atom_kJmol": round(end_per_atom, 3),
            "driving_force_per_atom_kJmol": round(F, 3),
            "driving_force_per_fu_kJmol": round(F * 4, 3),
        })
    # Pt: start = 0, end = weighted avg
    pt_start = mu_start["Pt"]
    pt_end = mu_end["Pt"]
    F_Pt = pt_start - pt_end
    dmu_rows.append({
        "element": "Pt",
        "mu_start_liquid_kJmol": round(pt_start, 3),
        "start_per_atom_kJmol": round(pt_start, 3),
        "mu_HEI_per_atom_kJmol": round(pt_end, 3),
        "driving_force_per_atom_kJmol": round(F_Pt, 3),
        "driving_force_per_fu_kJmol": round(F_Pt * 4, 3),
    })
    dmu_df = pd.DataFrame(dmu_rows)

    # Total delta_G_rxn: weighted sum of transfer events
    # Each beta element contributes y_i fraction of transfers
    # Per f.u.: delta_G_rxn = -sum_i y_i * F_i * 4  (F_i in per-atom)
    total_F_per_fu = sum(
        Y_COMP[el] * dmu_df[dmu_df["element"] == el]["driving_force_per_fu_kJmol"].values[0]
        for el in BETA_ELEMENTS
    )
    total_F_per_atom = total_F_per_fu / 4.0

    dG_summary = {
        "delta_G_rxn_per_fu_kJmol": round(-total_F_per_fu, 3),
        "delta_G_rxn_per_atom_kJmol": round(-total_F_per_atom, 3),
        "driving_force_per_fu_kJmol": round(total_F_per_fu, 3),
        "driving_force_per_atom_kJmol": round(total_F_per_atom, 3),
        "notes": "delta_G_rxn<0 => spontaneous; driving_force>0 => spontaneous",
    }

    # Tier summary for plotting
    tier_rows: list[dict] = []
    for el in BETA_ELEMENTS:
        start_pa = (3 * mu_start["Pt"] + mu_start[el]) / 4.0
        tier_rows.append({
            "tier": "liquid_top",
            "element": el,
            "y_or_x": Y_COMP[el],
            "mu_per_atom_kJmol": round(start_pa, 3),
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
            "mu_per_atom_kJmol": round(mu_end[el], 3),
            "phase_label": "Pt3(Ga,In,Sn,Zn) L1_2",
        })
    tier_df = pd.DataFrame(tier_rows)

    # Write
    dmu_df.to_csv(out / "delta_mu_v2_0K.csv", index=False, encoding="utf-8")
    pd.DataFrame([dG_summary]).to_csv(
        out / "delta_G_rxn_v2_summary.csv", index=False, encoding="utf-8")
    tier_df.to_csv(out / "panel_a_v2_tier_summary.csv", index=False, encoding="utf-8")
    meta = {
        "script": "script_C_drive_force_v2.py",
        "convention": "per atom of f.u.; SER reference; 0K enthalpic",
        "driving_force_summary": dG_summary,
        "v2_change": "SER liquid + end-member-resolved HEI = element-specific driving forces",
    }
    (out / "script_C_v2_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    # Print
    print("=" * 70)
    print("Script C v2: element-resolved driving forces (per atom, SER)")
    print("=" * 70)
    print(dmu_df.to_string(index=False))
    print()
    print("-- Total reaction --")
    for k, v in dG_summary.items():
        if isinstance(v, str):
            print(f"  {k:40s}: {v}")
        else:
            print(f"  {k:40s}: {v:+.3f} kJ/mol")
    print()
    print("-- Tier summary for Panel a --")
    print(tier_df.to_string(index=False))
    print(f"\nOutputs -> {out}")


if __name__ == "__main__":
    main()
