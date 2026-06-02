"""
regen_FigA_data.py
==================
Regenerate `data_FigA_v2_FamilyOrdered_Origin_regen.csv` from the Miedema
binary mixing-enthalpy formula using the element parameters in
`shared/data_periodic_table.py`.

Formula (manuscript Methods eq. 1):
    ΔH_AB(50:50) = −P (Δφ*)^2  +  Q (Δn_WS^(1/3))^2
        P = 14.1 kJ V^-2,  Q = 9.4 kJ (d.u.)^-2/3

Convention note: the `n_ws` field in `shared/data_periodic_table.py`
already stores n_WS^(1/3) directly (Pt = 1.78 = n_WS^(1/3), not the
raw electron density). Therefore Δ(n_WS^(1/3)) is computed as a plain
difference of the tabulated values, with no further cube root.

This script supersedes the externally-precomputed CSV that previously
shipped with this panel and was generated with an extra cube root on
the n_ws column (which shifted every binary ΔH by 1–2 kJ/mol). The
Pt–Ga value reproduces the manuscript Methods quote of −32.1 kJ/mol.
"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
from data_periodic_table import periodic_table_data as TBL  # noqa: E402

P_KJ_V2 = 14.1
Q_KJ_DU = 9.4

HOSTS = ["Pt", "Ir", "Pd", "Rh", "Ru", "Os",
         "Ni", "Co", "Fe", "Cr",
         "Re", "Au", "Zr", "W", "Hf", "Mo",
         "Ce", "La", "Y", "Sc"]
PARTNERS = ["Ga", "In", "Sn", "Zn", "Bi", "Hg"]


def dH_miedema(host: str, partner: str) -> float:
    a, b = TBL[host], TBL[partner]
    dphi = a["Phi"] - b["Phi"]
    dn = a["n_ws"] - b["n_ws"]            # tabulated n_ws IS n_WS^(1/3); use directly
    return -P_KJ_V2 * dphi ** 2 + Q_KJ_DU * dn ** 2


def main() -> None:
    rows = {p: {"Target_Element": p} for p in PARTNERS}
    for p in PARTNERS:
        for h in HOSTS:
            rows[p][h] = round(dH_miedema(h, p), 2)
    df = pd.DataFrame([rows[p] for p in PARTNERS])
    df = df[["Target_Element"] + HOSTS]
    out = HERE / "data_FigA_v2_FamilyOrdered_Origin_regen.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {out}")
    print("Canonical reference remains data_FigA_v2_FamilyOrdered_Origin.csv")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
