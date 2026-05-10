"""
20260419_FigG_Origin_Ready.py
=============================
One-shot Origin-ready wide CSV for Panel g, mirroring the proven layout of
  Fig1h_Origin_AllPanels_PreciseNames_Hf.csv
so the user can drop it into an existing Origin template and plot directly.

Layout (per Ga% bin block):
  {GA}%_{CAT}_X, {GA}%_{CAT}_Y        -- one pair per category present in that bin
  {GA}%_Mean_X, {GA}%_Mean_Y, {GA}%_Mean_Std   -- one mean + std per bin

  - X jitter: uniform(-0.35, 0.35) centered on 0 (deterministic, seed=42)
  - Y:        ElementRef_Hf_kJ_mol (element-referenced formation enthalpy)
  - Mean_X = 0.0, Mean_Y = bin mean, Mean_Std = bin std (population, ddof=0)
  - Category slugs: GaSn, GaIn, InSn, GaInSn, FiveElem (Pt+Ga+In+Sn+Zn),
    PartZn (has Zn but not all 4 non-Pt elements), Pure (single non-Pt element)
    -- identical to Fig1h_..._PreciseNames_Hf.csv so existing Origin templates
    keyed on those slugs work without editing.

Input:  data_FigG_165_ElementReferenced_Hf.csv
Output: written to the panel directory
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR

SRC = SCRIPT_DIR / "data_FigG_165_ElementReferenced_Hf.csv"
OUT = SCRIPT_DIR / "data_FigG_Origin_Ready_regen.csv"

JITTER_HALFWIDTH = 0.35
SEED = 42

# Category slugs must match Fig1h_..._PreciseNames_Hf.csv for template reuse
CAT_ORDER = ["GaSn", "GaIn", "InSn", "GaInSn", "FiveElem", "PartZn", "Pure"]


def classify(ng: int, ni: int, ns: int, nz: int) -> str:
    """Map (Ga,In,Sn,Zn) atom counts on the 8-site B-sublattice to a slug.
    Matches the scheme used by Fig1h_..._PreciseNames_Hf.csv.
    """
    has = {"Ga": ng > 0, "In": ni > 0, "Sn": ns > 0, "Zn": nz > 0}
    k = sum(has.values())

    if k == 1:
        return "Pure"
    if k == 4:
        return "FiveElem"       # Pt + Ga + In + Sn + Zn -> 5 elements total
    if k == 2:
        if has["Ga"] and has["Sn"]: return "GaSn"
        if has["Ga"] and has["In"]: return "GaIn"
        if has["In"] and has["Sn"]: return "InSn"
        if has["Zn"]:                return "PartZn"     # Ga-Zn, In-Zn, or Sn-Zn
    if k == 3:
        if has["Ga"] and has["In"] and has["Sn"]: return "GaInSn"
        if has["Zn"]:                               return "PartZn"  # any triple with Zn but not all 4
    raise ValueError(f"Unclassified composition: Ga={ng} In={ni} Sn={ns} Zn={nz}")


def main() -> None:
    df = pd.read_csv(SRC)
    if len(df) != 165:
        raise ValueError(f"Expected 165 rows, got {len(df)}")

    df["CatShort"] = [
        classify(int(r.Ga_count), int(r.In_count), int(r.Sn_count), int(r.Zn_count))
        for r in df.itertuples(index=False)
    ]

    # Integer-percent label to match reference ("0%", "12%", "25%", ...)
    df["GaPctLabel"] = df["Ga_pct"].apply(lambda v: f"{int(v)}%")
    bin_order = sorted(df["Ga_pct"].unique())   # 0, 12.5, ..., 100

    # Determine max rows needed (= largest bin population per single category)
    max_rows = 0
    for ga in bin_order:
        for cat in CAT_ORDER:
            n = len(df[(df["Ga_pct"] == ga) & (df["CatShort"] == cat)])
            max_rows = max(max_rows, n)

    rng = np.random.default_rng(SEED)
    all_cols: dict[str, list] = {}

    for ga in bin_order:
        label = f"{int(ga)}%"
        sub = df[df["Ga_pct"] == ga].copy()
        sub = sub.sort_values("ElementRef_Hf_kJ_mol").reset_index(drop=True)
        n_bin = len(sub)

        # Deterministic jitter matched 1:1 to sorted rows in this bin
        sub["X"] = np.round(rng.uniform(-JITTER_HALFWIDTH, JITTER_HALFWIDTH, size=n_bin), 3)

        for cat in CAT_ORDER:
            cat_rows = sub[sub["CatShort"] == cat]
            if len(cat_rows) == 0:
                continue
            x_col = f"{label}_{cat}_X"
            y_col = f"{label}_{cat}_Y"
            xs = cat_rows["X"].tolist() + [np.nan] * (max_rows - len(cat_rows))
            ys = cat_rows["ElementRef_Hf_kJ_mol"].round(4).tolist() + [np.nan] * (max_rows - len(cat_rows))
            all_cols[x_col] = xs
            all_cols[y_col] = ys

        # Mean + std for the bin
        mean_y = round(float(sub["ElementRef_Hf_kJ_mol"].mean()), 4)
        std_y  = round(float(sub["ElementRef_Hf_kJ_mol"].std(ddof=0)), 4) if n_bin > 1 else ""

        all_cols[f"{label}_Mean_X"]   = [0.0]    + [np.nan] * (max_rows - 1)
        all_cols[f"{label}_Mean_Y"]   = [mean_y] + [np.nan] * (max_rows - 1)
        all_cols[f"{label}_Mean_Std"] = [std_y]  + [np.nan] * (max_rows - 1)

    out = pd.DataFrame(all_cols)
    out.to_csv(OUT, index=False)
    print(f"[Panel g Origin] wrote {OUT}")
    print(f"  rows={len(out)}, cols={out.shape[1]}, bins={len(bin_order)}")

    # Sanity: count categories present per bin for the doc
    print("\n-- Categories present per Ga% bin --")
    for ga in bin_order:
        present = [c for c in CAT_ORDER
                   if len(df[(df["Ga_pct"] == ga) & (df["CatShort"] == c)]) > 0]
        n = len(df[df["Ga_pct"] == ga])
        print(f"  {int(ga):>3d}%: N={n:2d}, cats={present}")

    print("\n-- Per-bin mean Hf (kJ/mol) --")
    for ga in bin_order:
        sub = df[df["Ga_pct"] == ga]
        mean = float(sub["ElementRef_Hf_kJ_mol"].mean())
        std  = float(sub["ElementRef_Hf_kJ_mol"].std(ddof=0)) if len(sub) > 1 else 0.0
        print(f"  {int(ga):>3d}%: N={len(sub):2d}  mean={mean:+7.3f}  std={std:5.3f}")


if __name__ == "__main__":
    main()
