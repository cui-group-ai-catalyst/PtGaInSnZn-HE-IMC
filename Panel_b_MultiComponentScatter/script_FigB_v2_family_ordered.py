"""
script_FigB_v2_family_ordered.py
================================
Panel b v2 — multi-component (cocktail) only, 20 hosts in family-grouped order.

This release reads the shipped long-format CSV `data_FigB_v2_Long.csv`
(Host, Family, c_host, X_plot, dH_multi_kJmol, dG_multi_kJmol per scatter
point; 40 points per host, 800 rows total) and regenerates the preview
PNG plus the Origin-friendly wide CSVs. The thermodynamic core that
originally generated `data_FigB_v2_Long.csv` from Miedema parameters lives
outside this release; the values are bundled here directly.

Outputs (all in SCRIPT_DIR):
- preview_FigB_scatter.png
- data_FigB_v2_Origin_dH.csv (wide ({Host}_X, {Host}_Y) for ΔH)
- data_FigB_v2_Origin_dG.csv (same for ΔG)
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR

LONG_CSV = SCRIPT_DIR / "data_FigB_v2_Long.csv"

FAMILY_COLOR = {
    "PGM": "#C62828",
    "3d_TM": "#1565C0",
    "Other_TM": "#2E7D32",
    "RE_Group3": "#6A1B9A",
}
FAMILY_LABELS = {"PGM": "PGMs", "3d_TM": "3d TM", "Other_TM": "Other TM", "RE_Group3": "Group 3 / RE"}


def main() -> None:
    long_df = pd.read_csv(LONG_CSV)
    # Preserve the family-grouped order as it appears in the shipped CSV
    hosts_df = (long_df[["Host", "Host_Index", "Family"]]
                .drop_duplicates()
                .sort_values("Host_Index")
                .reset_index(drop=True))
    hosts_df["FamilyColor"] = hosts_df["Family"].map(FAMILY_COLOR)

    # Origin wide format: {Host}_X, {Host}_Y per host, separately for ΔH and ΔG
    def to_origin_wide(value_col: str) -> pd.DataFrame:
        cols = {}
        for host in hosts_df["Host"]:
            sub = long_df[long_df["Host"] == host]
            cols[f"{host}_X"] = sub["X_plot"].values
            cols[f"{host}_Y"] = sub[value_col].values
        return pd.DataFrame(cols)

    wide_dh = to_origin_wide("dH_multi_kJmol")
    wide_dg = to_origin_wide("dG_multi_kJmol")
    csv_wide_dh = RESULTS_DIR / "data_FigB_v2_Origin_dH.csv"
    csv_wide_dg = RESULTS_DIR / "data_FigB_v2_Origin_dG.csv"
    wide_dh.to_csv(csv_wide_dh, index=False)
    wide_dg.to_csv(csv_wide_dg, index=False)
    print(f"Saved Origin wide ΔH: {csv_wide_dh}")
    print(f"Saved Origin wide ΔG: {csv_wide_dg}")

    # Preview PNG: two stacked subplots (ΔH on top, ΔG on bottom)
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    for ax, ycol, ylabel in [
        (axes[0], "dH_multi_kJmol",
            r"Multi-component mixing enthalpy, $\Delta H_{mix}$ (kJ mol$^{-1}$)"),
        (axes[1], "dG_multi_kJmol",
            r"Multi-component mixing free energy at 500 K, $\Delta G_{mix}$ (kJ mol$^{-1}$)"),
    ]:
        for _, row in hosts_df.iterrows():
            host = row["Host"]
            color = row["FamilyColor"]
            sub = long_df[long_df["Host"] == host]
            size = 55 if host == "Pt" else 22
            edge = "black" if host == "Pt" else "none"
            ax.scatter(sub["X_plot"], sub[ycol], c=color, s=size, alpha=0.65,
                       edgecolors=edge, linewidths=0.6, zorder=3 if host == "Pt" else 2)

        ax.axhline(0, color="grey", linewidth=0.5, linestyle="--", zorder=1)

        # Family separators
        boundaries, last_fam = [], hosts_df.iloc[0]["Family"]
        for i, fam in enumerate(hosts_df["Family"]):
            if fam != last_fam:
                boundaries.append(i)
                last_fam = fam
        for b in boundaries:
            ax.axvline(b - 0.5, color="black", linewidth=0.8, linestyle=":", alpha=0.5)

        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(axis="y", alpha=0.25, linestyle=":")

    axes[1].set_xticks(range(len(hosts_df)))
    axes[1].set_xticklabels(hosts_df["Host"], fontsize=10)
    for tick, color in zip(axes[1].get_xticklabels(), hosts_df["FamilyColor"]):
        tick.set_color(color)

    legend_handles = [Patch(facecolor=FAMILY_COLOR[f], label=FAMILY_LABELS[f])
                      for f in ["PGM", "3d_TM", "Other_TM", "RE_Group3"]]
    axes[0].legend(handles=legend_handles, loc="lower right", ncol=4, frameon=True, fontsize=9)

    fig.suptitle("Panel b v2 — Multi-component cocktail thermodynamics across 20 hosts "
                 "(family-grouped, 40 perturbation points each)",
                 fontsize=12, y=0.995)
    plt.tight_layout()
    out_png = RESULTS_DIR / "preview_FigB_scatter.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved preview PNG: {out_png}")


if __name__ == "__main__":
    main()
