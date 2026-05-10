"""
script_FigA_v2_family_ordered.py
================================
Panel a v2 — heatmap of binary Miedema mixing enthalpies, restricted to the
20-host family-grouped set.

This release reads the shipped Origin-format CSV
(`data_FigA_v2_FamilyOrdered_Origin.csv`, rows = liquid partners, columns =
hosts already in family-grouped order) and regenerates the long-format CSV
plus the heatmap PNG. Family colours and family ordering are encoded
inline so the script is fully self-contained inside the release tree.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SOURCE_DATA = SCRIPT_DIR / "data_FigA_v2_FamilyOrdered_Origin.csv"

LIQUID_PARTNERS = ["Ga", "In", "Sn", "Zn", "Bi", "Hg"]

# Family assignment (matches the shipped column order in the Origin CSV)
FAMILY_OF = {
    "Pt": "PGM", "Ir": "PGM", "Pd": "PGM", "Rh": "PGM", "Ru": "PGM", "Os": "PGM",
    "Ni": "3d_TM", "Co": "3d_TM", "Fe": "3d_TM", "Cr": "3d_TM",
    "Re": "Other_TM", "Au": "Other_TM", "Zr": "Other_TM",
    "W": "Other_TM", "Hf": "Other_TM", "Mo": "Other_TM",
    "Ce": "RE_Group3", "La": "RE_Group3", "Y": "RE_Group3", "Sc": "RE_Group3",
}
FAMILY_COLOR = {
    "PGM": "#C62828",
    "3d_TM": "#1565C0",
    "Other_TM": "#2E7D32",
    "RE_Group3": "#6A1B9A",
}
FAMILY_LABELS = {"PGM": "PGMs", "3d_TM": "3d transition metals",
                 "Other_TM": "Other transition metals", "RE_Group3": "Group 3 / Rare earths"}


def main() -> None:
    src = pd.read_csv(SOURCE_DATA)
    hosts = [c for c in src.columns if c != "Target_Element"]

    # Long format for matplotlib / seaborn / re-imports
    long_rows = []
    for host in hosts:
        family = FAMILY_OF[host]
        for p in LIQUID_PARTNERS:
            v = float(src.loc[src["Target_Element"] == p, host].iloc[0])
            long_rows.append({"Host": host, "Family": family, "Partner": p, "dH_kJmol": round(v, 4)})
    long_df = pd.DataFrame(long_rows)
    csv_long = RESULTS_DIR / "data_FigA_v2_FamilyOrdered_Long.csv"
    long_df.to_csv(csv_long, index=False)
    print(f"Saved long-format CSV: {csv_long}")

    # Heatmap PNG
    matrix = src.set_index("Target_Element")[hosts].values  # (6, 20)
    families_in_order = [FAMILY_OF[h] for h in hosts]

    fig, ax = plt.subplots(figsize=(11, 4.4))

    # Custom diverging colormap centred at 0
    cmap = LinearSegmentedColormap.from_list(
        "blue_white",
        ["#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#FFFFFF"],
        N=256,
    )
    vmin = float(np.nanmin(matrix))
    vmax = 0.0
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")

    # Tick labels
    ax.set_xticks(range(len(hosts)))
    ax.set_xticklabels(hosts, fontsize=10)
    for tick, host in zip(ax.get_xticklabels(), hosts):
        tick.set_color(FAMILY_COLOR[FAMILY_OF[host]])
    ax.set_yticks(range(len(LIQUID_PARTNERS)))
    ax.set_yticklabels(LIQUID_PARTNERS, fontsize=10)

    # Family separators (vertical lines)
    boundaries = []
    last_fam = families_in_order[0]
    for i, fam in enumerate(families_in_order):
        if fam != last_fam:
            boundaries.append(i)
            last_fam = fam
    for b in boundaries:
        ax.axvline(b - 0.5, color="black", linewidth=1.2)

    # Cell annotations
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    fontsize=8, color="black" if v > -20 else "white")

    # Colour bar
    cbar = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.02)
    cbar.set_label(r"$\Delta H_{mix}$ (kJ mol$^{-1}$)", fontsize=10)

    # Family legend
    legend_handles = [Patch(facecolor=FAMILY_COLOR[f], label=FAMILY_LABELS[f])
                      for f in ["PGM", "3d_TM", "Other_TM", "RE_Group3"]]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=4, frameon=False, fontsize=9)

    ax.set_xlabel("")
    ax.set_ylabel("Liquid-metal partner", fontsize=10)
    ax.set_title("Panel a v2 — Binary mixing enthalpy of 20 hosts against six liquid metals "
                 "(family-grouped)", fontsize=11, pad=12)

    plt.tight_layout()
    out_png = RESULTS_DIR / "preview_FigA_heatmap.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved heatmap PNG: {out_png}")


if __name__ == "__main__":
    main()
