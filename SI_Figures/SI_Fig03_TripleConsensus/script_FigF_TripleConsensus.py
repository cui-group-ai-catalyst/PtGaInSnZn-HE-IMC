"""
script_FigF_TripleConsensus.py
==============================
SI Fig 3 — three-method consensus on 12 M-Ga binary formation enthalpies:
  (1) Corrected Miedema (at the actual host fraction c_host)
  (2) Materials Project DFT ground truth (formation_energy_per_atom)
  (3) UMA-s-1p1 Fairchem single-point, element-referenced
Plus a CHGNet column for SI-level comparison.

Default mode (no flags): post-processing only — read the shipped
`data_FigF_TripleConsensus_Data.csv` (4-method consensus already computed)
and regenerate the bar-chart PNG plus Origin-ready wide CSV.

`--rerun-uma`: re-run the full UMA + Miedema pipeline from CIFs. Requires
fairchem-core, ASE, and the local CIF tree which is NOT included in this
release. This branch is preserved for archival reference only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_CSV    = SCRIPT_DIR / "data_FigF_TripleConsensus_Data.csv"
SUMMARY_CSV = SCRIPT_DIR / "data_FigF_TripleConsensus_Summary.csv"
OUT_ORIGIN  = SCRIPT_DIR / "data_FigF_TripleConsensus_Origin_regen.csv"
OUT_PLOT    = SCRIPT_DIR / "preview_FigF_TripleConsensus_regen.png"

THREE_ABS = [
    "MP_DFT_abs_Hf_kJ_mol_atom",
    "UMA_abs_Hf_kJ_mol_atom",
    "CHGNet_abs_Hf_kJ_mol_atom",
]


def load_rho_dict() -> dict:
    if not SUMMARY_CSV.exists():
        return {}
    sdf = pd.read_csv(SUMMARY_CSV)
    out: dict = {}
    for _, row in sdf.iterrows():
        try:
            out[row["Metric"]] = float(row["Value"])
        except (TypeError, ValueError):
            pass
    return out


def save_origin_wide(df: pd.DataFrame) -> None:
    df = df.copy()
    df["Bar_Mean_absHf_kJ_mol_atom"] = df[THREE_ABS].mean(axis=1).round(4)
    df["Bar_Std_absHf_kJ_mol_atom"]  = df[THREE_ABS].std(axis=1, ddof=1).round(4)
    df["Consensus_Rank_Mean"] = df[["MP_DFT_Rank", "UMA_Rank", "CHGNet_Rank"]].mean(axis=1).round(2)
    df = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)

    cols = [
        "Host", "Formula", "c_host", "Size_Pass",
        "MP_DFT_Hf_kJ_mol_atom",   "MP_DFT_Rank",
        "UMA_Hf_kJ_mol_atom",      "UMA_Rank",
        "CHGNet_Hf_kJ_mol_atom",   "CHGNet_Rank",
        "MP_DFT_abs_Hf_kJ_mol_atom",
        "UMA_abs_Hf_kJ_mol_atom",
        "CHGNet_abs_Hf_kJ_mol_atom",
        "Bar_Mean_absHf_kJ_mol_atom",
        "Bar_Std_absHf_kJ_mol_atom",
        "Consensus_Rank_Mean",
    ]
    df[cols].to_csv(OUT_ORIGIN, index=False, float_format="%.4f")


def plot_consensus(df: pd.DataFrame, rho: dict) -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })

    order = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)
    x = np.arange(len(order))

    bar_mean = order[THREE_ABS].mean(axis=1).to_numpy()
    bar_std  = order[THREE_ABS].std(axis=1, ddof=1).to_numpy()

    fig, ax = plt.subplots(figsize=(12.0, 5.8))

    colors  = ["#C62828" if h == "Pt" else "#455A64" for h in order["Host"]]
    hatches = ["" if sp else "///" for sp in order["Size_Pass"]]

    for xi, mean_i, std_i, color, hatch in zip(x, bar_mean, bar_std, colors, hatches):
        ax.bar(xi, mean_i, width=0.72, color=color, alpha=0.90,
               edgecolor="black", linewidth=0.8, hatch=hatch, zorder=2)
        ax.errorbar(xi, mean_i, yerr=std_i, fmt="none",
                    ecolor="black", elinewidth=1.1, capsize=4, capthick=1.1,
                    zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{h}\n{f}" for h, f in zip(order["Host"], order["Formula"])],
        fontsize=8.5,
    )
    ax.set_ylabel(r"$|\Delta H_f|$  (kJ mol$^{-1}$ atom$^{-1}$)",
                  fontsize=11, fontweight="bold")
    ax.set_title(
        "SI Fig 3 | Binary M-Ga formation enthalpy magnitude by host\n"
        "Bar = mean over {MP-DFT, UMA, CHGNet};  error bar = std across the three methods",
        fontsize=10.5, fontweight="bold", pad=10,
    )
    ax.grid(axis="y", linestyle=":", alpha=0.35)
    y_max = float((bar_mean + bar_std).max()) * 1.12
    ax.set_ylim(0, y_max)
    ax.set_xlim(-0.6, len(order) - 0.4)

    handles = [
        Patch(facecolor="#C62828", edgecolor="black", label="Pt (focus host)"),
        Patch(facecolor="#455A64", edgecolor="black", label="Other hosts"),
        Patch(facecolor="white",   edgecolor="black", hatch="///",
              label="Size_Pass = False"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, frameon=True)

    if rho:
        info_lines = [f"Spearman rho (size-pass only, N={int(df['Size_Pass'].sum())}):"]
        if "MP_DFT_vs_UMA_sizepass" in rho:
            info_lines.append(f"  MP-DFT  vs UMA    = {rho['MP_DFT_vs_UMA_sizepass']:+.3f}")
        if "MP_DFT_vs_CHGNet_all" in rho:
            info_lines.append(f"  MP-DFT  vs CHGNet = {rho['MP_DFT_vs_CHGNet_all']:+.3f} (all-15)")
        if "UMA_vs_CHGNet_all" in rho:
            info_lines.append(f"  UMA     vs CHGNet = {rho['UMA_vs_CHGNet_all']:+.3f} (all-15)")
        ax.text(0.012, 0.98, "\n".join(info_lines), transform=ax.transAxes, fontsize=8.5,
                color="#333333", va="top", ha="left",
                bbox={"facecolor": "white", "alpha": 0.92, "edgecolor": "#9E9E9E"})

    plt.tight_layout()
    fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main_postprocess() -> None:
    if not DATA_CSV.exists():
        raise FileNotFoundError(f"Missing shipped data CSV: {DATA_CSV}")

    df = pd.read_csv(DATA_CSV)
    rho = load_rho_dict()

    df["Consensus_Rank_Mean"] = (
        df[["MP_DFT_Rank", "UMA_Rank", "CHGNet_Rank"]].mean(axis=1).round(2)
    )
    df = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)

    save_origin_wide(df)
    print(f"[Panel f post] Origin wide -> {OUT_ORIGIN}")

    plot_consensus(df, rho)
    print(f"[Panel f post] preview PNG -> {OUT_PLOT}")


def main_rerun_uma() -> None:
    """Full UMA + Miedema rebuild from CIFs.

    This branch requires fairchem-core, ASE, and the local CIF tree which
    is not included in this release. Kept here for archival reference; if
    you need to re-derive the consensus numbers from scratch, restore the
    CIF inputs alongside this script and adapt the data paths below.
    """
    raise NotImplementedError(
        "UMA rerun branch is not bundled with this release. "
        "Restore the binary/element CIF tree and the original "
        "this script with --rerun-uma to rebuild "
        "data_FigF_TripleConsensus_Data.csv from scratch."
    )


if __name__ == "__main__":
    if "--rerun-uma" in sys.argv:
        main_rerun_uma()
    else:
        main_postprocess()
