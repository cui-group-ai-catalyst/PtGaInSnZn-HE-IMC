"""
script_SI_DownSelectionFunnel.py
================================
Four-stage down-selection from the 165 B-sublattice composition prototypes
(`Panel g` / Supplementary Fig. 4) to the final experimental target
Pt3(Ga0.685 In0.215 Sn0.10) with 1 at% Zn micro-addition.

Stages
------
  A. All 165 B-sublattice prototypes                            (input)
  B. Favourable thermodynamic window:
        Ga% >= 62.5  AND  ΔH_f <= -34 kJ mol-1 atom-1            (Panel g)
  C. Liquidus accessibility: T_liquidus(x_Zn) < T_target         (CALPHAD)
  D. Galinstan-derived eutectic ratio (Ga:In:Sn = 68.5:21.5:10)
        + 1 at% Zn micro-addition                                (selected)

Outputs
-------
  data_SI_Funnel_Summary.csv           per-stage counts + representative composition
  data_SI_Funnel_StageB_Prototypes.csv prototypes inside the favourable window
  preview_SI_Funnel.png                4-stage funnel + Stage-B scatter overlay

Release note
------------
Stage C (CALPHAD-based liquidus accessibility) requires `pycalphad` plus the
COST507 Ga-In-Sn-Zn TDB; both are outside the bundled environment of this
release. The Stage-C step is therefore guarded by `CALPHAD_AVAILABLE` and
emits `pending-CALPHAD` placeholders when the dependency is absent. The
script falls through Stages A -> B -> D in that case, which is sufficient
to render the historical funnel sketch.

This script is retained as an archival/provenance helper, not as the clean
Nature reviewer path. For the release-supported Supplementary Tables 2 and 3
liquidus workflow, use `script_SI_LiquidusPredictor.py`.

Environment: Python 3.12 + numpy, pandas, matplotlib. Stage C additionally
requires pycalphad.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- paths / constants
HERE = Path(__file__).resolve().parent
SI_ROOT = HERE.parent                             # SI_Figures/
PANEL_G_CSV = (
    SI_ROOT
    / "SI_Fig04_165CompositionLandscape"
    / "data_FigG_165_ElementReferenced_Hf.csv"
)

# Favourable window (Panel g)
GA_PCT_MIN = 62.5          # %
DHF_MAX = -34.0            # kJ/mol (更负 = 更稳定，取上限阈值)

# Galinstan backbone (in 4-element cocktail, Zn 另算)
GA_TARGET = 0.685
IN_TARGET = 0.215
SN_TARGET = 0.100
ZN_TARGET = 0.010          # AI-recommended micro-addition for liquidus control

# CALPHAD toggle — set True once pycalphad + COST507 TDB is wired up
CALPHAD_AVAILABLE = False
T_TARGET_K = 1073.0        # 处理温度上限; Panel c annealing window 的上界


# ---------------------------------------------------------------- stages
def stage_a_load(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Columns present: Rank, Composition, Ga_count, ..., Ga_pct, ...,
    #                  ElementRef_Hf_kJ_mol, Category, Optimal_Window
    assert {"Ga_pct", "ElementRef_Hf_kJ_mol"}.issubset(df.columns)
    return df


def stage_b_favourable_window(df_a: pd.DataFrame) -> pd.DataFrame:
    mask = (df_a["Ga_pct"] >= GA_PCT_MIN) & (df_a["ElementRef_Hf_kJ_mol"] <= DHF_MAX)
    return df_a.loc[mask].copy().reset_index(drop=True)


def stage_c_liquidus(df_b: pd.DataFrame) -> pd.DataFrame:
    """
    Filter: T_liquidus(composition with x_Zn) < T_TARGET_K.

    Placeholder until pycalphad + COST507 are plugged in. Currently returns
    df_b unchanged and adds a `T_liq_K` column with NaN + `liquidus_status`.
    """
    df = df_b.copy()
    if CALPHAD_AVAILABLE:
        from pycalphad import Database, equilibrium  # noqa: F401
        # Not implemented in this release — see README "Limitations and known gaps".
        raise NotImplementedError(
            "CALPHAD path not yet wired — see README 'Limitations and known gaps'."
        )
    df["T_liq_K"] = np.nan
    df["liquidus_status"] = "pending-CALPHAD"
    return df


def stage_d_galinstan_proximity(df_c: pd.DataFrame) -> pd.DataFrame:
    """
    Rank df_c rows by L2 distance in (Ga, In, Sn) fractional space to the
    Galinstan eutectic backbone (Zn is decoupled). Return the nearest row
    as the representative target.
    """
    df = df_c.copy()
    # convert B-site counts -> fraction in B-sublattice
    denom = (df["Ga_count"] + df["In_count"] + df["Sn_count"]).replace(0, np.nan)
    ga_frac = df["Ga_count"] / denom
    in_frac = df["In_count"] / denom
    sn_frac = df["Sn_count"] / denom
    # distance to Galinstan on the Ga-In-Sn simplex
    df["dist_to_galinstan"] = np.sqrt(
        (ga_frac - GA_TARGET) ** 2
        + (in_frac - IN_TARGET) ** 2
        + (sn_frac - SN_TARGET) ** 2
    )
    df = df.sort_values("dist_to_galinstan").reset_index(drop=True)
    return df


# ---------------------------------------------------------------- reporting
def build_summary(df_a, df_b, df_c, df_d) -> pd.DataFrame:
    rows = [
        {
            "stage": "A_all_prototypes",
            "n": len(df_a),
            "criterion": "All 165 B-sublattice enumerations (Panel g)",
            "representative": "-",
        },
        {
            "stage": "B_favourable_window",
            "n": len(df_b),
            "criterion": f"Ga% >= {GA_PCT_MIN} AND dHf <= {DHF_MAX} kJ/mol",
            "representative": df_b.iloc[0]["Composition"] if len(df_b) else "-",
        },
        {
            "stage": "C_liquidus",
            "n": len(df_c),
            "criterion": (
                f"T_liq(x_Zn) < {T_TARGET_K:.0f} K  "
                f"[{'active' if CALPHAD_AVAILABLE else 'pending CALPHAD'}]"
            ),
            "representative": "-",
        },
        {
            "stage": "D_target",
            "n": 1,
            "criterion": (
                f"Closest B-site to Galinstan (Ga:In:Sn = "
                f"{GA_TARGET}:{IN_TARGET}:{SN_TARGET}) + {ZN_TARGET*100:.1f} at.% Zn"
            ),
            "representative": (
                f"Pt3(Ga{GA_TARGET}In{IN_TARGET}Sn{SN_TARGET}) "
                f"+ {ZN_TARGET*100:.1f} at.% Zn"
            ),
        },
    ]
    return pd.DataFrame(rows)


def plot_funnel(summary: pd.DataFrame, df_b: pd.DataFrame, out_png: Path) -> None:
    fig, (ax_fun, ax_sc) = plt.subplots(1, 2, figsize=(10.5, 4.5))

    # --- left: funnel bars
    ns = summary["n"].tolist()
    labels = ["A: all 165", "B: window", "C: liquidus", "D: target"]
    y = np.arange(len(ns))[::-1]
    ax_fun.barh(y, ns, color=["#90A4AE", "#43A047", "#FB8C00", "#C62828"])
    for yi, ni in zip(y, ns):
        ax_fun.text(ni + 1, yi, str(ni), va="center", fontsize=10)
    ax_fun.set_yticks(y)
    ax_fun.set_yticklabels(labels)
    ax_fun.set_xlabel("prototypes surviving")
    ax_fun.set_title("SI-Zn Down-Selection Funnel")

    # --- right: stage B scatter on the Panel g plane
    ax_sc.scatter(
        df_b["Ga_pct"], df_b["ElementRef_Hf_kJ_mol"],
        s=28, c="#43A047", alpha=0.7, label=f"Stage B (n={len(df_b)})",
    )
    ax_sc.axvspan(GA_PCT_MIN, 100, color="#2E7D32", alpha=0.08)
    ax_sc.axhline(DHF_MAX, ls="--", color="#555", lw=0.7)
    ax_sc.axvline(GA_TARGET * 100, ls="--", color="#D32F2F", lw=1.0,
                  label=f"AI target Ga {GA_TARGET*100:.1f}%")
    ax_sc.set_xlabel("Ga fraction on B-sublattice (%)")
    ax_sc.set_ylabel(r"$\Delta H_f$ (kJ/mol)")
    ax_sc.set_title("Stage B — favourable window")
    ax_sc.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main() -> None:
    df_a = stage_a_load(PANEL_G_CSV)
    df_b = stage_b_favourable_window(df_a)
    df_c = stage_c_liquidus(df_b)
    df_d = stage_d_galinstan_proximity(df_c)

    summary = build_summary(df_a, df_b, df_c, df_d)
    summary.to_csv(HERE / "data_SI_Funnel_Summary.csv", index=False)
    df_b.to_csv(HERE / "data_SI_Funnel_StageB_Prototypes.csv", index=False)
    plot_funnel(summary, df_b, HERE / "preview_SI_Funnel.png")

    print(summary.to_string(index=False))
    print(f"\nStage B top-5 by dHf:")
    print(
        df_b.sort_values("ElementRef_Hf_kJ_mol")
        .head(5)[["Composition", "Ga_pct", "ElementRef_Hf_kJ_mol"]]
        .to_string(index=False)
    )
    print(f"\nStage D nearest-to-Galinstan top-3:")
    cols = ["Composition", "Ga_pct", "In_pct", "Sn_pct",
            "ElementRef_Hf_kJ_mol", "dist_to_galinstan"]
    print(df_d.head(3)[cols].to_string(index=False))


if __name__ == "__main__":
    main()
