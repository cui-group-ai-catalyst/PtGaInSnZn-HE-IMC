"""
20260418_FigF_TripleConsensus.py
================================
Panel f 修订版 (2026-04-18) --- 三方 rank consensus:
  (1) Corrected Miedema (fixed at each binary's TRUE host fraction c_host)
  (2) Materials Project DFT ground truth (formation_energy_per_atom)
  (3) UMA-s-1p1 Fairchem single-point, element-referenced
External bonus: CHGNet values reloaded from 20260415_FigF_Validation.csv
as a 4th column for SI-level comparison.

Fix vs 20260415_FigF:
  - Miedema used c=0.5 equimolar regardless of real structure. Here we
    compute c_host = N_host_atoms / N_total_atoms from the localized CIF
    composition and use the regular-solution scaling dH_mix = 4 dH_AB c (1-c)
    so it is on the same per-mol-atom basis as DFT / UMA / CHGNet.
  - CHGNet-only ML cross-validation is upgraded to Miedema + MP-DFT + UMA.
    UMA column anchors the main figures; MP-DFT is the independent DFT ground
    truth; Miedema is the independent empirical method. CHGNet stays in the
    SI column.

Data sources
------------
- 02_data/20260415_FigF_LocalStructures/binary_manifest.csv
- 02_data/20260415_FigF_LocalStructures/element_manifest.csv
- 02_data/20260415_FigF_LocalStructures/binaries/*.cif
- 02_data/20260415_FigF_LocalStructures/elements/*.cif
- 03_results/20260418_FigF_MP_DFT_References.csv   (pulled once via MP API)
- 03_results/20260415_FigF_Validation.csv          (reloaded CHGNet numbers)

Outputs
-------
- 03_results/20260418_FigF_TripleConsensus_Data.csv   long form, per-host
- 03_results/20260418_FigF_TripleConsensus_Origin.csv Origin-ready wide
- 03_results/20260418_FigF_TripleConsensus_Summary.csv Spearman rho + key stats
- 03_results/20260418_FigF_TripleConsensus_Plot.png   preview figure
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.io import read
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "inputs" / "structures"
RESULTS_DIR = SCRIPT_DIR

BINARY_MANIFEST  = DATA_DIR / "binary_manifest.csv"
ELEMENT_MANIFEST = DATA_DIR / "element_manifest.csv"
MP_DFT_CSV       = RESULTS_DIR / "data_FigF_MP_DFT_References.csv"
CHGNET_CSV       = RESULTS_DIR / "data_FigF_CHGNet_References.csv"
CANONICAL_DATA   = RESULTS_DIR / "data_FigF_TripleConsensus_Data.csv"

OUT_DATA    = RESULTS_DIR / "data_FigF_TripleConsensus_Data_uma_regen.csv"
OUT_ORIGIN  = RESULTS_DIR / "data_FigF_TripleConsensus_Origin_uma_regen.csv"
OUT_SUMMARY = RESULTS_DIR / "data_FigF_TripleConsensus_Summary_uma_regen.csv"
OUT_PLOT    = RESULTS_DIR / "preview_FigF_TripleConsensus_uma_regen.png"
OUT_REPORT  = RESULTS_DIR / "validation_FigF_UMA_Rerun.json"

HOSTS = [
    "Pt", "Ir", "Pd", "Rh", "Ru", "Re", "Os", "Ni",
    "Au", "Co", "Ce", "La", "Fe", "Y", "Cu",
]

# Miedema Phi/n_ws — same set as 20260415_FigF
PARAMS = {
    "Ga": {"Phi": 4.10, "n_ws": 1.34},
    "Pt": {"Phi": 5.65, "n_ws": 1.78},
    "Ir": {"Phi": 5.55, "n_ws": 1.83},
    "Pd": {"Phi": 5.45, "n_ws": 1.67},
    "Rh": {"Phi": 5.40, "n_ws": 1.76},
    "Ru": {"Phi": 5.40, "n_ws": 1.83},
    "Re": {"Phi": 5.40, "n_ws": 1.96},
    "Os": {"Phi": 5.40, "n_ws": 1.98},
    "Ni": {"Phi": 5.20, "n_ws": 1.75},
    "Au": {"Phi": 5.15, "n_ws": 1.57},
    "Co": {"Phi": 5.10, "n_ws": 1.75},
    "Ce": {"Phi": 3.05, "n_ws": 1.13},
    "La": {"Phi": 3.05, "n_ws": 1.09},
    "Fe": {"Phi": 4.93, "n_ws": 1.77},
    "Y":  {"Phi": 3.20, "n_ws": 1.11},
    "Cu": {"Phi": 4.45, "n_ws": 1.47},
}
P_MIED, Q_MIED = 14.1, 9.4
EV_TO_KJ = 96.485


# -----------------------------------------------------------------------------
# UMA calculator
# -----------------------------------------------------------------------------
def load_uma_calc() -> FAIRChemCalculator:
    default_checkpoint = (
        Path(os.environ["USERPROFILE"])
        / ".cache" / "fairchem"
        / "models--facebook--UMA"
        / "snapshots"
        / "38529caa2c51a9a8a0d71f0b56b79ac33bc9eceb"
        / "checkpoints" / "uma-s-1p1.pt"
    )
    ckpt = Path(os.environ.get("UMA_CHECKPOINT", default_checkpoint)).expanduser()
    if not ckpt.exists():
        raise FileNotFoundError(f"UMA checkpoint not found at {ckpt}")
    return FAIRChemCalculator.from_model_checkpoint(
        str(ckpt), task_name="oc20", device=os.environ.get("UMA_DEVICE", "cpu")
    )


def manifest_path(value: str) -> Path:
    """Resolve historical Windows-style manifest paths on every platform."""
    return DATA_DIR.joinpath(*str(value).replace("\\", "/").split("/"))


def uma_energy_per_atom(calc: FAIRChemCalculator, cif_path: Path) -> tuple[float, int]:
    atoms = read(cif_path)
    atoms.calc = calc
    e_tot = float(atoms.get_potential_energy())
    n = int(len(atoms))
    return e_tot / n, n


# -----------------------------------------------------------------------------
# Miedema (fixed at real c_host)
# -----------------------------------------------------------------------------
def miedema_dH_AB(host: str, solvent: str = "Ga") -> float:
    """Binary pair interaction energy at c=0.5 (kJ/mol-atom)."""
    d_phi = PARAMS[host]["Phi"] - PARAMS[solvent]["Phi"]
    d_n   = PARAMS[host]["n_ws"] ** (1/3) - PARAMS[solvent]["n_ws"] ** (1/3)
    return -P_MIED * d_phi ** 2 + Q_MIED * d_n ** 2


def miedema_dH_at_composition(host: str, c_host: float) -> float:
    """Regular-solution scaled mixing enthalpy at the real host fraction."""
    dH_AB = miedema_dH_AB(host)
    return 4.0 * dH_AB * c_host * (1.0 - c_host)


def c_host_from_cif(cif_path: Path, host: str) -> float:
    """Read atom fractions from the binary CIF, return c_host = n_host/n_total."""
    atoms = read(cif_path)
    syms = atoms.get_chemical_symbols()
    n_host = sum(1 for s in syms if s == host)
    return n_host / len(syms)


# -----------------------------------------------------------------------------
# Main consensus build
# -----------------------------------------------------------------------------
def build_consensus_table() -> tuple[pd.DataFrame, dict]:
    binary_df  = pd.read_csv(BINARY_MANIFEST)
    element_df = pd.read_csv(ELEMENT_MANIFEST)
    mp_df      = pd.read_csv(MP_DFT_CSV)
    chgnet_df  = pd.read_csv(CHGNET_CSV)

    binary_df  = binary_df.set_index("Host")
    element_df = element_df.set_index("Element")
    mp_df      = mp_df.set_index("Host")
    chgnet_df  = chgnet_df.set_index("Host")

    missing = sorted(set(HOSTS) - set(binary_df.index))
    if missing:
        raise ValueError(f"Missing hosts from binary_manifest: {missing}")

    print("[UMA] Loading checkpoint (CPU)...")
    calc = load_uma_calc()

    # UMA element references (15 hosts + Ga)
    uma_elem = {}
    needed_elements = sorted(set(HOSTS) | {"Ga"})
    print(f"[UMA] Scoring {len(needed_elements)} element references:")
    for el in needed_elements:
        row = element_df.loc[el]
        e_pa, n = uma_energy_per_atom(calc, manifest_path(row["CIF_File"]))
        uma_elem[el] = e_pa
        print(f"   {el:3s}  N={n:3d}  E={e_pa:+.6f} eV/atom")

    # Miedema + UMA binary + consolidate
    records = []
    print(f"\n[UMA] Scoring {len(HOSTS)} binary CIFs:")
    for host in HOSTS:
        bin_row = binary_df.loc[host]
        cif = manifest_path(bin_row["CIF_File"])
        e_alloy_pa, n_atoms = uma_energy_per_atom(calc, cif)

        c = c_host_from_cif(cif, host)
        e_ref_pa = c * uma_elem[host] + (1.0 - c) * uma_elem["Ga"]
        uma_dHf_kJ = (e_alloy_pa - e_ref_pa) * EV_TO_KJ

        mied_dHmix = miedema_dH_at_composition(host, c)

        mp_row = mp_df.loc[host]
        chg_row = chgnet_df.loc[host]

        records.append({
            "Host":              host,
            "Formula":           bin_row["Formula"],
            "MP_ID":             bin_row["MP_ID"],
            "Binary_N_sites":    int(bin_row["N_sites"]),
            "c_host":            round(c, 6),
            "Size_Pass":         bool(bin_row["Size_Pass"]),

            # Miedema (fixed at real c)
            "Miedema_dH_kJ_mol_atom":           round(mied_dHmix, 4),
            "Miedema_abs_dH_kJ_mol_atom":       round(abs(mied_dHmix), 4),

            # MP DFT (ground truth)
            "MP_DFT_Hf_kJ_mol_atom":            round(float(mp_row["MP_DFT_Hf_kJ_mol_atom"]), 4),
            "MP_DFT_abs_Hf_kJ_mol_atom":        round(abs(float(mp_row["MP_DFT_Hf_kJ_mol_atom"])), 4),
            "MP_DFT_EAboveHull_eV_atom":        round(float(mp_row["MP_DFT_EAboveHull_eV_atom"]), 6),

            # UMA (main-line)
            "UMA_Alloy_E_eV_atom":              round(e_alloy_pa, 6),
            "UMA_ElementRef_E_eV_atom":         round(e_ref_pa, 6),
            "UMA_Hf_kJ_mol_atom":               round(uma_dHf_kJ, 4),
            "UMA_abs_Hf_kJ_mol_atom":           round(abs(uma_dHf_kJ), 4),

            # CHGNet (SI compare)
            "CHGNet_Hf_kJ_mol_atom":            round(float(chg_row["CHGNet_Hf_kJ_mol_atom"]), 4),
            "CHGNet_abs_Hf_kJ_mol_atom":        round(abs(float(chg_row["CHGNet_Hf_kJ_mol_atom"])), 4),
        })
        print(f"   {host:3s}  c={c:.4f}  UMA Hf={uma_dHf_kJ:+.3f}  "
              f"MP-DFT Hf={mp_row['MP_DFT_Hf_kJ_mol_atom']:+.3f}  "
              f"Miedema={mied_dHmix:+.3f}  CHGNet={chg_row['CHGNet_Hf_kJ_mol_atom']:+.3f}")

    df = pd.DataFrame(records)

    # Ranks: 1 = strongest binding (largest |dH|)
    for prefix in ["Miedema", "MP_DFT", "UMA", "CHGNet"]:
        key = f"{prefix}_abs_" + ("dH_kJ_mol_atom" if prefix == "Miedema" else "Hf_kJ_mol_atom")
        df[f"{prefix}_Rank"] = df[key].rank(ascending=False, method="min").astype(int)

    # Pairwise Spearman rho (rank correlation is the standard metric here)
    def spearman(a, b):
        return float(np.corrcoef(
            pd.Series(a).rank(), pd.Series(b).rank()
        )[0, 1])

    full = df
    size = df[df["Size_Pass"]]

    rho = {
        "Miedema_vs_MP_DFT_all":      spearman(full["Miedema_abs_dH_kJ_mol_atom"], full["MP_DFT_abs_Hf_kJ_mol_atom"]),
        "Miedema_vs_UMA_all":         spearman(full["Miedema_abs_dH_kJ_mol_atom"], full["UMA_abs_Hf_kJ_mol_atom"]),
        "MP_DFT_vs_UMA_all":          spearman(full["MP_DFT_abs_Hf_kJ_mol_atom"],  full["UMA_abs_Hf_kJ_mol_atom"]),
        "MP_DFT_vs_CHGNet_all":       spearman(full["MP_DFT_abs_Hf_kJ_mol_atom"],  full["CHGNet_abs_Hf_kJ_mol_atom"]),
        "UMA_vs_CHGNet_all":          spearman(full["UMA_abs_Hf_kJ_mol_atom"],     full["CHGNet_abs_Hf_kJ_mol_atom"]),
        "Miedema_vs_CHGNet_all":      spearman(full["Miedema_abs_dH_kJ_mol_atom"], full["CHGNet_abs_Hf_kJ_mol_atom"]),
        "Miedema_vs_MP_DFT_sizepass": spearman(size["Miedema_abs_dH_kJ_mol_atom"], size["MP_DFT_abs_Hf_kJ_mol_atom"]),
        "Miedema_vs_UMA_sizepass":    spearman(size["Miedema_abs_dH_kJ_mol_atom"], size["UMA_abs_Hf_kJ_mol_atom"]),
        "MP_DFT_vs_UMA_sizepass":     spearman(size["MP_DFT_abs_Hf_kJ_mol_atom"],  size["UMA_abs_Hf_kJ_mol_atom"]),
    }

    # Manuscript consensus = mean of the three atomistic methods.
    df["Consensus_Rank_Mean"] = df[["MP_DFT_Rank", "UMA_Rank", "CHGNet_Rank"]].mean(axis=1).round(2)
    df = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank", "MP_DFT_Rank"]).reset_index(drop=True)

    return df, rho


# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
def save_origin_wide(df: pd.DataFrame) -> None:
    """Origin-ready wide table: hosts as rows, per-method numeric columns plus
    the bar-chart quantities (mean/std of DFT+UMA+CHGNet |Hf|) and the
    Miedema triangle-marker value.
    """
    df = df.copy()
    three = ["MP_DFT_abs_Hf_kJ_mol_atom",
             "UMA_abs_Hf_kJ_mol_atom",
             "CHGNet_abs_Hf_kJ_mol_atom"]
    df["Bar_Mean_absHf_kJ_mol_atom"] = df[three].mean(axis=1).round(4)
    df["Bar_Std_absHf_kJ_mol_atom"]  = df[three].std(axis=1, ddof=1).round(4)
    df["Miedema_Marker_absdH_kJ_mol_atom"] = df["Miedema_abs_dH_kJ_mol_atom"].round(4)

    cols = [
        "Host", "Formula", "c_host", "Size_Pass",
        # signed per method
        "Miedema_dH_kJ_mol_atom",  "Miedema_Rank",
        "MP_DFT_Hf_kJ_mol_atom",   "MP_DFT_Rank",
        "UMA_Hf_kJ_mol_atom",      "UMA_Rank",
        "CHGNet_Hf_kJ_mol_atom",   "CHGNet_Rank",
        # absolute per method (bar chart inputs)
        "MP_DFT_abs_Hf_kJ_mol_atom",
        "UMA_abs_Hf_kJ_mol_atom",
        "CHGNet_abs_Hf_kJ_mol_atom",
        # bar + marker + consensus
        "Bar_Mean_absHf_kJ_mol_atom",
        "Bar_Std_absHf_kJ_mol_atom",
        "Miedema_Marker_absdH_kJ_mol_atom",
        "Consensus_Rank_Mean",
    ]
    df[cols].to_csv(OUT_ORIGIN, index=False, float_format="%.4f")


def save_summary(rho: dict, df: pd.DataFrame) -> None:
    rows = [{"Metric": k, "Value": round(v, 4)} for k, v in rho.items()]
    rows.append({"Metric": "N_hosts_all", "Value": len(df)})
    rows.append({"Metric": "N_hosts_sizepass", "Value": int(df["Size_Pass"].sum())})
    pt = df[df["Host"] == "Pt"].iloc[0]
    rows.extend([
        {"Metric": "Pt_Miedema_Rank",   "Value": int(pt["Miedema_Rank"])},
        {"Metric": "Pt_MP_DFT_Rank",    "Value": int(pt["MP_DFT_Rank"])},
        {"Metric": "Pt_UMA_Rank",       "Value": int(pt["UMA_Rank"])},
        {"Metric": "Pt_CHGNet_Rank",    "Value": int(pt["CHGNet_Rank"])},
        {"Metric": "Pt_Consensus_Rank", "Value": float(pt["Consensus_Rank_Mean"])},
    ])
    pd.DataFrame(rows).to_csv(OUT_SUMMARY, index=False)


def plot_consensus(df: pd.DataFrame, rho: dict) -> None:
    """Horizontal-hosts, vertical |Hf| bar chart:
      bar = mean(|MP-DFT|, |UMA|, |CHGNet|)
      err = std of the three
      triangle overlay per host = Miedema |dH| (independent empirical check)
      hosts sorted by consensus rank (strongest binder on the left)
      Pt bar red-highlighted; Size_Pass=False hosts hatched.
    """
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })

    order = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)
    x = np.arange(len(order))

    three = ["MP_DFT_abs_Hf_kJ_mol_atom",
             "UMA_abs_Hf_kJ_mol_atom",
             "CHGNet_abs_Hf_kJ_mol_atom"]
    bar_mean = order[three].mean(axis=1).to_numpy()
    bar_std  = order[three].std(axis=1, ddof=1).to_numpy()
    miedema  = order["Miedema_abs_dH_kJ_mol_atom"].to_numpy()

    fig, ax = plt.subplots(figsize=(12.5, 6.2))

    colors = ["#C62828" if h == "Pt" else "#455A64" for h in order["Host"]]
    hatches = ["" if sp else "///" for sp in order["Size_Pass"]]

    # Main bar: mean |Hf| over DFT + UMA + CHGNet, error bar = std
    for xi, mean_i, std_i, color, hatch in zip(x, bar_mean, bar_std, colors, hatches):
        ax.bar(xi, mean_i, width=0.72, color=color, alpha=0.88,
               edgecolor="black", linewidth=0.8, hatch=hatch, zorder=2)
        ax.errorbar(xi, mean_i, yerr=std_i, fmt="none",
                    ecolor="black", elinewidth=1.1, capsize=4, capthick=1.1,
                    zorder=4)

    # Per-method transparent markers on top of the bar (so the bar shows DFT+UMA+CHGNet spread)
    method_markers = {
        "MP_DFT_abs_Hf_kJ_mol_atom":  ("o", "#43A047", "MP-DFT"),
        "UMA_abs_Hf_kJ_mol_atom":     ("s", "#E65100", "UMA"),
        "CHGNet_abs_Hf_kJ_mol_atom":  ("D", "#8E24AA", "CHGNet"),
    }
    for col, (marker, mc, label) in method_markers.items():
        ax.scatter(x, order[col].to_numpy(),
                   marker=marker, s=38, facecolor="white",
                   edgecolor=mc, linewidths=1.4, zorder=5, label=label)

    # Miedema triangle overlay (independent empirical method)
    ax.scatter(x, miedema, marker="^", s=110,
               facecolor="#FDD835", edgecolor="black", linewidths=1.1,
               zorder=6, label="Miedema (empirical)")

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{h}\n{f}"
         for h, f in zip(order["Host"], order["Formula"])],
        fontsize=8.5,
    )
    ax.set_ylabel(r"$|\Delta H_f|$  (kJ mol$^{-1}$ atom$^{-1}$)",
                  fontsize=11, fontweight="bold")
    ax.set_title(
        "Panel f (2026-04-18) | Binary M-Ga formation enthalpy magnitude by host\n"
        "Bar = mean(|MP-DFT|, |UMA|, |CHGNet|), error = std; "
        "triangle = Miedema (independent empirical)",
        fontsize=10.5, fontweight="bold", pad=10,
    )
    ax.grid(axis="y", linestyle=":", alpha=0.35)
    y_max = max(bar_mean.max() + bar_std.max(), miedema.max()) * 1.12
    ax.set_ylim(0, y_max)
    ax.set_xlim(-0.6, len(order) - 0.4)

    # Legend: methods + size-fail swatch + Pt swatch
    from matplotlib.patches import Patch
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([
        Patch(facecolor="#C62828", edgecolor="black", label="Pt (focus host)"),
        Patch(facecolor="white", edgecolor="black", hatch="///",
              label="Size_Pass = False"),
    ])
    ax.legend(handles=handles, loc="upper right", fontsize=8.5,
              frameon=True, ncol=2)

    info = (
        f"Spearman rho (size-pass only, N={int(df['Size_Pass'].sum())}):\n"
        f"  Miedema vs MP-DFT = {rho['Miedema_vs_MP_DFT_sizepass']:+.3f}\n"
        f"  Miedema vs UMA    = {rho['Miedema_vs_UMA_sizepass']:+.3f}\n"
        f"  MP-DFT  vs UMA    = {rho['MP_DFT_vs_UMA_sizepass']:+.3f}"
    )
    ax.text(0.012, 0.98, info, transform=ax.transAxes, fontsize=8.5,
            color="#333333", va="top", ha="left",
            bbox={"facecolor": "white", "alpha": 0.92, "edgecolor": "#9E9E9E"})

    plt.tight_layout()
    fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def validate_against_canonical(df: pd.DataFrame) -> dict:
    canonical = pd.read_csv(CANONICAL_DATA).set_index("Host")
    regenerated = df.set_index("Host")
    if set(canonical.index) != set(regenerated.index):
        raise AssertionError("Canonical and regenerated host sets differ")
    differences = {
        column: float((canonical[column] - regenerated[column]).abs().max())
        for column in (
            "UMA_Alloy_E_eV_atom",
            "UMA_ElementRef_E_eV_atom",
            "UMA_Hf_kJ_mol_atom",
        )
    }
    passed = (
        differences["UMA_Alloy_E_eV_atom"] <= 1.0e-5
        and differences["UMA_ElementRef_E_eV_atom"] <= 1.0e-5
        and differences["UMA_Hf_kJ_mol_atom"] <= 1.0e-3
    )
    report = {
        "status": "passed" if passed else "failed",
        "model": "UMA-s-1p1",
        "n_binary_structures": int(len(regenerated)),
        "n_element_reference_structures": 16,
        "max_abs_differences": differences,
        "tolerances": {
            "energy_eV_atom": 1.0e-5,
            "formation_enthalpy_kJ_mol_atom": 1.0e-3,
        },
    }
    OUT_REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise AssertionError("SI Fig. 3 UMA regeneration exceeded tolerance")
    return report


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if not MP_DFT_CSV.exists():
        raise FileNotFoundError(f"MP DFT references missing: {MP_DFT_CSV}")
    if not CHGNET_CSV.exists():
        raise FileNotFoundError(f"CHGNet validation CSV missing: {CHGNET_CSV}")

    df, rho = build_consensus_table()

    df.to_csv(OUT_DATA, index=False, float_format="%.6f")
    save_origin_wide(df)
    save_summary(rho, df)
    plot_consensus(df, rho)
    report = validate_against_canonical(df)

    print(f"\n[Panel f] data (per host, all methods) -> {OUT_DATA}")
    print(f"[Panel f] Origin-ready wide            -> {OUT_ORIGIN}")
    print(f"[Panel f] summary (Spearman + ranks)   -> {OUT_SUMMARY}")
    print(f"[Panel f] preview plot                 -> {OUT_PLOT}")
    print(f"[Panel f] validation report            -> {OUT_REPORT}")
    print(json.dumps(report, indent=2))

    print("\n-- Spearman rho (rank correlation, size-pass hosts only) --")
    print(f"  Miedema vs MP-DFT  = {rho['Miedema_vs_MP_DFT_sizepass']:+.4f}")
    print(f"  Miedema vs UMA     = {rho['Miedema_vs_UMA_sizepass']:+.4f}")
    print(f"  MP-DFT  vs UMA     = {rho['MP_DFT_vs_UMA_sizepass']:+.4f}")

    print("\n-- Pt ranks --")
    pt = df[df["Host"] == "Pt"].iloc[0]
    print(f"  Miedema  rank = {int(pt['Miedema_Rank'])}")
    print(f"  MP-DFT   rank = {int(pt['MP_DFT_Rank'])}")
    print(f"  UMA      rank = {int(pt['UMA_Rank'])}")
    print(f"  CHGNet   rank = {int(pt['CHGNet_Rank'])}")
    print(f"  3-way consensus rank (mean) = {pt['Consensus_Rank_Mean']:.2f}")


if __name__ == "__main__":
    main()
