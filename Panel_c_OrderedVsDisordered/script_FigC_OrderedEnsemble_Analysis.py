"""Compare 30 ordered B-sublattice motifs with 30 random full-lattice motifs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


SCRIPT_DIR = Path(__file__).resolve().parent
RELEASE_ROOT = SCRIPT_DIR.parent
DEFAULT_ORDERED = SCRIPT_DIR / "data_FigC_OrderedEnsemble_Raw_UMA_Energies_regen.csv"
DEFAULT_DISORDERED = SCRIPT_DIR / "data_FigC_Raw_UMA_Energies.csv"
DEFAULT_REFS = RELEASE_ROOT / "shared" / "UMA_Element_Reference_Energies.csv"
DEFAULT_LONG = SCRIPT_DIR / "data_FigC_OrderedVsDisordered_Ensembles_regen.csv"
DEFAULT_SUMMARY_CSV = SCRIPT_DIR / "data_FigC_OrderedVsDisordered_EnsembleSummary_regen.csv"
DEFAULT_SUMMARY_JSON = SCRIPT_DIR / "validation_FigC_OrderedVsDisordered_Ensembles_regen.json"
DEFAULT_FIGURE = SCRIPT_DIR / "preview_FigC_OrderedVsDisordered_Ensembles_regen"

EV_TO_KJMOL = 96.485
EXPECTED_COMPOSITION = {"Pt": 24, "Ga": 2, "In": 2, "Sn": 2, "Zn": 2}
HISTORICAL_ORDERED_E_EV_ATOM = -4.945261


def reference_energy(refs_path: Path) -> float:
    refs = pd.read_csv(refs_path).set_index("Element")["UMA_E_eV_atom"].to_dict()
    missing = set(EXPECTED_COMPOSITION) - set(refs)
    if missing:
        raise ValueError(f"Missing element references: {sorted(missing)}")
    return sum(EXPECTED_COMPOSITION[element] * refs[element] for element in refs) / 32.0


def load_data(ordered_path: Path, disordered_path: Path, refs_path: Path) -> pd.DataFrame:
    ordered = pd.read_csv(ordered_path)
    disordered_raw = pd.read_csv(disordered_path)
    disordered = disordered_raw[disordered_raw["Type"] == "Disordered_Random"].copy()
    if len(ordered) != 30 or len(disordered) != 30:
        raise ValueError(f"Expected 30 ordered and 30 disordered rows, got {len(ordered)}, {len(disordered)}")
    if ordered["Canonical_Class"].nunique() != 30:
        raise ValueError("Ordered rows are not from 30 distinct symmetry classes")
    if ordered["Assignment_SHA256"].nunique() != 30:
        raise ValueError("Ordered rows contain duplicate occupancy fingerprints")

    reference = reference_energy(refs_path)
    ordered_long = pd.DataFrame(
        {
            "State": "Ordered_B_sublattice",
            "Config_ID": ordered["Config_ID"].astype(int),
            "Seed": ordered["Selection_Seed"],
            "Assignment": ordered["Assignment_B0_to_B7"],
            "Class_Degeneracy": ordered["Class_Degeneracy"].astype(int),
            "Energy_eV_atom": ordered["Energy_eV_atom"].astype(float),
        }
    )
    disordered_long = pd.DataFrame(
        {
            "State": "Disordered_full_lattice",
            "Config_ID": disordered["Config_ID"].astype(int),
            "Seed": disordered["Seed"].astype(int),
            "Assignment": "",
            "Class_Degeneracy": np.nan,
            "Energy_eV_atom": disordered["Energy_eV_atom"].astype(float),
        }
    )
    combined = pd.concat([ordered_long, disordered_long], ignore_index=True)
    combined["ElementRef_Hf_kJ_mol_atom"] = (
        combined["Energy_eV_atom"] - reference
    ) * EV_TO_KJMOL
    return combined


def summarize(combined: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    ordered = combined.loc[
        combined["State"] == "Ordered_B_sublattice", "ElementRef_Hf_kJ_mol_atom"
    ].to_numpy(dtype=float)
    disordered = combined.loc[
        combined["State"] == "Disordered_full_lattice", "ElementRef_Hf_kJ_mol_atom"
    ].to_numpy(dtype=float)

    def group_stats(values):
        n = len(values)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        sem = sd / math.sqrt(n)
        ci = stats.t.interval(0.95, n - 1, loc=mean, scale=sem)
        return {
            "n": n,
            "mean": mean,
            "sd": sd,
            "sem": sem,
            "ci95_low": float(ci[0]),
            "ci95_high": float(ci[1]),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
        }

    ordered_stats = group_stats(ordered)
    disordered_stats = group_stats(disordered)
    gap = disordered_stats["mean"] - ordered_stats["mean"]
    se_gap = math.sqrt(
        ordered_stats["sd"] ** 2 / ordered_stats["n"]
        + disordered_stats["sd"] ** 2 / disordered_stats["n"]
    )
    numerator = (
        ordered_stats["sd"] ** 2 / ordered_stats["n"]
        + disordered_stats["sd"] ** 2 / disordered_stats["n"]
    ) ** 2
    denominator = (
        (ordered_stats["sd"] ** 2 / ordered_stats["n"]) ** 2
        / (ordered_stats["n"] - 1)
        + (disordered_stats["sd"] ** 2 / disordered_stats["n"]) ** 2
        / (disordered_stats["n"] - 1)
    )
    welch_df = numerator / denominator
    t_critical = float(stats.t.ppf(0.975, welch_df))
    gap_ci = (gap - t_critical * se_gap, gap + t_critical * se_gap)
    pooled_sd = math.sqrt(
        (
            (ordered_stats["n"] - 1) * ordered_stats["sd"] ** 2
            + (disordered_stats["n"] - 1) * disordered_stats["sd"] ** 2
        )
        / (ordered_stats["n"] + disordered_stats["n"] - 2)
    )
    nonoverlap_margin = disordered_stats["min"] - ordered_stats["max"]

    summary_rows = []
    for state, values in (
        ("Ordered_B_sublattice", ordered_stats),
        ("Disordered_full_lattice", disordered_stats),
    ):
        summary_rows.append({"State": state, **values})
    summary_df = pd.DataFrame(summary_rows)
    report = {
        "status": "passed",
        "provenance_mode": "measured_computational",
        "energy_protocol": "fixed-cell fixed-coordinate UMA-s-1p1 single-point",
        "ordered": ordered_stats,
        "disordered": disordered_stats,
        "mean_gap_disordered_minus_ordered_kJ_mol_atom": gap,
        "gap_sem_kJ_mol_atom": se_gap,
        "gap_welch_df": welch_df,
        "gap_ci95_kJ_mol_atom": [float(gap_ci[0]), float(gap_ci[1])],
        "nonoverlap_margin_kJ_mol_atom": nonoverlap_margin,
        "cohen_d_pooled": gap / pooled_sd,
        "historical_anchor_abs_difference_eV_atom": float(
            abs(
                combined.loc[
                    (combined["State"] == "Ordered_B_sublattice")
                    & (combined["Config_ID"] == 0),
                    "Energy_eV_atom",
                ].iloc[0]
                - HISTORICAL_ORDERED_E_EV_ATOM
            )
        ),
    }
    return summary_df, report


def create_figure(combined: pd.DataFrame, report: dict, output_base: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    palette = ["#0072B2", "#D55E00"]
    labels = ["Ordered\nB-sublattice", "Disordered\nfull lattice"]
    states = ["Ordered_B_sublattice", "Disordered_full_lattice"]
    values = [
        combined.loc[combined["State"] == state, "ElementRef_Hf_kJ_mol_atom"].to_numpy()
        for state in states
    ]

    fig, ax = plt.subplots(figsize=(4.45, 3.35), constrained_layout=True)
    violin = ax.violinplot(values, positions=[0, 1], widths=0.64, showextrema=False)
    for body, color in zip(violin["bodies"], palette):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.16)

    rng = np.random.default_rng(20260727)
    for position, (group, color) in enumerate(zip(values, palette)):
        jitter = rng.uniform(-0.16, 0.16, size=len(group))
        ax.scatter(
            np.full(len(group), position) + jitter,
            group,
            s=22,
            facecolor=color,
            edgecolor="white",
            linewidth=0.45,
            alpha=0.88,
            zorder=3,
        )
        mean = float(np.mean(group))
        ci = stats.t.interval(0.95, len(group) - 1, loc=mean, scale=stats.sem(group))
        ax.errorbar(
            position,
            mean,
            yerr=[[mean - ci[0]], [ci[1] - mean]],
            fmt="D",
            markersize=5.2,
            color="#202020",
            markerfacecolor="white",
            markeredgewidth=1.0,
            capsize=3,
            linewidth=1.1,
            zorder=5,
        )

    gap = report["mean_gap_disordered_minus_ordered_kJ_mol_atom"]
    ci_low, ci_high = report["gap_ci95_kJ_mol_atom"]
    upper = max(max(group) for group in values) + 1.8
    ax.plot([0, 0, 1, 1], [upper - 0.35, upper, upper, upper - 0.35], color="#303030", lw=0.8)
    ax.text(
        0.5,
        upper + 0.25,
        f"mean gap = {gap:.2f} kJ mol$^{{-1}}$ atom$^{{-1}}$\n95% CI: {ci_low:.2f} to {ci_high:.2f}",
        ha="center",
        va="bottom",
        fontsize=8.2,
    )

    ax.set_xticks([0, 1], labels)
    ax.set_ylabel(r"Element-referenced $\Delta H_f$ (kJ mol$^{-1}$ atom$^{-1}$)")
    ax.set_xlim(-0.48, 1.48)
    ax.set_ylim(min(min(group) for group in values) - 2.2, upper + 4.0)
    ax.axhline(0, color="#808080", linewidth=0.6, linestyle="--", alpha=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="in", length=3)
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 600}),
        (".tiff", {"dpi": 600}),
    ):
        fig.savefig(output_base.with_suffix(suffix), bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ordered", type=Path, default=DEFAULT_ORDERED)
    parser.add_argument("--disordered", type=Path, default=DEFAULT_DISORDERED)
    parser.add_argument("--refs", type=Path, default=DEFAULT_REFS)
    parser.add_argument("--long-output", type=Path, default=DEFAULT_LONG)
    parser.add_argument("--summary-csv", type=Path, default=DEFAULT_SUMMARY_CSV)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY_JSON)
    parser.add_argument("--figure-base", type=Path, default=DEFAULT_FIGURE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    combined = load_data(args.ordered, args.disordered, args.refs)
    combined.to_csv(args.long_output, index=False, float_format="%.9f", quoting=csv.QUOTE_MINIMAL)
    summary_df, report = summarize(combined)
    summary_df.to_csv(args.summary_csv, index=False, float_format="%.9f")
    args.summary_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    create_figure(combined, report, args.figure_base)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
