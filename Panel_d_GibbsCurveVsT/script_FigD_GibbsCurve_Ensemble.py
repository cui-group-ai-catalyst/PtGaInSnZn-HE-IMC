"""Build Panel-d Gibbs curves from 30 ordered and 30 disordered UMA structures.

The script reads the measured fixed-cell UMA single-point ensembles used by
the extended Panel c. Historical Panel-d canonical files are not overwritten.
Configurational-entropy curves are model bounds, not direct finite-temperature
free-energy calculations.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parent
PANEL_C = RELEASE_ROOT / "Panel_c_OrderedVsDisordered"
DEFAULT_ORDERED = PANEL_C / "data_FigC_OrderedEnsemble_Raw_UMA_Energies_regen.csv"
DEFAULT_DISORDERED = PANEL_C / "data_FigC_Raw_UMA_Energies.csv"
DEFAULT_REFS = RELEASE_ROOT / "shared" / "UMA_Element_Reference_Energies.csv"
DEFAULT_CURVES = HERE / "data_FigD_GibbsCurve_Ensemble_regen.csv"
DEFAULT_KEY_POINTS = HERE / "data_FigD_KeyPoints_Ensemble_regen.csv"
DEFAULT_REPORT = HERE / "validation_FigD_GibbsCurve_Ensemble_regen.json"
DEFAULT_FIGURE = HERE / "preview_FigD_GibbsCurve_Ensemble_regen"

COMPOSITION = {"Pt": 24, "Ga": 2, "In": 2, "Sn": 2, "Zn": 2}
EV_TO_KJMOL = 96.485
R_KJ_MOL_K = 8.31446261815324e-3
N_ATOMS = 32
B_SUBLATTICE_FRACTION = 8 / 32
KEY_TEMPERATURES = (500, 1000, 1500)


def element_reference_energy(path: Path) -> float:
    refs = pd.read_csv(path).set_index("Element")["UMA_E_eV_atom"].to_dict()
    missing = set(COMPOSITION) - set(refs)
    if missing:
        raise ValueError(f"Missing elemental reference energies: {sorted(missing)}")
    return sum(COMPOSITION[element] * refs[element] for element in COMPOSITION) / N_ATOMS


def load_ensembles(
    ordered_path: Path,
    disordered_path: Path,
    refs_path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    ordered = pd.read_csv(ordered_path)
    disordered_raw = pd.read_csv(disordered_path)
    disordered = disordered_raw.loc[
        disordered_raw["Type"] == "Disordered_Random"
    ].copy()
    if len(ordered) != 30 or len(disordered) != 30:
        raise ValueError(
            f"Expected 30 ordered and 30 disordered structures, got {len(ordered)} and {len(disordered)}"
        )
    if ordered["Canonical_Class"].nunique() != 30:
        raise ValueError("Ordered structures do not represent 30 symmetry-distinct classes")
    if ordered["Assignment_SHA256"].nunique() != 30:
        raise ValueError("Ordered structures contain duplicate occupancy fingerprints")
    if disordered["Seed"].astype(int).tolist() != list(range(100, 130)):
        raise ValueError("Disordered seeds must be exactly 100-129")

    reference = element_reference_energy(refs_path)
    ordered_hf = (ordered["Energy_eV_atom"].to_numpy(float) - reference) * EV_TO_KJMOL
    disordered_hf = (
        disordered["Energy_eV_atom"].to_numpy(float) - reference
    ) * EV_TO_KJMOL
    return ordered_hf, disordered_hf


def sample_summary(values: np.ndarray) -> dict[str, float | int]:
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


def entropy_disordered_full_lattice() -> float:
    fractions = np.array([24 / 32, 2 / 32, 2 / 32, 2 / 32, 2 / 32])
    return float(-R_KJ_MOL_K * np.sum(fractions * np.log(fractions)))


def entropy_ordered_b_sublattice_upper_bound() -> float:
    return float(B_SUBLATTICE_FRACTION * R_KJ_MOL_K * math.log(4.0))


def build_curves(
    ordered: np.ndarray,
    disordered: np.ndarray,
    temperatures: np.ndarray,
) -> tuple[pd.DataFrame, dict]:
    ordered_stats = sample_summary(ordered)
    disordered_stats = sample_summary(disordered)
    s_disordered = entropy_disordered_full_lattice()
    s_ordered_upper = entropy_ordered_b_sublattice_upper_bound()

    g_disordered = disordered_stats["mean"] - temperatures * s_disordered
    g_ordered_bmix = ordered_stats["mean"] - temperatures * s_ordered_upper
    g_ordered_frozen = np.full_like(temperatures, ordered_stats["mean"], dtype=float)
    gap_bmix = g_disordered - g_ordered_bmix
    gap_frozen = g_disordered - g_ordered_frozen

    curves = pd.DataFrame(
        {
            "T_K": temperatures,
            "G_disordered_mean_kJ_mol_atom": g_disordered,
            "G_disordered_minus1SD": g_disordered - disordered_stats["sd"],
            "G_disordered_plus1SD": g_disordered + disordered_stats["sd"],
            "G_ordered_Bmix_mean_kJ_mol_atom": g_ordered_bmix,
            "G_ordered_Bmix_minus1SD": g_ordered_bmix - ordered_stats["sd"],
            "G_ordered_Bmix_plus1SD": g_ordered_bmix + ordered_stats["sd"],
            "G_ordered_frozen_mean_kJ_mol_atom": g_ordered_frozen,
            "Gap_disordered_minus_ordered_Bmix": gap_bmix,
            "Gap_disordered_minus_ordered_frozen": gap_frozen,
        }
    )

    gap_0k = disordered_stats["mean"] - ordered_stats["mean"]
    gap_sem = math.sqrt(
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

    report = {
        "status": "passed",
        "provenance_mode": "measured_computational_plus_ideal_entropy_model",
        "energy_protocol": "fixed-cell fixed-coordinate UMA-s-1p1 single-point",
        "ordered": ordered_stats,
        "disordered": disordered_stats,
        "entropy_disordered_kJ_mol_atom_K": s_disordered,
        "entropy_ordered_Bmix_upper_bound_kJ_mol_atom_K": s_ordered_upper,
        "gap_0K_kJ_mol_atom": gap_0k,
        "gap_0K_welch_ci95_kJ_mol_atom": [
            gap_0k - t_critical * gap_sem,
            gap_0k + t_critical * gap_sem,
        ],
        "crossover_Bmix_K": gap_0k / (s_disordered - s_ordered_upper),
        "crossover_frozen_K": gap_0k / s_disordered,
        "interpretation_boundary": (
            "Temperature curves add ideal configurational-entropy bounds to 0 K UMA ensemble means; "
            "they are not phonon-inclusive or directly sampled finite-temperature free energies."
        ),
    }
    return curves, report


def create_figure(
    curves: pd.DataFrame,
    ordered: np.ndarray,
    disordered: np.ndarray,
    report: dict,
    output_base: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    ordered_color = "#0072B2"
    disordered_color = "#D55E00"
    fig, ax = plt.subplots(figsize=(5.4, 3.75), constrained_layout=True)
    t = curves["T_K"].to_numpy()

    ax.axvspan(500, 1500, color="#6A8F3A", alpha=0.08, zorder=0)
    ax.fill_between(
        t,
        curves["G_disordered_minus1SD"],
        curves["G_disordered_plus1SD"],
        color=disordered_color,
        alpha=0.13,
        linewidth=0,
        label=r"Disordered ensemble ($\pm 1$ SD)",
    )
    ax.fill_between(
        t,
        curves["G_ordered_Bmix_minus1SD"],
        curves["G_ordered_Bmix_plus1SD"],
        color=ordered_color,
        alpha=0.16,
        linewidth=0,
        label=r"Ordered ensemble ($\pm 1$ SD)",
    )
    ax.plot(
        t,
        curves["G_disordered_mean_kJ_mol_atom"],
        color=disordered_color,
        linewidth=1.8,
        label="Disordered: full-lattice ideal mixing",
    )
    ax.plot(
        t,
        curves["G_ordered_Bmix_mean_kJ_mol_atom"],
        color=ordered_color,
        linewidth=1.8,
        label="Ordered: B-sublattice mixing bound",
    )
    ax.plot(
        t,
        curves["G_ordered_frozen_mean_kJ_mol_atom"],
        color=ordered_color,
        linewidth=1.0,
        linestyle="--",
        label="Ordered: frozen-occupancy bound",
    )

    for temperature in KEY_TEMPERATURES:
        row = curves.loc[curves["T_K"] == temperature].iloc[0]
        gap = float(row["Gap_disordered_minus_ordered_Bmix"])
        g_dis = float(row["G_disordered_mean_kJ_mol_atom"])
        g_ord = float(row["G_ordered_Bmix_mean_kJ_mol_atom"])
        ax.plot([temperature, temperature], [g_ord, g_dis], color="#4A4A4A", lw=0.65)
        ax.text(
            temperature + 45,
            0.5 * (g_ord + g_dis),
            f"{gap:.1f}",
            va="center",
            ha="left",
            fontsize=7.2,
            color="#303030",
        )

    ax.text(
        1000,
        -11.2,
        "synthesis window",
        ha="center",
        va="top",
        fontsize=7.2,
        color="#496525",
    )
    ax.set_xlim(0, 3600)
    ax.set_ylim(-43.5, -8.0)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"Modeled $G_f$ (kJ mol$^{-1}$ atom$^{-1}$)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="in", length=3)
    ax.legend(loc="lower left", frameon=False, ncol=1)

    inset = ax.inset_axes([0.68, 0.54, 0.28, 0.36])
    rng = np.random.default_rng(20260727)
    for position, values, color in (
        (0, ordered, ordered_color),
        (1, disordered, disordered_color),
    ):
        jitter = rng.uniform(-0.13, 0.13, len(values))
        inset.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=10,
            color=color,
            edgecolor="white",
            linewidth=0.25,
            alpha=0.85,
        )
        inset.plot(
            [position - 0.18, position + 0.18],
            [np.mean(values), np.mean(values)],
            color="#202020",
            lw=1.0,
        )
    inset.set_xticks([0, 1], ["Ord.", "Dis."])
    inset.set_ylabel(r"$Delta H_f$ at 0 K", fontsize=6.8)
    inset.tick_params(labelsize=6.4, direction="in", length=2)
    inset.spines["top"].set_visible(False)
    inset.spines["right"].set_visible(False)
    inset.text(
        0.5,
        1.04,
        "n = 30 per group",
        transform=inset.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.5,
    )

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
    parser.add_argument("--curves-output", type=Path, default=DEFAULT_CURVES)
    parser.add_argument("--key-points-output", type=Path, default=DEFAULT_KEY_POINTS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--figure-base", type=Path, default=DEFAULT_FIGURE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ordered, disordered = load_ensembles(args.ordered, args.disordered, args.refs)
    temperatures = np.arange(0.0, 3600.0 + 10.0, 10.0)
    curves, report = build_curves(ordered, disordered, temperatures)
    key_points = curves.loc[curves["T_K"].isin(KEY_TEMPERATURES)].copy()
    args.curves_output.parent.mkdir(parents=True, exist_ok=True)
    curves.to_csv(args.curves_output, index=False)
    key_points.to_csv(args.key_points_output, index=False)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    create_figure(curves, ordered, disordered, report, args.figure_base)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
