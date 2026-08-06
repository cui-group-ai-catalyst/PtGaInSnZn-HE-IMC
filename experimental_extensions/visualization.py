"""Generate static validation figures from manifest-selected modules."""
from __future__ import annotations

import math
from pathlib import Path


def _plot_manifold(ax, result: dict, panel_label: str, module_id: str) -> None:
    import numpy as np

    baseline = result["endmember_only_nonendmember_metrics"]["RMSE"]
    loocv = result["nonendmember_LOOCV_metrics"]["RMSE"]
    labels = ["Endmember-only\nbaseline", "Pairwise CEF\nLOOCV"]
    values = [baseline, loocv]
    bars = ax.bar(
        np.arange(2), values, width=0.62, color=["#6d7780", "#00897b"],
        edgecolor="black", linewidth=0.8,
    )
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + max(values) * 0.03,
            f"{value:.3f}",
            ha="center",
            va="bottom",
        )
    group = result.get("group_holdout")
    if group:
        ax.axhspan(
            group["RMSE_min"], group["RMSE_max"], color="#f9a825", alpha=0.22,
            label=f"{group['element']}-count group holdout range",
        )
        ax.legend(loc="upper right", frameon=False, fontsize=8)
    ax.set_xticks(np.arange(2), labels)
    unit = result.get("energy_unit", "reported energy unit")
    ax.set_ylabel(f"RMSE ({unit})")
    ax.set_ylim(0, max(values) * 1.25)
    ax.set_title(
        f"{panel_label}  Fixed-manifold interpolation", loc="left", fontweight="bold"
    )
    ax.text(
        0.02, 0.97, f"Internal validation only | {module_id}",
        transform=ax.transAxes, va="top", ha="left", color="#455a64", fontsize=8,
    )


def _plot_backend_comparison(
    ax, result: dict, panel_label: str, module_id: str
) -> None:
    import numpy as np

    backends = result["backends"]
    if len(backends) < 2:
        raise ValueError("An energy-backend comparison requires at least two backends")
    rows = result["ranking_rows"]
    reference = backends[0]
    reference_id = reference["id"]
    x = np.array([row[f"{reference_id}_rank"] for row in rows], dtype=float)
    all_subset = result["pairwise_results"][0]["subset_id"]
    rho_by_pair = {
        row["pair_id"]: row["spearman_rho"]
        for row in result["pairwise_results"]
        if row["subset_id"] == all_subset
    }
    colors = ["#c62828", "#1565c0", "#6a1b9a", "#00897b", "#ef6c00"]
    markers = ["o", "s", "D", "^", "v"]
    maximum_rank = float(max(x))
    for index, backend in enumerate(backends[1:]):
        backend_id = backend["id"]
        y = np.array([row[f"{backend_id}_rank"] for row in rows], dtype=float)
        maximum_rank = max(maximum_rank, float(max(y)))
        pair_id = f"{reference_id}__{backend_id}"
        rho = rho_by_pair.get(pair_id)
        label = backend.get("label", backend_id)
        if rho is not None:
            label += f" (rho={rho:.3f})"
        ax.scatter(
            x, y, s=36, c=colors[index % len(colors)],
            marker=markers[index % len(markers)], edgecolor="black",
            linewidth=0.45, label=label,
        )
    ax.plot(
        [1, maximum_rank], [1, maximum_rank], color="#777777",
        linestyle="--", linewidth=0.9,
    )
    ax.set_xlim(maximum_rank + 0.7, 0.3)
    ax.set_ylim(maximum_rank + 0.7, 0.3)
    ax.set_xlabel(f"{reference.get('label', reference_id)} rank")
    ax.set_ylabel("Comparison-backend rank")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.set_title(
        f"{panel_label}  Reference-set rank consistency",
        loc="left", fontweight="bold",
    )
    ax.text(
        0.02, 0.97, f"Reference-set comparison only | {module_id}",
        transform=ax.transAxes, va="top", ha="left", color="#455a64", fontsize=8,
    )


def write_validation_figure(results: dict, output_stem: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modules = list(results["modules"].items())
    if not modules:
        raise ValueError("No enabled modules were available for visualization")
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )
    ncols = min(2, len(modules))
    nrows = math.ceil(len(modules) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5.6 * ncols, 4.5 * nrows), squeeze=False
    )
    flat_axes = list(axes.flat)
    for index, ((module_id, module), ax) in enumerate(zip(modules, flat_axes)):
        panel_label = chr(ord("a") + index)
        if module["kind"] == "manifold_regression":
            _plot_manifold(ax, module["result"], panel_label, module_id)
        elif module["kind"] == "energy_backend_comparison":
            _plot_backend_comparison(ax, module["result"], panel_label, module_id)
        else:
            raise ValueError(f"Unsupported visualization module: {module['kind']}")
    for ax in flat_axes[len(modules):]:
        ax.set_visible(False)
    scope = str(results.get("scientific_scope", ""))
    if "synthetic" in scope.lower():
        fig.suptitle(
            "SYNTHETIC INTERFACE TEST - NO MATERIAL CLAIM",
            fontsize=10, color="#a61b1b", fontweight="bold",
        )
        fig.tight_layout(w_pad=2.0, rect=(0, 0, 1, 0.96))
    else:
        fig.tight_layout(w_pad=2.0)
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
