"""
20260418_FigG_OriginReady_ElementRef.py
========================================
Panel g 修订版（2026-04-18）——切到 element-reference 轴 + Origin 宽表。

修订动因
--------
- 旧 Panel g 的 Y 轴基于 Calibrated_Display = (E + 4.356) * 96.485，
  其中 4.356 eV/atom 的基线没有明确物理出处。
- 切换为 element-referenced formation enthalpy：
      ΔH_f = E_alloy - Σ x_i · E_i^elem
  element reference 来自 20260417_UMA_Element_Reference_Energies.csv
- 切换后全局最低点从 Ga5Sn3 变为 Ga8 Pt3Ga-baseline (-43.38)，
  但 Ga-rich 端依然是整条低焓脊；结论仍为"Ga-rich 窗口有利"。

数据源
------
- data_FigG_165_ElementReferenced_Hf.csv  (165 prototypes, both columns present)

输出
----
03_results/20260418_FigG_OriginReady_Scatter_byCategory.csv
    Origin 散点宽表：每个 Category 占一对 X/Y 列
    （保留和旧 20260415_FigB_Origin_Compact 类似的多分组布局）
03_results/20260418_FigG_OriginReady_Stats.csv
    按 Ga_pct 的统计：Mean/Std/Min/Max/N，用于 Origin 误差棒序列
03_results/20260418_FigG_OriginReady_RefPoints.csv
    关键标注点（Ga8 baseline / Ga6InSn model-selected / Ga5In2Sn experimental-rep）
03_results/20260418_FigG_OriginReady_ZnSweep.csv
    Zn 副轴子图用：Ga_{8-k}Zn_k 一条，沿 Zn% 单调恶化
03_results/20260418_FigG_OriginReady_LowBand.csv
    Ga-rich 低焓带底色区的边界（[-43, -35] kJ/mol，Ga% 62.5-100%）
03_results/20260418_FigG_Preview_Plot.png
    快速预览图（供核对，最终请在 Origin 重绘）
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SRC = SCRIPT_DIR / "data_FigG_165_ElementReferenced_Hf.csv"

# 分类顺序：与旧图保持一致
CATEGORY_ORDER = [
    "Pure Pt3Ga",
    "Pure Pt3Sn",
    "Pure Pt3In",
    "Pure Pt3Zn",
    "Ga-Sn binary",
    "Ga-In binary",
    "In-Sn binary",
    "Ga-In-Sn ternary",
    "Zn-binary/ternary",
    "Quaternary",
]

# Origin 友好的列名前缀（ASCII only）
CATEGORY_SLUG = {
    "Pure Pt3Ga":          "PurePt3Ga",
    "Pure Pt3Sn":          "PurePt3Sn",
    "Pure Pt3In":          "PurePt3In",
    "Pure Pt3Zn":          "PurePt3Zn",
    "Ga-Sn binary":        "GaSn",
    "Ga-In binary":        "GaIn",
    "In-Sn binary":        "InSn",
    "Ga-In-Sn ternary":    "GaInSn",
    "Zn-binary/ternary":   "ZnBinTer",
    "Quaternary":          "Quaternary",
}

# 源 CSV 里旧标签 -> 新标签（Zn 亚组定义不变，只改名）
CATEGORY_RENAME = {
    "Contains Zn": "Zn-binary/ternary",
}


def build_scatter_wide(df: pd.DataFrame) -> pd.DataFrame:
    """
    把 165 组按 Category 拆成成对 X/Y 列，NaN 填充对齐。
    X 列已经叠加了 deterministic jitter（±2.3 Ga%），
    Origin 直接用 XY 散点就能看到横向散布，不需要再额外抖动。
    同时保留 X_raw（真实 Ga%）和 Composition，方便核对。
    """
    rng = np.random.default_rng(20260418)  # 固定种子，可复现
    cols = {}
    max_len = 0
    for cat in CATEGORY_ORDER:
        sub = df[df["Category"] == cat].reset_index(drop=True)
        slug = CATEGORY_SLUG[cat]
        ga_raw = sub["Ga_pct"].to_numpy(dtype=float)
        jitter = rng.uniform(-2.3, 2.3, size=len(sub))
        cols[f"{slug}_Ga_pct_jitter"] = ga_raw + jitter
        cols[f"{slug}_Hf_kJ_mol"]     = sub["ElementRef_Hf_kJ_mol"].to_numpy(dtype=float)
        cols[f"{slug}_Ga_pct_raw"]    = ga_raw
        cols[f"{slug}_Composition"]   = sub["Composition"].to_numpy(dtype=object)
        max_len = max(max_len, len(sub))

    padded = {}
    for k, v in cols.items():
        pad = max_len - len(v)
        if pad > 0:
            if v.dtype == object:
                v = np.concatenate([v, np.full(pad, "", dtype=object)])
            else:
                v = np.concatenate([v, np.full(pad, np.nan)])
        padded[k] = v

    return pd.DataFrame(padded)


def build_stats(df: pd.DataFrame) -> pd.DataFrame:
    """按 Ga_pct 聚合 ElementRef_Hf 的 mean/std/min/max/N。"""
    stats = (
        df.groupby("Ga_pct")["ElementRef_Hf_kJ_mol"]
        .agg(["mean", "std", "count", "min", "max"])
        .reset_index()
        .rename(columns={
            "mean":  "Mean_Hf_kJ_mol",
            "std":   "Std_Hf_kJ_mol",
            "count": "N_configs",
            "min":   "Min_Hf_kJ_mol",
            "max":   "Max_Hf_kJ_mol",
        })
    )
    stats["Range_kJ_mol"] = stats["Max_Hf_kJ_mol"] - stats["Min_Hf_kJ_mol"]
    return stats


def build_refpoints(df: pd.DataFrame) -> pd.DataFrame:
    """
    挑 4 个关键点做图上标注：
      1) Ga8 Pt3Ga-baseline              --- 新的 element-ref 全局最低
      2) Ga5Sn3 GLOBAL-MIN (历史)        --- 老 calibrated 轴的全局最低，对比用
      3) Ga6InSn MODEL-SELECTED          --- 原先模型选点
      4) Ga5In2Sn experimental-rep       --- 32 原子最接近实验 Ga5.42In1.70Sn0.78 的代表点
    """
    picks = []
    # (Ga_count, In_count, Sn_count, Zn_count, label)
    targets = [
        ((8, 0, 0, 0), "Ga8 Pt3Ga-baseline (ElementRef global min)"),
        ((5, 0, 3, 0), "Ga5Sn3 (historical calibrated global min)"),
        ((6, 1, 1, 0), "Ga6InSn (model-selected)"),
        ((5, 2, 1, 0), "Ga5In2Sn (experimental-representative)"),
    ]
    for (ng, ni, ns, nz), label in targets:
        row = df[
            (df["Ga_count"] == ng)
            & (df["In_count"] == ni)
            & (df["Sn_count"] == ns)
            & (df["Zn_count"] == nz)
        ]
        if len(row) != 1:
            raise ValueError(
                f"Reference point lookup failed for ({ng},{ni},{ns},{nz}); "
                f"got {len(row)} rows"
            )
        r = row.iloc[0]
        picks.append({
            "Label":          label,
            "Composition":    r["Composition"],
            "Ga_pct":         float(r["Ga_pct"]),
            "Hf_kJ_mol":      float(r["ElementRef_Hf_kJ_mol"]),
            "Rank_overall":   int(r["Rank"]),
            "Ga_count":       int(r["Ga_count"]),
            "In_count":       int(r["In_count"]),
            "Sn_count":       int(r["Sn_count"]),
            "Zn_count":       int(r["Zn_count"]),
        })
    return pd.DataFrame(picks)


def build_zn_sweep(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ga_{8-k}Zn_k 沿 Zn 含量的 ΔH_f 序列（9 点，含 Ga8 和 Zn8 两个端点）。
    供叙事"少量 Zn 代价小、多量 Zn 代价大"的插图使用。
    """
    rows = []
    for k in range(9):
        sub = df[
            (df["Ga_count"] == 8 - k)
            & (df["In_count"] == 0)
            & (df["Sn_count"] == 0)
            & (df["Zn_count"] == k)
        ]
        if len(sub) != 1:
            raise ValueError(
                f"Zn-sweep lookup failed at Ga_{8-k}Zn_{k}; got {len(sub)} rows"
            )
        r = sub.iloc[0]
        rows.append({
            "Zn_count":         int(k),
            "Zn_pct":           float(r["Zn_pct"]),
            "Ga_pct":           float(r["Ga_pct"]),
            "Composition":      r["Composition"],
            "Hf_kJ_mol":        float(r["ElementRef_Hf_kJ_mol"]),
        })
    dfr = pd.DataFrame(rows)
    dfr["Penalty_vs_Ga8_kJ_mol"] = dfr["Hf_kJ_mol"] - dfr.iloc[0]["Hf_kJ_mol"]
    return dfr


def build_lowband() -> pd.DataFrame:
    """Ga-rich 低焓带底色，矩形 4 顶点。"""
    # [Ga% low, Ga% high, Hf low (更负), Hf high (更浅)]
    return pd.DataFrame([
        {"X_Ga_pct": 62.5,  "Y_Hf_kJ_mol": -43.5, "Vertex": "bottom-left"},
        {"X_Ga_pct": 100.0, "Y_Hf_kJ_mol": -43.5, "Vertex": "bottom-right"},
        {"X_Ga_pct": 100.0, "Y_Hf_kJ_mol": -34.5, "Vertex": "top-right"},
        {"X_Ga_pct": 62.5,  "Y_Hf_kJ_mol": -34.5, "Vertex": "top-left"},
    ])


def plot_preview(df: pd.DataFrame, stats: pd.DataFrame, refs: pd.DataFrame,
                 out_path: Path) -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 9,
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })

    cat_styles = {
        "Pure Pt3Ga":        ("#111111", "*", 85),
        "Pure Pt3Sn":        ("#1565C0", "*", 85),
        "Pure Pt3In":        ("#2E7D32", "*", 85),
        "Pure Pt3Zn":        ("#7B1FA2", "*", 85),
        "Ga-Sn binary":      ("#1E88E5", "o", 28),
        "Ga-In binary":      ("#43A047", "s", 28),
        "In-Sn binary":      ("#66BB6A", "d", 28),
        "Ga-In-Sn ternary":  ("#FB8C00", "^", 30),
        "Zn-binary/ternary": ("#AB47BC", "v", 26),
        "Quaternary":        ("#78909C", "D", 26),
    }

    fig, ax = plt.subplots(figsize=(11, 6.8))

    # 低焓带底色
    ax.axhspan(-43.5, -34.5, xmin=0.625, xmax=1.0, alpha=0.10,
               color="#2E7D32", zorder=0)
    # 注：axhspan 的 xmin/xmax 是 axes 比例，不是数据坐标 —— 后面会 set_xlim

    # 散点，按 Ga% 水平抖动（与 CSV 里的 jitter 种子一致，保证一一对应）
    rng = np.random.default_rng(20260418)
    for cat in CATEGORY_ORDER:
        sub = df[df["Category"] == cat].reset_index(drop=True)
        if len(sub) == 0:
            continue
        c, m, s = cat_styles.get(cat, ("#999999", "o", 22))
        jitter = rng.uniform(-2.3, 2.3, size=len(sub))
        ax.scatter(
            sub["Ga_pct"].to_numpy() + jitter,
            sub["ElementRef_Hf_kJ_mol"].to_numpy(),
            c=c, marker=m, s=s, alpha=0.78, edgecolors="none", zorder=3,
            label=cat,
        )

    # 误差棒（Mean ± Std per Ga%）
    for _, row in stats.iterrows():
        if row["N_configs"] > 1 and not np.isnan(row["Std_Hf_kJ_mol"]):
            ax.errorbar(
                row["Ga_pct"], row["Mean_Hf_kJ_mol"],
                yerr=row["Std_Hf_kJ_mol"],
                fmt="D", color="#C62828", markersize=6, capsize=4,
                capthick=1.2, elinewidth=1.2,
                markerfacecolor="white", markeredgecolor="#C62828",
                markeredgewidth=1.3, zorder=10,
            )
        else:
            ax.plot(row["Ga_pct"], row["Mean_Hf_kJ_mol"], "D",
                    color="#C62828", markersize=6,
                    markerfacecolor="white", markeredgecolor="#C62828",
                    markeredgewidth=1.3, zorder=10)

    # 参考点（SI 用，不画在主预览图上 —— 最终图只保留 ensemble trend + 窗口）
    # 如需画出，参见 RefPoints.csv
    del refs  # 主图不再使用，仅 CSV 留档

    # 实验 Ga% 垂直线
    ax.axvline(68.5, color="#D32F2F", linestyle="--", linewidth=1.3, alpha=0.8)
    ax.text(68.5, ax.get_ylim()[0] if False else -45.5,
            "Expt. 68.5% Ga", color="#D32F2F", fontsize=8,
            fontweight="bold", rotation=0, ha="center", va="top")

    ax.set_xlim(-5, 105)
    ax.set_ylim(-46, -18)
    # Origin 版本 X 轴用 jitter 后的连续 Ga%，这里预览图也用连续 X 展示散布
    ax.set_xticks([0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100])
    ax.set_xticklabels(
        [f"{g:.1f}%" if g % 1 else f"{int(g)}%"
         for g in [0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100]],
        fontsize=8,
    )
    ax.set_xlabel("Ga fraction on the 8-site B-sublattice (%)",
                  fontsize=11, fontweight="bold")
    ax.set_ylabel(r"Formation enthalpy, $\Delta H_f$ (kJ mol$^{-1}$)",
                  fontsize=11, fontweight="bold")
    ax.set_title(
        "Panel g preview (2026-04-18) | Element-referenced formation enthalpy\n"
        "165 stoichiometric B-sublattice prototypes (Pt3X8), UMA-s-1p1 single-point",
        fontsize=10, fontweight="bold", pad=8,
    )
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    ax.legend(loc="upper right", fontsize=7.3, framealpha=0.94,
              edgecolor="#CCCCCC", ncol=1)

    # 低焓带注释
    ax.text(
        81.0, -40.5,
        "Ga-rich favourable window\n(62.5% - 100% Ga)",
        color="#2E7D32", fontsize=8.5, fontweight="bold",
        ha="center", va="center",
        bbox=dict(facecolor="white", alpha=0.88, edgecolor="#2E7D32", lw=0.8),
    )

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    if not SRC.exists():
        raise FileNotFoundError(f"Missing source CSV: {SRC}")

    df = pd.read_csv(SRC)
    if len(df) != 165:
        raise ValueError(f"Expected 165 rows, got {len(df)}")

    # 旧标签 -> 新标签 (Zn 亚组只改名，定义不变)
    df["Category"] = df["Category"].replace(CATEGORY_RENAME)

    # 1) Category-split scatter wide
    scatter = build_scatter_wide(df)
    p1 = RESULTS_DIR / "data_FigG_OriginReady_Scatter_byCategory_regen.csv"
    scatter.to_csv(p1, index=False, float_format="%.6f")
    print(f"[Panel g] scatter (by category)    -> {p1}")

    # 2) Stats per Ga_pct
    stats = build_stats(df)
    p2 = RESULTS_DIR / "data_FigG_OriginReady_Stats_regen.csv"
    stats.to_csv(p2, index=False, float_format="%.6f")
    print(f"[Panel g] stats per Ga_pct         -> {p2}")

    # 3) Reference points
    refs = build_refpoints(df)
    p3 = RESULTS_DIR / "data_FigG_OriginReady_RefPoints_regen.csv"
    refs.to_csv(p3, index=False, float_format="%.6f")
    print(f"[Panel g] ref points               -> {p3}")

    # 4) Zn sweep (Ga_{8-k}Zn_k)
    zn = build_zn_sweep(df)
    p4 = RESULTS_DIR / "data_FigG_OriginReady_ZnSweep_regen.csv"
    zn.to_csv(p4, index=False, float_format="%.6f")
    print(f"[Panel g] Zn sweep along GaZn axis -> {p4}")

    # 5) Low band rectangle vertices
    band = build_lowband()
    p5 = RESULTS_DIR / "data_FigG_OriginReady_LowBand_regen.csv"
    band.to_csv(p5, index=False, float_format="%.2f")
    print(f"[Panel g] low-band rectangle       -> {p5}")

    # 6) Preview plot
    p6 = RESULTS_DIR / "preview_FigG_GaSweep_regen.png"
    plot_preview(df, stats, refs, p6)
    print(f"[Panel g] preview plot             -> {p6}")

    # Summary print
    print("\n-- Panel g key numbers (element-reference axis) --")
    print(f"  Global min (overall)        : Rank 1 = "
          f"{df.iloc[0]['Composition']} at "
          f"{df.iloc[0]['ElementRef_Hf_kJ_mol']:.3f} kJ/mol")
    for _, r in refs.iterrows():
        print(f"  {r['Label']:55s} : "
              f"{r['Composition']} @ Ga%={r['Ga_pct']:.1f}, "
              f"Hf={r['Hf_kJ_mol']:.3f} kJ/mol (rank {int(r['Rank_overall'])})")
    print("\n-- Zn sweep Ga_{8-k}Zn_k (penalty vs Ga8) --")
    print(zn[["Composition", "Zn_pct", "Hf_kJ_mol",
              "Penalty_vs_Ga8_kJ_mol"]].to_string(index=False))


if __name__ == "__main__":
    main()
