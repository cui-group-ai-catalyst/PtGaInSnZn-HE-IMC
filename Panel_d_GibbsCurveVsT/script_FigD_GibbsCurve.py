"""
SI — HEI vs HEA Gibbs-energy temperature curves
================================================

把 Panel c 的 0 K ordering gap (16.04 kJ/mol) 沿温度轴推广为连续的
ΔG(T) 曲线，直接展示在 manuscript 合成温度区间 (500-1500 K) 内，
HEI (L1₂ 有序 Pt24Ga2In2Sn2Zn2) 的 Gibbs 自由能**始终低于** HEA
(相同组成 disordered 固溶体) → 支撑 "HEI 更倾向被合成" 的论点。

Inputs
------
Panel c 已有的 UMA 0 K element-referenced 形成焓 (32-原子 supercell)：
  - Ordered   L1₂            Pt24Ga2In2Sn2Zn2  :  ΔH_f = -30.008 kJ/mol/atom
  - Disordered random (N=30) Pt24Ga2In2Sn2Zn2  :  ΔH_f = -13.964 ± 2.262 kJ/mol/atom
Read from `Panel_c_OrderedVsDisordered/data_FigC_Long.csv` (no UMA re-run).

Thermodynamic model
-------------------
对相同组成 x = (x_Pt=0.75, x_Ga=x_In=x_Sn=x_Zn=0.0625)，
以纯元素为参考态：

  ΔG^HEA(T) = <ΔH^HEA>  −  T · S_config^HEA
  ΔG^HEI(T) =  ΔH^HEI   −  T · S_config^HEI

其中：
  - S_config^HEA = -R Σ_i x_i ln x_i        (全 32 位点完全无序)
  - S_config^HEI (two bounds):
       * 保守上界 (sublattice model): B 子晶格 8/32 位点上 4 元素等分的残余
         熵 → (8/32)·(-R Σ_j (1/4) ln(1/4)) = 0.25 R ln 4
       * 严格下界: HEI 取 Panel c 的具体有序占位 (single frozen config)
         → S_config^HEI = 0

Outputs
-------
  data_SI_HEIvsHEA_GibbsCurve.csv    — T × {ΔG_HEA, ±1σ band, ΔG_HEI(upper), ΔG_HEI(lower)}
  data_SI_HEIvsHEA_KeyPoints.csv     — T = 500 / 1000 / 1500 K 的 ΔG 差
  preview_SI_HEIvsHEA_GibbsCurve.png — T vs ΔG 曲线 + 交叉温度标注

Env: py312 (numpy, pandas, matplotlib)
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Windows: force UTF-8 stdout so subscripts / Δ / ² 等字符不触发 GBK 编码错误
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------------------------------------------------------------- paths
HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parent                    # code_release_v2/
PANEL_C_LONG = (
    RELEASE_ROOT
    / "Panel_c_OrderedVsDisordered"
    / "data_FigC_Long.csv"
)

# ---------------------------------------------------------------- composition
# 32-atom supercell: Pt24 Ga2 In2 Sn2 Zn2
X_PT = 24 / 32
X_X = 2 / 32           # each of Ga / In / Sn / Zn
COMPOSITION_LABEL = r"Pt$_{24}$(Ga,In,Sn,Zn)$_{2,2,2,2}$"

R_GAS = 8.314_462_618          # J / (mol·K)
R_kJ = R_GAS / 1000.0          # kJ / (mol·K)

# B-sublattice sites per atom (8 B-sites out of 32 total)
B_SUBLATTICE_FRAC = 8 / 32


def s_config_hea() -> float:
    """Full random solid solution: S = -R Σ x_i ln x_i  [kJ/(mol·K)]."""
    xs = np.array([X_PT] + [X_X] * 4)
    return -R_kJ * np.sum(xs * np.log(xs))


def s_config_hei_sublattice() -> float:
    """
    Conservative upper bound for HEI configurational entropy:
    B-sublattice carries residual 4-way mixing (Ga,In,Sn,Zn equimolar).

        S = (n_B / n_total) · (-R Σ_j y_j ln y_j)
          = (8/32) · R ln 4
    """
    return B_SUBLATTICE_FRAC * R_kJ * np.log(4.0)


def s_config_hei_zero() -> float:
    """HEI taken as single frozen occupancy (Panel c's representative config)."""
    return 0.0


# ---------------------------------------------------------------- load Panel c
def load_panel_c_enthalpies() -> tuple[float, float, float, int, np.ndarray]:
    df = pd.read_csv(PANEL_C_LONG)
    ord_mask = df["Type"] == "Ordered_L12_Equimolar"
    dis_mask = df["Type"] == "Disordered_Random_Equimolar"
    dh_ord = float(df.loc[ord_mask, "ZeroK_ElementRef_Hf_kJ_mol"].iloc[0])
    dh_dis_vals = df.loc[dis_mask, "ZeroK_ElementRef_Hf_kJ_mol"].to_numpy(float)
    return (
        dh_ord,
        float(dh_dis_vals.mean()),
        float(dh_dis_vals.std(ddof=0)),
        int(len(dh_dis_vals)),
        dh_dis_vals,
    )


# ---------------------------------------------------------------- curves
def build_curves(
    dh_hei: float,
    dh_hea_mean: float,
    dh_hea_std: float,
    t_grid: np.ndarray,
) -> pd.DataFrame:
    s_hea = s_config_hea()
    s_hei_up = s_config_hei_sublattice()      # lower bound on ΔG_HEI (more negative)
    s_hei_lo = s_config_hei_zero()            # upper bound on ΔG_HEI (less negative)

    dg_hea_mean = dh_hea_mean - t_grid * s_hea
    dg_hea_plus = (dh_hea_mean + dh_hea_std) - t_grid * s_hea
    dg_hea_minus = (dh_hea_mean - dh_hea_std) - t_grid * s_hea
    dg_hei_sub = dh_hei - t_grid * s_hei_up
    dg_hei_frozen = dh_hei - t_grid * s_hei_lo

    return pd.DataFrame({
        "T_K": t_grid,
        "dG_HEA_mean_kJmol_atom": dg_hea_mean,
        "dG_HEA_plus1sigma": dg_hea_plus,
        "dG_HEA_minus1sigma": dg_hea_minus,
        "dG_HEI_sublattice_upper_bound_entropy": dg_hei_sub,
        "dG_HEI_frozen_zero_entropy": dg_hei_frozen,
        "Gap_HEA_minus_HEI_sublattice": dg_hea_mean - dg_hei_sub,
        "Gap_HEA_minus_HEI_frozen": dg_hea_mean - dg_hei_frozen,
    })


def crossover_temperature(dh_gap: float, ds_gap: float) -> float:
    """T* where ΔG_HEA = ΔG_HEI, assuming linear curves."""
    if ds_gap <= 0:
        return np.inf
    return dh_gap / ds_gap


# ---------------------------------------------------------------- plot
def plot_curves(
    df: pd.DataFrame,
    t_star_sub: float,
    t_star_frozen: float,
    out_png: Path,
    hea_raw_configs: np.ndarray,   # 30 × ΔH_f at 0 K
    dh_hei_single: float,          # 1 ordered value at 0 K
) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.6))

    # --- synthesis-window shading (500-1500 K)
    ax.axvspan(500, 1500, color="#43A047", alpha=0.08, zorder=0)
    ax.text(
        1000, ax.get_ylim()[0] if False else -12.5,  # placeholder; adjusted after ylim
        "manuscript synthesis window\n500–1500 K",
        ha="center", va="top", fontsize=8.5, color="#1B5E20",
        style="italic",
    )

    # --- HEA band + mean
    ax.fill_between(
        df["T_K"],
        df["dG_HEA_minus1sigma"],
        df["dG_HEA_plus1sigma"],
        color="#90A4AE", alpha=0.28, label="HEA ±1σ (N=30 configs)",
    )
    ax.plot(
        df["T_K"], df["dG_HEA_mean_kJmol_atom"],
        color="#37474F", lw=2.2, label=r"HEA disordered $\langle\Delta G\rangle$",
    )

    # --- HEI two bounds
    ax.plot(
        df["T_K"], df["dG_HEI_sublattice_upper_bound_entropy"],
        color="#C62828", lw=2.2,
        label=r"HEI L1$_2$ (with B-sublattice $S_{\rm conf}$)",
    )
    ax.plot(
        df["T_K"], df["dG_HEI_frozen_zero_entropy"],
        color="#C62828", lw=1.3, ls="--",
        label=r"HEI L1$_2$ (frozen occupancy, $S=0$)",
    )

    # --- Panel c T=0 单点：HEI 只有 1 个结构（红菱形）；HEA 的 30 个
    #     disordered 散点已经由 mean 线 + ±1σ 灰带代表，不再重复画
    ax.scatter(
        [0], [dh_hei_single],
        s=110, marker="D", color="#C62828", edgecolors="white", linewidths=1.0,
        zorder=4, label=f"HEI ordered (Panel c, N=1; HEA N={len(hea_raw_configs)} → band)",
    )
    # Panel c 锚点 — 标 +16.04 @ 0 K 和 +10.43 @ 1200 K 的 gap
    for t_anchor, gap_anchor in [(0, 16.04), (1200, 10.43)]:
        g_hea = float(np.interp(t_anchor, df["T_K"], df["dG_HEA_mean_kJmol_atom"]))
        g_hei = float(np.interp(t_anchor, df["T_K"],
                                df["dG_HEI_sublattice_upper_bound_entropy"]))
        ax.annotate(
            f"Panel c anchor\nT={t_anchor} K,  Δ=+{gap_anchor:.2f}",
            xy=(t_anchor, 0.5 * (g_hea + g_hei)),
            xytext=(t_anchor + 200, 0.5 * (g_hea + g_hei) + 6),
            fontsize=8.5, color="#000000",
            arrowprops=dict(arrowstyle="->", color="#000", lw=0.8),
        )

    # --- key T points (500 / 1000 / 1500 K) with gap annotations
    key_colors = {500: "#1565C0", 1000: "#6A1B9A", 1500: "#E65100"}
    key_rows = df[df["T_K"].isin([500, 1000, 1500])]
    for _, row in key_rows.iterrows():
        t = float(row["T_K"])
        g_hea = float(row["dG_HEA_mean_kJmol_atom"])
        g_hei = float(row["dG_HEI_sublattice_upper_bound_entropy"])
        gap = float(row["Gap_HEA_minus_HEI_sublattice"])
        c = key_colors[int(t)]
        # vertical connector
        ax.plot([t, t], [g_hea, g_hei], color=c, lw=1.4, alpha=0.9, zorder=3)
        # endpoint dots
        ax.scatter([t, t], [g_hea, g_hei], s=40, color=c, zorder=4,
                   edgecolors="white", linewidths=0.8)
        # gap label at midpoint
        y_mid = 0.5 * (g_hea + g_hei)
        ax.annotate(
            f"T={int(t)} K\nΔ = +{gap:.1f}",
            xy=(t, y_mid), xytext=(t + 60, y_mid),
            color=c, fontsize=9, fontweight="bold",
            va="center", ha="left",
        )

    # --- crossover markers
    if np.isfinite(t_star_frozen):
        ax.axvline(t_star_frozen, color="#1B5E20", lw=1.0, ls=":",
                   label=f"T* (frozen) ≈ {t_star_frozen:.0f} K")
    if np.isfinite(t_star_sub):
        ax.axvline(t_star_sub, color="#2E7D32", lw=1.0, ls="-.",
                   label=f"T* (sublattice) ≈ {t_star_sub:.0f} K")

    # --- axes
    ax.set_xlim(-80, 3600)
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$\Delta G_f$ (kJ mol$^{-1}$ atom$^{-1}$)")
    ax.set_title(f"HEI vs HEA Gibbs-energy curves — {COMPOSITION_LABEL}")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.93)
    ax.grid(alpha=0.25)

    # re-position the green-window label now that ylim is known
    ymin, ymax = ax.get_ylim()
    ax.texts[0].set_position((1000, ymax - 0.5))

    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main() -> None:
    dh_hei, dh_hea_mean, dh_hea_std, n_dis, hea_raw = load_panel_c_enthalpies()
    print(f"Panel c inputs (kJ/mol/atom):")
    print(f"  HEI (ordered L1₂)        ΔH_f = {dh_hei:+.4f}")
    print(f"  HEA (disordered, N={n_dis}) ΔH_f = {dh_hea_mean:+.4f} ± {dh_hea_std:.4f}")
    print(f"  Gap (HEA - HEI)          ΔH_gap = {dh_hea_mean - dh_hei:+.4f}")

    s_hea = s_config_hea()
    s_hei_sub = s_config_hei_sublattice()
    s_hei_fro = s_config_hei_zero()
    ds_sub = s_hea - s_hei_sub
    ds_fro = s_hea - s_hei_fro
    print(f"\nConfigurational entropies (kJ/mol/atom/K):")
    print(f"  S_HEA              = {s_hea:.6f}")
    print(f"  S_HEI (sublattice) = {s_hei_sub:.6f}")
    print(f"  S_HEI (frozen)     = {s_hei_fro:.6f}")

    dh_gap = dh_hea_mean - dh_hei
    t_star_sub = crossover_temperature(dh_gap, ds_sub)
    t_star_fro = crossover_temperature(dh_gap, ds_fro)
    print(f"\nCrossover temperatures (where ΔG_HEA = ΔG_HEI):")
    print(f"  T* (sublattice entropy on HEI) ≈ {t_star_sub:.0f} K")
    print(f"  T* (frozen occupancy HEI)      ≈ {t_star_fro:.0f} K")

    t_grid = np.linspace(0, 3600, 361)
    df = build_curves(dh_hei, dh_hea_mean, dh_hea_std, t_grid)
    df.to_csv(HERE / "data_FigD_GibbsCurve_regen.csv", index=False)

    key = df[df["T_K"].isin([500, 1000, 1500])].reset_index(drop=True)
    key.to_csv(HERE / "data_FigD_KeyPoints_regen.csv", index=False)
    print(f"\nKey points at T = 500 / 1000 / 1500 K:")
    print(key[[
        "T_K",
        "dG_HEA_mean_kJmol_atom",
        "dG_HEI_sublattice_upper_bound_entropy",
        "Gap_HEA_minus_HEI_sublattice",
    ]].to_string(index=False))

    plot_curves(
        df, t_star_sub, t_star_fro,
        HERE / "preview_FigD_GibbsCurve_regen.png",
        hea_raw_configs=hea_raw,
        dh_hei_single=dh_hei,
    )
    print(f"\nWritten:")
    print(f"  {HERE / 'data_FigD_GibbsCurve_regen.csv'}")
    print(f"  {HERE / 'data_FigD_KeyPoints_regen.csv'}")
    print(f"  {HERE / 'preview_FigD_GibbsCurve_regen.png'}")


if __name__ == "__main__":
    main()
