# Panel g 讲解 — 165 种 Pt₃X₈ B-亚晶格化学计量扫描

**面板角色**: 把 "宿主选择" 推进到 "多元化学计量景观"；找出 Ga-rich 低焓窗口并验证实验配比落在其中。
**最新修订**: 2026-04-19（Element-reference 形成焓 + 7-slug 分类 + Origin 单表 + Ga-rich 绿带 / 少量 Zn 叙述）。
**完整脚注来源**: [legacy_FigG_Technical_Report.md](legacy_FigG_Technical_Report.md)

---

## 1. 这张图回答什么问题

**在固定 Pt 主骨架下，当 B 位的 Ga / In / Sn / Zn 化学计量发生变化时，哪一类组分组合在原子级能量上更稳定？**

这是 Fig.1 的 "多元化学计量景观" 面板。Panel c 证明了有序化可以发生；Panel g 证明哪个化学计量点是最优的。

---

## 2. 输入数据

| 来源 | 内容 |
|:---|:---|
| UMA-s-1p1 (Fairchem) | `165` 个 fixed-lattice single-point |
| 32-原子 L1₂ supercell | `a = 3.903 Å`，24 个 Pt 主骨架，8 个 B-site |
| B-site 化学计量枚举 | 所有满足 $n_{\text{Ga}} + n_{\text{In}} + n_{\text{Sn}} + n_{\text{Zn}} = 8$ 的整数组合，共 $\binom{11}{3} = 165$ 个点 |
| 每个化学计量点 | 一个 deterministic representative occupancy（不是随机），保证 provenance |
| 元素参考能 (2026-04-17 重建) | Pt `mp-126`; Ga `mp-142`; In 生成 `I4/mmm`; Sn 生成 diamond-α; Zn 生成 hcp。详见 [../04_reports/20260417_UMA_ElementReference_Strategy.md] |

---

## 3. 计算方法与公式

### 3.1 Element-referenced formation enthalpy（manuscript 图面使用）

$$
\boxed{\Delta H_f = \left( E_{alloy} - \sum_i x_i E_i^{elem} \right) \times 96.485 \text{ kJ/mol/atom}}
$$

**这是 2026-04-17 起采用的严格定义**，替代了之前的 `E_ref = -4.356 eV/atom` 历史校准。之前的 `Calibrated_Display_kJ_mol` 列作为遗留供参考，但 manuscript 图面 Y 轴 = `ElementRef_Hf_kJ_mol`。

### 3.2 7-类别 slug 分类（匹配 Origin 模板）

根据 8-site B-sublattice 上 (Ga, In, Sn, Zn) 的占用数 $(n_g, n_i, n_s, n_z)$ 分类:

| 非零元素数 | 额外条件 | Slug |
|:---:|:---|:---|
| 1 | 单一非 Pt 元素 | `Pure` |
| 2 | (Ga, Sn) | `GaSn` |
| 2 | (Ga, In) | `GaIn` |
| 2 | (In, Sn) | `InSn` |
| 2 | 含 Zn | `PartZn` |
| 3 | (Ga, In, Sn) | `GaInSn` |
| 3 | 含 Zn | `PartZn` |
| 4 | 全部 4 种 | `FiveElem` (Pt + 4 种 = 5 元素) |

### 3.3 Ga% binning

$$
\text{Ga\%} = \frac{n_{\text{Ga}}}{8} \times 100\% \in \{0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100\}
$$

### 3.4 Jitter (视觉用，不改变数据)

X 方向 uniform $[-0.35, +0.35]$，seed = 42（与历史 `Fig1h_Origin_AllPanels_PreciseNames_Hf.csv` 保持一致，方便 Origin 模板复用）。

---

## 4. 判断标准

- **Ga-rich favourable window**（绿色半透明矩形底纹）:
  - $\text{Ga\%} \in [62.5, 100]$
  - $\Delta H_f \in [-46, -34]$ kJ/mol
  - 填充色 `#2E7D32` @ 10% α，Z-order 最底
- **Target Ga fraction**（垂直虚线）:
  - $x = 68.5\%$ — AI-identified target composition，位于 favourable window 内；经 SI-Zn 的 down-selection (Panel g → CALPHAD liquidus → Galinstan eutectic accessibility) 得到
  - 色 `#D32F2F`, 1.3 pt dashed
  - 贴底标注 `Target 68.5 % Ga`
- **Bin mean ± 1σ overlay**（误差棒）:
  - 每个 Ga% bin 一个 open diamond (`#C62828`, 9 pt, white fill)
  - Error bar = $\pm 1\sigma$（population，ddof=0），cap 3 pt

---

## 5. 生成脚本与运行

| 文件 | 作用 |
|:---|:---|
| [script_FigG_ElementRef.py](script_FigG_ElementRef.py) | (2026-04-18) 枚举 165 个点 + UMA single-point + 元素参考能 → 生成 `data_FigG_165_ElementReferenced_Hf.csv` |
| [script_FigG_Origin_Ready.py](script_FigG_Origin_Ready.py) | (2026-04-19) 读取 165 行 ElementRef CSV → 重排成 Origin 单表 + 每 Ga% bin 的 Mean/Std overlay |

**Python 环境**:
- `ElementRef.py`: **hea_ai**（需要 fairchem UMA）
- `Origin_Ready.py`: **py312**（仅 pandas/numpy）

**运行**:
```bash
# 1) 首次生成（耗时 ~15 min）
python script_FigG_ElementRef.py

# 2) 纯重排到 Origin 格式（< 1 s）
python script_FigG_Origin_Ready.py
```

---

## 6. 输出数据表

| 文件 | 行 × 列 | 用途 |
|:---|:---|:---|
| [data_FigG_165_ElementReferenced_Hf.csv](data_FigG_165_ElementReferenced_Hf.csv) | 165 × ~16 | **源数据** (long format): Rank, Composition, Ga/In/Sn/Zn_count, Ga_pct, Energy_eV_atom, `ElementRef_Hf_kJ_mol`, Category, Optimal_Window |
| [data_FigG_Origin_Ready.csv](data_FigG_Origin_Ready.csv) | 35 × 99 | **Origin 直接导入**: 7 slug × 9 bin 的 XY 列 + 每 bin 的 Mean_X/Y/Std overlay |
| [data_FigG_GaPct_Statistics.csv](data_FigG_GaPct_Statistics.csv) | 9 行 | 每 Ga% bin 的 N, Mean, Std, Min, Max |
| [data_FigG_KeyPoint_Relax_Sensitivity.csv](data_FigG_KeyPoint_Relax_Sensitivity.csv) | 5 行 | 2026-04-17 cell-relaxation 敏感性关键点试算 |

**Origin Ready CSV 列结构**:

```
{GA}%_{CAT}_X, {GA}%_{CAT}_Y        ← 每个 bin 里出现的每个 category 一对
{GA}%_Mean_X, {GA}%_Mean_Y, {GA}%_Mean_Std    ← 每个 bin 一个 mean overlay
```

例: `0%_InSn_X, 0%_InSn_Y, 0%_PartZn_X, 0%_PartZn_Y, 0%_Pure_X, 0%_Pure_Y, 0%_Mean_X, 0%_Mean_Y, 0%_Mean_Std, 12%_GaSn_X, ...`

---

## 7. Origin 绘制要点

完整操作 + 符号样式总表见顶层 [99_Origin_Plotting_Guide_All_Panels.md § Panel g](../99_Origin_Plotting_Guide_All_Panels.md#panel-g).

每 Ga% bin 一层 layer（共 9 层），7 slug + Mean 共 8 种样式，**9 层共享同一套配色**:

| Slug | Marker | 颜色 | Size |
|:---|:---|:---|:---:|
| GaSn | Circle ● | `#1E88E5` | 6 pt |
| GaIn | Square ■ | `#43A047` | 6 pt |
| InSn | Diamond ◆ | `#66BB6A` | 6 pt |
| GaInSn | Up-Triangle ▲ | `#FB8C00` | 6 pt |
| FiveElem | Hexagon ⬢ | `#78909C` | 6 pt |
| PartZn | Down-Triangle ▼ | `#AB47BC` | 7 pt |
| Pure | Star ★ | `#111111` | 11 pt |
| Mean overlay | Open Diamond ◇ + error bar | `#C62828` | 9 pt |

**额外元素**:

- 绿色矩形底纹: $x \in [62.5, 100]$, $y \in [-46, -34]$
- 红色虚线 $x = 68.5\%$ + 底标 `Expt. 68.5 % Ga`
- 绿带内文字: `Ga-rich favourable window (62.5 – 100 % Ga)`

---

## 8. 结果小结

### 8.1 各 Ga% bin 的统计（element-referenced ΔH_f, kJ/mol）

| Ga% | N | Mean | Std | Min | Max |
|---:|---:|---:|---:|---:|---:|
| 0 | 45 | −24.82 | 3.17 | −31.26 | −19.27 |
| 12.5 | 36 | −27.21 | 2.74 | −32.62 | −23.10 |
| 25 | 28 | −29.55 | 2.32 | −33.99 | −26.33 |
| 37.5 | 21 | −31.87 | 1.94 | −35.46 | −29.15 |
| 50 | 15 | −34.19 | 1.57 | −37.04 | −32.00 |
| 62.5 | 10 | −36.52 | 1.24 | −38.69 | −34.84 |
| 75 | 6 | −38.83 | 0.94 | −40.29 | −37.70 |
| 87.5 | 3 | −41.10 | 0.65 | −41.82 | −40.53 |
| 100 | 1 | **−43.39** | — | −43.39 | −43.39 |

**单调下降**: Ga% 越高 → 平均 $\Delta H_f$ 越负。

### 8.2 全局最低点

**Pt₃Ga (Ga₈, 100% Ga) = −43.39 kJ/mol** 是 element-referenced 尺度下 165 个点的最低能量（原始值 `-43.384818`，2 位小数四舍五入）。

### 8.3 Ga-rich favourable window

- $\text{Ga\%} \in [62.5, 100]$, $\Delta H_f \lesssim -34$ kJ/mol
- 覆盖了 Pt₃Ga 全族 + 大多数 Ga-Sn / Ga-In / Ga-In-Sn 原型
- 实验目标配比 Ga = 68.5% 由本窗口配合 SI-Zn 的 CALPHAD 液相面约束与 Galinstan eutectic accessibility 共同 down-select 得到 (详见 SI-Zn)

### 8.4 少量 Zn 的角色

- Pure Pt₃Zn 的 $\Delta H_f$ 明显浅于 Pt₃Ga（~5 kJ/mol 差）
- PartZn 点（含 Zn 但非五元）在 0 – 62.5% Ga 各 bin 都**系统性偏上 ~5 kJ/mol**
- FiveElem 点（全 4 元）在相同 Ga% 下也比 Ga-Sn / Ga-In-Sn 浅 ~5 kJ/mol

**结论**: Zn 不深化 $\Delta H_f$。据此 AI 推断合成目标应尽量少加 Zn（仅作构型熵 / 液相温度调节微量添加），不应作为主元。manuscript 最终目标配比中 ~1 at.% Zn 的具体数值由 SI-Zn 的 CALPHAD 液相约束 + host-ranking 鲁棒性 confirmed（Panel b/d/f 对 Zn ∈ [0, 5 at.%] 不敏感）；多加 Zn 会用焓换熵，不推荐。

### 8.5 读者应看到

右半部（$\text{Ga\%} \ge 62.5\%$）绿色带内部密集堆着 Pure Pt₃Ga 黑星 + Ga-Sn / Ga-In / Ga-In-Sn 点；左半部（$\text{Ga\%} \le 37.5\%$）红色误差棒 diamond 明显走高，紫色 PartZn + 灰色 FiveElem 散布在 $\Delta H_f > -30$ kJ/mol 的区域。红色虚线穿过绿带左半侧 — AI-identified 目标配比落在安全窗口内（非最优点），为保证液相可及性做了 Ga% 下调的权衡（见 SI-Zn）。

---

## 9. 这张图不能单独说明什么

- **Panel g 是 fixed-cell representative screen**，不是完全弛豫的全局化学计量相图
- 2026-04-17 cell-relaxation 敏感性试算（见 [data_FigG_KeyPoint_Relax_Sensitivity.csv](data_FigG_KeyPoint_Relax_Sensitivity.csv)）:
  - `Ga8 (Pt3Ga)`: $-5.162 \to -5.178$ eV/atom（微变）
  - `Sn8 (Pt3Sn)`: $-5.142 \to -5.336$ eV/atom（**大幅下移 ~9 kJ/mol**）
  - `In8 (Pt3In)`: $-4.837 \to -4.999$ eV/atom（**大幅下移 ~16 kJ/mol**）
  - `Ga5Sn3`: 20 步内未完全收敛
- **方法学边界**: Sn / In 端元对 cell relaxation 非常敏感；Ga-rich 核心相对稳健。因此当前 Panel 的定位是 "fixed-cell representative stoichiometric screen"，全球化学计量基态景观**仍需**完整弛豫计算（超出当前工作范围）
- 它不能回答 "一个随机的 Ga₆InSn 是否稳定" — 每个化学计量点只试了一个 representative occupancy，而非所有排列

---

## 10. 参考

- 完整 Technical Report: [legacy_FigG_Technical_Report.md](legacy_FigG_Technical_Report.md)
- Element-reference 重建决策: `../04_reports/20260417_UMA_ElementReference_Strategy.md`
- 2026-04-18 Ga-rich 绿带 + "少量 Zn" 叙述: 本日 / 昨日会话交互记录
- Origin 符号总表: [99_Origin_Plotting_Guide_All_Panels.md § Panel g](../99_Origin_Plotting_Guide_All_Panels.md#panel-g)
