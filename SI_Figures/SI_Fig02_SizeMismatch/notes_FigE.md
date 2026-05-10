# Panel e 讲解 — 原子尺寸失配 × 化学驱动力的二维筛选

**面板角色**: 把化学吸引 (Panel b/d) 和几何匹配 (Panel e 本身) 放在同一张图里讨论。
**生成日期**: 2026-04-15（本版沿用；无后续方法学修改）。
**完整脚注来源**: [legacy_FigE_Technical_Report.md](legacy_FigE_Technical_Report.md)

---

## 1. 这张图回答什么问题

**即使某个宿主在热力学和界面润湿上都表现良好，它在固相中是否仍然具备足够的几何兼容性，不会因为尺寸失配破坏晶格稳定？**

Panel e 扮演 "几何约束" 角色 — 筛掉那些虽然化学吸引强，但尺寸过大或过小、难以稳定进入目标晶格的宿主（典型如稀土 La / Ce / Y）。

---

## 2. 输入数据

| 来源 | 内容 |
|:---|:---|
| `01_scripts/data_periodic_table.py` | 统一参数: `r` (原子半径), `Phi`, `n_ws` |
| 脚本内固定 liquid cocktail | Ga 0.65 / In 0.20 / Sn 0.10 / Zn 0.05 |

使用固定 cocktail 的理由同 Panel b/d，保证跨 panel 解释不矛盾。

---

## 3. 计算方法与公式

### 3.1 等效液体半径

$$
r_{cocktail} = \sum_k w_k r_k
$$

### 3.2 尺寸失配

$$
\delta = \left| \frac{r_{cocktail} - r_{host}}{r_{host}} \right| \times 100\%
$$

### 3.3 化学驱动力（Corrected Miedema）

$$
\text{Enthalpy\_Drive} = \sum_k w_k \left[ -P (\Delta \Phi_k)^2 + Q (\Delta n_{ws,k}^{1/3})^2 \right]
$$

**注意**: 这里的 `Enthalpy_Drive` 是 corrected Miedema 驱动力，**不是** 第一性原理的 "真实形成焓"。避免混用术语。

---

## 4. 判断标准（启发式高亮规则，非物理相界）

$$
\delta < 15\% \qquad \text{AND} \qquad \text{Enthalpy\_Drive} < -20 \text{ kJ/mol}
$$

这两条线**不是**严格物理相边界，而是 manuscript 讨论中的启发式高亮规则。用于把 "同时满足几何和化学" 的 host 圈出来。

---

## 5. 生成脚本与运行

| 文件 | 作用 |
|:---|:---|
| [script_FigE_compute_and_plot.py](script_FigE_compute_and_plot.py) | 完整 compute + plot 管线 |

**Python 环境**: py312（仅需 `numpy` / `pandas` / `matplotlib`）
**运行**: `python script_FigE_compute_and_plot.py`
**运行时长**: < 1 s

---

## 6. 输出数据表

| 文件 | 行 × 列 | 用途 |
|:---|:---|:---|
| [data_FigE_Resistance_Ranked.csv](data_FigE_Resistance_Ranked.csv) | 28 | 主数据表，含 `Mismatch_Percent` + `Enthalpy_Drive` |
| [data_FigE_True_ThreeWay.csv](data_FigE_True_ThreeWay.csv) | 扩展三方对照版本 | 与预览图 `preview_FigE_ThreeWay.png` 对应 |

**主表关键列**:

- `Host` — 元素符号
- `r_host` — 原子半径 (Å)
- `Mismatch_Percent` — $\delta$ (%)
- `Enthalpy_Drive` — corrected Miedema 驱动力 (kJ/mol)
- `High_Mismatch` / `Low_Enthalpy` — 布尔高亮标记

---

## 7. Origin 绘制要点

完整操作见顶层 [99_Origin_Plotting_Guide_All_Panels.md § Panel e](../99_Origin_Plotting_Guide_All_Panels.md#panel-e).

核心:

- X 轴: `Mismatch_Percent` (%)，log 或 linear 均可
- Y 轴: `Enthalpy_Drive` (kJ/mol)，**更负在下**
- 两条高亮线: 垂直 $x = 15$、水平 $y = -20$
- 左下 "stable zone" 绿色半透明填充 + 内标文字 "Size + chemistry favourable"
- Pt 用 `#C62828` 大圆 + 黑边，其他 host 灰色小圆
- 稀土 (La, Ce, Y) 额外标 label，标注其为何被筛除

---

## 8. 结果小结

**Pt 的关键值**:

- `Mismatch_Percent = 2.05%`
- `Enthalpy_Drive = −35.38 kJ/mol`

**按 Mismatch 从小到大排序 Top-8**:

| Host | Mismatch (%) | Enthalpy_Drive (kJ/mol) |
|:---|---:|---:|
| Au | 1.49 | −16.65 |
| Ag | 1.49 | −1.23 |
| Mo | 2.05 | −4.79 |
| **Pt** | **2.05** | **−35.38** |
| W | 2.05 | −7.54 |
| Nb | 2.84 | −0.09 |
| Ta | 2.84 | −0.03 |
| Ti | 3.50 | −2.50 |

**满足两条高亮规则**（$\delta < 15\%$ AND $\text{Enthalpy\_Drive} < -20$ kJ/mol）**的 host 共 7 个**:

`Pt`, `Ir`, `Pd`, `Rh`, `Ru`, `Re`, `Os`

**结论**: **Pt 同时占据 "极低尺寸失配" + "最强 corrected Miedema 驱动力" 两个优势**，是最稳健的结构候选宿主。Au 虽然尺寸好但驱动力弱；Mo/W 虽然尺寸好但驱动力弱；稀土虽然驱动力强但尺寸差。只有 Pt/Ir/Pd/Rh/Ru/Re/Os 这 7 个 host 同时满足。

---

## 9. 这张图不能单独说明什么

- 两条 `15%` / `-20 kJ/mol` 线是**启发式高亮**，不是严格物理相界
- 只能说明几何匹配 + corrected Miedema 驱动力的联合趋势，**不能**说明:
  - 原子级有序构型的基态能量 → **Panel c**
  - 多元化学计量景观的全局最低点 → **Panel g**
  - 动力学路径或相变速率（不在本文范围）
- `Enthalpy_Drive` 不应被写成 "真实 Enthalpy" — 它是 corrected Miedema 的半经验驱动力

---

## 10. 参考

- 完整 Technical Report: [legacy_FigE_Technical_Report.md](legacy_FigE_Technical_Report.md)
- 与 Panel c / g 联合解释: 见 [../03_Panel_c_OrderedVsDisordered/notes_FigC.md] 和 [../07_Panel_g_GaSweep/notes_FigG.md]
