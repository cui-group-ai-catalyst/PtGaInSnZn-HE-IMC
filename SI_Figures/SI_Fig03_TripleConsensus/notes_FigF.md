# Panel f 讲解 — 15 种 M-Ga 二元的三方法共识

**面板角色**: 给 corrected Miedema 排序增加一个独立的原子级 (DFT/ML) 参照。
**最新修订**: 2026-04-18（Miedema **已从主图移除**；主图 = MP-DFT + UMA + CHGNet 三方法均值 + 误差棒）。
**完整脚注来源**: [legacy_FigF_Technical_Report.md](legacy_FigF_Technical_Report.md)

---

## 1. 这张图回答什么问题

**如果只看 `M-Ga` 二元参考体系，`MP-DFT` / `UMA` / `CHGNet` 三方法对每个 host 的形成焓排序是否一致？**

这一步的目的:

1. 给前面 panel 的 corrected Miedema 排序增加一个原子级独立参照
2. 证明 Miedema 的 ranking 在 Ga-系二元上和 DFT/ML 家族大体一致
3. 保留原子级数值而不是把不同量纲硬凑成一个虚假的 "均值 + 偏差"

---

## 2. 输入数据

| 来源 | 内容 |
|:---|:---|
| Materials Project (本地 CIF 归档) | 15 个 `M-Ga` 二元相 + 元素参考态；见 [data_FigF_MP_DFT_References.csv](data_FigF_MP_DFT_References.csv) |
| `02_data/20260415_FigF_LocalStructures/binary_manifest.csv` | 二元相 CIF 清单 + MP-ID |
| `02_data/20260415_FigF_LocalStructures/element_manifest.csv` | 元素参考态 CIF 清单 |
| UMA-s-1p1 (Fairchem) | 对上述本地 CIF 做 single-point，得 $\Delta H_f^{\text{UMA}}$ |
| CHGNet | 同上做 single-point，得 $\Delta H_f^{\text{CHGNet}}$ |

**关键**: 当前脚本**不依赖实时 Materials Project API**。所有结构和 MP-DFT 能量都已本地化，保证审查可离线复现。

---

## 3. 计算方法与公式

### 3.1 MP-DFT 参考值

直接从 Materials Project 拉取的 formation energy per atom，单位 eV/atom → 乘 `96.485` 得 kJ/mol。

### 3.2 UMA / CHGNet 形成焓

对每个二元相:

$$
E_{alloy} \xrightarrow{\text{UMA or CHGNet}} \text{potential energy}
$$

元素参考态加权:

$$
E_{ref} = x_M E_M + x_{Ga} E_{Ga}
$$

形成焓:

$$
\Delta H_f = (E_{alloy} - E_{ref}) \times 96.485 \text{ kJ/mol/atom}
$$

### 3.3 取绝对值做条形图

$$
|\Delta H_f|^{\text{(method)}} \quad \text{for method} \in \{\text{MP-DFT, UMA, CHGNet}\}
$$

### 3.4 三方法共识柱 + 误差棒

$$
\text{Bar\_Mean\_absHf} = \frac{1}{3} \sum_{m} |\Delta H_f|^{(m)}
$$

$$
\text{Bar\_Std\_absHf} = \text{std}\left( \{ |\Delta H_f|^{(m)} \},\ \text{ddof}=1 \right)
$$

### 3.5 Consensus rank

$$
\text{Consensus\_Rank\_Mean} = \frac{1}{3} \left( \text{Rank}_{\text{MP-DFT}} + \text{Rank}_{\text{UMA}} + \text{Rank}_{\text{CHGNet}} \right)
$$

host 按此值升序排列，rank 1 = 最强 binder 在最左。

---

## 4. 判断标准

- 主图**不显示**个别方法的点（仅 bar + error bar），避免读者误把 ranking 差异放大
- 图中只保留两类 host 视觉编码:
  - **Pt** = `#C62828` 红色条（焦点 host，跨 panel 一致）
  - 其他 = `#455A64` 深灰条
  - **Size_Pass = False** 的 host 在条上加斜线 hatch `///`（标记不满足 Panel e 几何约束）
- Miedema 值已被主图**剔除**，但仍保留在源数据 [data_FigF_TripleConsensus_Data.csv](data_FigF_TripleConsensus_Data.csv) 的 `Miedema_abs_dH_kJ_mol_atom` 列中，供 SI 或审稿人质疑时随时调出。理由: **Miedema 对非等摩尔组分 (Ga5Ir3, Ga3Ni2, Ga3Fe, Ga4Cu9, ...) 量级偏小 ~2× 倍**，因为正规溶液 $4 \Delta H_{AB} c(1-c)$ 在端元附近塌缩到零。

---

## 5. 生成脚本与运行

| 文件 | 作用 |
|:---|:---|
| [script_FigF_TripleConsensus.py](script_FigF_TripleConsensus.py) | 完整 pipeline: 读 CIF → UMA + CHGNet 单点 → 导出 Data/Origin/Summary CSVs |
| [script_FigF_TripleConsensus_Replot.py](script_FigF_TripleConsensus_Replot.py) | 纯重绘脚本: 读 `Data.csv` → 生成最新版 Origin 表和 preview PNG；**不依赖 fairchem**，可在普通 py312 环境跑 |

**Python 环境**:
- 完整 pipeline (`TripleConsensus.py`): **hea_ai** 环境（需要 fairchem + chgnet）
- 重绘 (`TripleConsensus_Replot.py`): **py312**（仅 pandas/matplotlib）

**运行**:
```bash
# 首次生成原始数据（耗时 ~5 min）
python script_FigF_TripleConsensus.py

# 样式调整 / 重绘（< 5 s）
python script_FigF_TripleConsensus_Replot.py
```

---

## 6. 输出数据表

| 文件 | 行 × 列 | 用途 |
|:---|:---|:---|
| [data_FigF_TripleConsensus_Data.csv](data_FigF_TripleConsensus_Data.csv) | 15 × 多列 | 所有三方法原始 signed/abs ΔH_f + Miedema（供 SI 审查） |
| [data_FigF_TripleConsensus_Origin.csv](data_FigF_TripleConsensus_Origin.csv) | 15 行 | **Origin 直接导入**: `Host`, `Formula`, `Bar_Mean_absHf`, `Bar_Std_absHf`, 三方法 signed/abs + Rank, `Consensus_Rank_Mean`, `Size_Pass` |
| [data_FigF_TripleConsensus_Summary.csv](data_FigF_TripleConsensus_Summary.csv) | 几行 | Spearman ρ 汇总，供图内文字框 |
| [data_FigF_MP_DFT_References.csv](data_FigF_MP_DFT_References.csv) | 15 行 | 每个 host 的 MP-ID、结构式、MP-DFT formation energy |

---

## 7. Origin 绘制要点

完整操作见顶层 [99_Origin_Plotting_Guide_All_Panels.md § Panel f](../99_Origin_Plotting_Guide_All_Panels.md#panel-f).

核心:

1. 导入 `data_FigF_TripleConsensus_Origin.csv`
2. `Plot → Column` 主 Y = `Bar_Mean_absHf_kJ_mol_atom`
3. 误差棒 = `Bar_Std_absHf_kJ_mol_atom`，upper+lower, cap 3 pt
4. Bar 颜色: Pt 行 = `#C62828`；其余 = `#455A64`
5. `Size_Pass = False` 的 bar 套 diagonal hatch `///`
6. X tick label 两行: `{Host}\n{Formula}`
7. Y 轴 $|\Delta H_f|$ (kJ mol⁻¹ atom⁻¹)，range `[0, 80]`
8. 左上文字框: Spearman ρ 汇总 (MP-DFT vs UMA, MP-DFT vs CHGNet, UMA vs CHGNet)

---

## 8. 结果小结

**Spearman rank correlation** (from `_Summary.csv`):

- MP-DFT vs UMA (size-pass subset): **+0.951** → 非常一致
- MP-DFT vs CHGNet (all 15): 中等偏强
- UMA vs CHGNet (all 15): 中等偏强

**Consensus rank top-6**:

`Rh #1` → `La #2` → `Pd #3` → `Y #4` → `Ce #5` → **`Pt #6`**

**读者应看到**: Rh / La / Pd / Y / Ce / **Pt** 聚在左侧 $|\Delta H_f| \approx 60 – 70$ kJ/mol/atom；Cu 和 Re 在最右 ~10 kJ/mol/atom。La 和 Ce 的误差棒明显更大，反映 DFT/ML 对稀土 gallide 的分歧 — 一个有用的 caveat。**Pt 位于 "strong-binder" 邻域**，与 Panel b/d/e 的结论一致。

**结论**: **MP-DFT / UMA / CHGNet 三方法在 size-compatible host 子集中呈现强一致性，足以作为 Ga-系二元参考验证**。Pt 的 Ga 亲和在这一独立原子级尺度上仍然稳居高位。

---

## 9. 这张图不能单独说明什么

- Panel f 是**二元** `M-Ga` 参考，不是多元 `Ga-In-Sn-Zn` 终态计算
- 它不能单独推出最终多元有序相的化学计量最优点 → **Panel g**
- Miedema 从主图移除**不代表**该方法错，只说明它对非等摩尔组分量级不可靠 (~2× 偏小)。它仍保留在源数据中
- 它不替代 Panel c 的 ordered-vs-disordered 比较

因此 Panel f 在证据链中的角色是 "验证"，不是 "终判"。

---

## 10. 参考

- 完整 Technical Report: [legacy_FigF_Technical_Report.md](legacy_FigF_Technical_Report.md)
- 2026-04-18 "Miedema 从主图移除" 决策记录: 本日会话交互
- 原始 MP 结构清单: `../../02_data/20260415_FigF_LocalStructures/` (外部)
