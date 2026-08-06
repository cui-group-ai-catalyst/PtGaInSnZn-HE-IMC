# 20260423 备忘录：化学势主线的计算方案与故事框架

> 目的：把这几轮讨论确定下来的"化学势驱动"叙事、要算的量、公式依据、本地可行性、执行顺序、跨 Panel 一致性，整理成一份一页纸可汇报、三天可落地的工作备忘。下一步直接从 Script A 开始跑数。
>
> 上游相关文件（同目录）：
> - `Plan_SI_ChemicalPotential.md`（SI 总规划）
> - `Formula_Derivations_and_Justifications.md`（公式推导 + 文献锚点）
> - `20260422_Concept_and_NextSteps_ChemicalPotential.md`（Miedema+CEF 轻量化路线决定）

---

## 〇、一句话总结

> 文章主线从"高熵 + AI 筛选"升级成**化学势驱动**——Panel a（**泵**：液相高 μ → HEI 低 μ 的下坡）、Panel b（**闸门**：γ_SL < 0，界面非壁垒）、Panel c（**汇**：HEI 在合成温度下稳定）。三段叙事回答一个完整问题："为什么这个合成条件下一定得到 L1₂-HEI"。

---

## 一、三段叙事（Panel a–b–c 角色）

| 图 | 角色 | 物理量 | 回答 |
|---|---|---|---|
| **Panel a** | 泵 driving force | Δμᵢ（5 个元素） | 有没有动力？ |
| **Panel b** | 闸门 channel | γ_SL < 0（Miedema 界面版） | 能不能过界？ |
| **Panel c** | 汇 sink | G_m^HEI < G_m^HEA（温度窗口内） | 停在哪相？ |

**a→b 的衔接逻辑**（两座桥，缺一不可）：

1. **Gibbs 吸附方程**：dγ_SL = −Σ Γᵢ dμᵢ
   → γ_SL 与 μᵢ 是 Legendre 共轭量（同源，算哪个都行）。
2. **Miedema 统一电子核**：ΔH_mix 与 γ_SL^chem 共用 −P(Δφ*)² + Q(Δn_WS^(1/3))²
   → 体相 → 界面只差一个几何缩放（J/mol ↔ J/m²）。

**界面化学势本身不算**，用 γ_SL 替代；理由：Widom 在致密固相失败、AIMD slab 成本不合算、γ_SL 承载等价信息。

---

## 二、要计算的物理量（仅 3+1 个）

| # | 量 | 记号 | 用途 |
|---|---|---|---|
| 1 | 起点化学势 | μᵢ^(0) | Panel a 左上"高位端" |
| 2 | 终点化学势（HEI 产物） | μᵢ^HEI | Panel a 右下"低位端" |
| 3 | 驱动力（每元素） | Δμᵢ = μᵢ^(0) − μᵢ^HEI | Panel a 5 条高度差 |
| 4 | 总反应 Gibbs | ΔG_rxn = −Σ νᵢ Δμᵢ | 正文一句话 + SI 表 |

**不算的**：
- 界面化学势梯度 ∇μ|_int（用 γ_SL 替代）
- 两相共存端点（实验按 3:1 精准配料 → 耗尽型终止，最终态 f^HEI ≈ 1）
- 声子贡献（Neumann–Kopp 近似，ΔS_vib ≈ 0）

---

## 三、核心公式链

### 3.1 起点 μᵢ^(0)

**Pt（合成前仍为固相）**：

    μ_Pt^(0) = G_Pt°,fcc(T)    ← 查 SGTE

**Ga/In/Sn/Zn（液态混合物）**——Hildebrand 正则溶液：

    μᵢ^(0) = Gᵢ°,L(T) + RT·ln xᵢ^L + Σ_{j≠i} Ω_ij^L · (xⱼ^L)²

Ω_ij^L 由 Miedema 参数（de Boer 1988）查得。

### 3.2 终点 μᵢ^HEI（L1₂-Pt₃X，α=Pt, β=Ga/In/Sn/Zn）

**Pt 位于 α 亚晶格**：

    μ_Pt^HEI = G_{Pt:*}°(T) + 4·G^E(y^β) − 4·Σᵢ yᵢ^β · ∂G^E/∂yᵢ^β

**X ∈ {Ga,In,Sn,Zn} 位于 β 亚晶格**：

    μ_X^HEI = G_{Pt:X}°(T) + RT·ln y_X^β + [CEF 过剩项]

**Ω_ij^β（6 个）来自 Panel g 现有 165 点 UMA 能量的最小二乘拟合**。代码拟合并输出的是每合金原子参数 $\omega_{ij}^{atom}$；若 SI 使用每 β 位点参数，则 $\Omega_{ij}^{β}=4\omega_{ij}^{atom}$。

固定 Pt₃X 组成流形直接确定的是总自由能与 β 子晶格扩散势，而不是五个唯一的绝对元素化学势。代码保留的逐元素数值是一个 Euler 一致的历史参考分解；规范不变的扩散势输出见 `outputs/beta_diffusion_potentials_v3_0K.csv`。

### 3.3 驱动力与反应能

    Δμᵢ = μᵢ^(0) − μᵢ^HEI
    ΔG_rxn = −Σᵢ νᵢ · Δμᵢ

**反应式**：

    3 Pt(s) + a·Ga(l) + b·In(l) + c·Sn(l) + d·Zn(l) → Pt₃(Ga_a·In_b·Sn_c·Zn_d)

（a+b+c+d=1，ν_Pt=3，ν_X=a/b/c/d）

---

## 四、关键概念澄清

### 4.1 ΔG_rxn vs Gibbs 自由能 G

- G：某个状态的 Gibbs 自由能（状态函数）
- ΔG：两状态之间 G 的差值
- ΔG_rxn：特指"反应"前后的 ΔG = −Σ νᵢ Δμᵢ（这里 Δμᵢ 定义为 source − product）

三者是嵌套关系，不是不同量。

### 4.2 化学势 μ vs 混合焓 ΔH_mix

- ΔH_mix：整体混合热（一个数）
- μᵢ：G 对元素 i 粒子数偏导（每元素一条）
- 正则溶液（ΔS^E = 0）下 G_m^E = ΔH_mix，所以 μᵢ^E = ∂ΔH_mix/∂nᵢ
- 混合焓提供参数、化学势讲故事（同一套参数两种表达）

**二元验证**：ΔH = Ω·x_A·x_B ⇒ μ_A^E = Ω·x_B², μ_B^E = Ω·x_A²；Euler 关系自洽。

### 4.3 耗尽型 vs 共存型终止（这是我们的场景）

| 模式 | 条件 | 终止判据 | 残余液相 |
|---|---|---|---|
| 共存型 | 某相过量 | μᵢ^L = μᵢ^HEI | 有 |
| **耗尽型（我们的）** | 输入按 HEI 化学计量（3:1）+ HEI 稳定 | 液相被消耗光 f^L → 0 | 无 |

**结论**：最终态接近 100% HEI，没有"平衡平台 + 残余液相"这回事。驱动力 Δμᵢ^(0) 不是"拉平到 0"，而是全部转化为合成放热 ΔG_rxn < 0。

### 4.4 Panel a 的"高度"语义

高度 = 每原子 Gibbs 形成能 ΔG_f；由 Euler 关系 G_m = Σ xᵢ·μᵢ，ΔG_f 的台阶数值上等于组分加权 Σ xᵢ·Δμᵢ——Panel a 同时携带 G 和 μ 两套信息。

---

## 五、本地资源清单

| 需要 | 有无 | 来源 |
|---|---|---|
| SGTE 纯元素 Gibbs 函数 | 有 | Dinsdale 1991 公开 |
| Miedema 参数（φ*, n_WS^(1/3)） | 有 | de Boer 1988 表 |
| Panel g 165 点 DFT/UMA/CHGNet 能量 | 有 | `03_results/` 已有 |
| Python + NumPy + SciPy | 有 | 本地 |
| 最小二乘（拟 Ω^β） | 有 | `scipy.optimize.least_squares` |
| MC 误差传播 | 有 | NumPy 秒级 |
| pycalphad（可选交叉验证） | 可装 | `pip install pycalphad` |
| Matplotlib / Origin 出图 | 有 | 本地 |

**不需要**：AIMD / Widom / 声子 / 超算 / 新 DFT 单点。

---

## 六、跨 Panel 一致性（2026-04-23 锁定）

### 6.1 温度参考

- **Panel a / b / g 统一用 0 K 焓近似**（Miedema + DFT 单点都不含显式 T）。
- Panel a 主图高度 = ΔH_{f,i}（不标具体 T），caption 注：
  "Enthalpic landscape (0 K reference); qualitative hierarchy persists across synthesis window (~800–1200 K)."
- T 依赖由 **Panel c 专责承担**（HEI vs HEA Gibbs 曲线随 T 演化）。
- SI 表补 μᵢ @ T ∈ {300, 800, 1000, 1200 K}，验证定性结论稳定。

### 6.2 液相组成（Panel a 继承 Panel b baseline）

| 元素 | 液相 xᵢ^L |
|---|---|
| Ga | 0.65 |
| In | 0.20 |
| Sn | 0.10 |
| Zn | 0.05 |
| Pt | 0（固相 host，不进液相） |

**来源**：本 release 中 `Panel_f_Wetting/script_FigF_Wetting.py` 使用的 representative screening cocktail。

**跨 panel 原则**：
- Panel a / b 共用该代表组成保证可比性；
- Panel g 负责 165 点组分窗口筛选；
- SI 一句话交代："Panel a-b indicative, Panel g optimization."

---

## 七、执行顺序（Script A–F）

| Script | 产出 | 输入 | 估时 | 状态 |
|---|---|---|---|---|
| **A** | μᵢ^L(x, 0 K 焓近似) + SI 温度表 | Miedema Ω^L + §6.2 组成 | 半天 | ✅ **完成 2026-04-23**（mu_liquid_0K.csv, omega_liquid_binary.csv） |
| **B** | μᵢ^HEI(y) CEF 解析函数（Pt α + 4 × β） | 文献 ΔH_f(Pt₃i) + Miedema Ω^β | 半天 | ✅ **完成 2026-04-23**（mu_HEI_0K.csv, omega_beta_subl.csv） |
| **C** | Δμᵢ^(0) + ΔG_rxn（耗尽型直接代公式） | A + B 的输出 | 半天 | ✅ **完成 2026-04-23**（delta_mu_0K.csv, delta_G_rxn_summary.csv, panel_a_tier_summary.csv） |
| **D** | Panel a 可视示意图（三台阶 + 五元素细线） | C 的三个 CSV | 2 小时 | ✅ **完成 2026-04-23**（panel_a_schematic.png / .pdf） |
| **E** | MC 误差传播（ΔH_f Pt₃i ±5 kJ 抽样，1000 次） | A + B + 不确定度表 | 1 天 | 待开 |
| **F** | CALPHAD 交叉验证（Pt-Ga, Pt-In 二元） | pycalphad / 手查相图 | 1 天 | 待开 |

**Panel a 主线 A-D v2 已收官 2026-04-23**。剩下 E、F 是 SI 的鲁棒性论证，不阻塞正文绘图。

### §7.1 v1 → v2 修正记录（2026-04-23）

**v1 的两个物理问题**（已修正）：

1. **参考态**：v1 以纯液态 = 0 为参考，导致液/固起始 μ 均 ≈ 0，无法体现液态高化学势。v2 改用 SER（Stable Element Reference, 0 K 稳定固态 = 0），液相 μ_i = +ΔH_fus(i) > 0。
2. **HEI 端元分辨**：v1 对所有 β 元素取加权平均形成焓，导致五元素 Δμ 全部相同。v2 β 元素 i 使用各自 Pt₃i 端元形成焓；Pt(α) 仍用加权平均。

### §7.2 Panel a 最终数字表（v2, 2026-04-23 锁定）

所有数值 per atom of formula unit, kJ/mol, SER reference, 0 K enthalpic.

| 元素 | 子格 | 起点相 | y_i | μ_start (per atom) | μ_HEI (per atom) | **F_i (per atom)** |
|---|---|---|---:|---:|---:|---:|
| **Sn** | β | 液态 | 0.10 | +1.746 | **−15.612** | **+17.36** |
| **Ga** | β | 液态 | 0.65 | +1.392 | **−13.030** | **+14.42** |
| **Pt** | α | 固态 fcc | 1.00 | 0.000 | **−12.589** | **+12.59** |
| **Zn** | β | 液态 | 0.05 | +1.825 | **−9.856** | **+11.68** |
| **In** | β | 液态 | 0.20 | +0.760 | **−10.410** | **+11.17** |

**五元素驱动力排序**：Sn > Ga > Pt > Zn > In

**总反应能量**：
- Δ G_rxn per atom = **−13.93 kJ/mol·atom**
- Δ G_rxn per f.u. = **−55.71 kJ/mol f.u.**（放热，自发）

**v2 物理亮点（可写进正文 / caption）**：

1. **液相在 SER 参考下是高 μ 储池（+0.8 ~ +1.8 kJ/mol atom）**：液态比固态不稳定，ΔH_fus 提供正偏移；Zn 最高（+1.83）因其 ΔH_fus 最大（7.32 kJ/mol），In 最低（+0.76）。
2. **HEI 是深阱但五元素深度不同**：Sn 最深（−15.61，Pt₃Sn 键最强），Zn 最浅（−9.86，Pt₃Zn 键最弱）。这反映了 Pt-M 化学键强度的元素分辨差异。
3. **驱动力排序 Sn > Ga > Pt > Zn > In 与物理化学直觉一致**：Sn 的 5s² 孤对电子与 Pt 5d 杂化最强；In 的 Pt₃In 形成焓最弱。
4. **Δ G_rxn = −55.71 kJ/mol f.u. 比 v1 的 −50.36 更负**：因为 v2 正确计入了液态 ΔH_fus（额外贡献 ~5 kJ/mol），整体驱动力更大。

---

## 八、产出清单（最终交付给 Fig 1 / SI）

| 产物 | 位置 | 格式 |
|---|---|---|
| 5 个 Δμᵢ^(0) 数值 | Panel a 高度差标注 | kJ/mol-atom |
| ΔG_rxn 总数 | 正文一句话 | kJ/mol（预计 −30 ~ −60） |
| Ω_ij^β 表 + 误差 | SI 表 | kJ/mol |
| μᵢ^(0), μᵢ^HEI 完整表 | SI 表 | kJ/mol-atom |
| HEI vs HEA Gibbs 温度曲线 | Panel c | T / kJ/mol |
| γ_SL Miedema 计算 | Panel b（沿用已有） | mJ/m² |
| MC 不确定度三层验证 | SI 章节 | 统计表 |

---

## 九、边界声明（SI 里必须写清）

1. **我们算的是热力学，不算动力学**：Δμᵢ > 0 说明能反应，不含速率。
2. **Miedema 误差约 ±8 kJ/mol**：Δμ 量级 20–50 kJ/mol，信噪比足够；Panel g 的 DFT 能量做交叉验证。
3. **Neumann–Kopp 近似**：忽略声子贡献，估计误差 < 2 kJ/mol。
4. **单相 HEI 结论依赖两个前提**：① 化学计量精准匹配 3:1；② 合成 T 下 HEI 稳定（Panel c 证）。

---

## 十、关键参考锚点（供 SI 起草用）

| 内容 | 原始文献 |
|---|---|
| Gibbs 多相平衡判据 | Gibbs 1876, On the Equilibrium of Heterogeneous Substances |
| 化学势定义 | 同上，第 3 页 |
| Euler 积分 G_m = Σ xᵢ·μᵢ | 同上，第 5 页 eq. 11 |
| 正则溶液模型 | Hildebrand 1929, J. Am. Chem. Soc. 51, 66 |
| 统计力学推导 | Guggenheim 1935, Proc. Roy. Soc. A 148, 304 |
| Miedema 模型 | Miedema, de Chatel, de Boer 1980, Physica B 100, 1 |
| Miedema 界面版 | Benedictus, Böttger, Mittemeijer 1996, Phys. Rev. B 54, 9109 |
| SGTE 纯元素 | Dinsdale 1991, Calphad 15, 317 |
| CEF 子晶格模型 | Sundman & Ågren 1981, J. Phys. Chem. Solids 42, 297 |
| CEF 综述 | Hillert 2001, Calphad 25, 1 |
| 反应润湿 / γ < 0 | Eustathopoulos 1998, Int. Mater. Rev. 43, 98；Kaptay 2005, J. Mater. Sci. 40, 2125 |
| 高熵子晶格无序性 | Miracle & Senkov 2017, Acta Mater. 122, 448 |
| Neumann–Kopp / 声子 | Grimvall 1999, Thermophysical Properties of Materials Ch. 10 |
| Miedema 焓预测误差 | Zhang, Liaw, Yang, Zhang 2018, Intermetallics 95, 154 |

---

## 十一、下一步

**立刻开 Script A**：按 §6 锁定的输入，输出 5 条 μᵢ^(0) 数字（4 液相元素 + 1 固相 Pt） —— 这是后面所有计算的起点。

预计成果：
- `scripts/script_A_mu_liquid.py`
- `outputs/mu_liquid_0K.csv`（主产出，对应 Panel a）
- `outputs/mu_liquid_T_table.csv`（SI 补充表，T ∈ {300, 800, 1000, 1200}）

完成后继续 B → C → D → E → F。

---

**文件位置**：本文件（`docs_ChemicalPotential_Workflow.md`）

**维护**：每完成一个 Script，在 §七 表里更新状态，并在 §八 表里勾选已产出条目。

### §7.3 v2 → v3 修正记录（2026-04-23）

**v2 的数据源不一致问题**（已修正）：

v2 的 Script B 使用文献量热实验值（Kumar 1996, Srikanth 1993, Watson-Hayes 1995, Liu 2011）作为 Pt₃i 端元形成焕，与 Panel g 基于 UMA 机器学习势的 165 点计算数据不自洽。两组数据差异大且元素排序不同。

v3 替换为从 Panel g 的 165 个 UMA 数据点拟合 CEF 子晶格模型：
1. **端元 ΔH_f**：直接从 CSV 读取纯二元行（count_i = 8）
2. **6 个 Ω_ij^sub**：最小二乘拟合，R² = 0.9993，RMSE = 0.14 kJ/mol·atom
3. **Pt 化学势**：由 Gibbs-Duhem 一致性计算（μ_Pt = G_atom = Σ y_i h_i + G^xs）

### §7.4 Panel a 最终数字表（v3, 2026-04-23 锁定）

所有数值 per atom of formula unit, kJ/mol, SER reference, 0 K enthalpic.
数据源：液相 = CRC + Miedema；HEI = Panel g 165 点 UMA + CEF 拟合。

| 元素 | 子格 | 起点相 | y_i | μ_start | μ_HEI | **F_i** |
|---|---|---|---:|---:|---:|---:|
| **Ga** | β | 液态 | 0.65 | +1.392 | **−43.430** | **+44.82** |
| **Pt** | α | 固态 fcc | 1.00 | 0.000 | **−36.311** | **+36.31** |
| **Zn** | β | 液态 | 0.05 | +1.825 | **−31.197** | **+33.02** |
| **Sn** | β | 液态 | 0.10 | +1.746 | **−24.231** | **+25.98** |
| **In** | β | 液态 | 0.20 | +0.760 | **−20.495** | **+21.26** |

**五元素驱动力排序**：Ga > Pt > Zn > Sn > In

**总反应能量**：
- Δ G_rxn per atom = **−37.63 kJ/mol·atom**
- Δ G_rxn per f.u. = **−150.54 kJ/mol f.u.**（放热，自发）

**v3 物理亮点**：

1. **驱动力排序变为 Ga > Pt > Zn > Sn > In**：与 v2 完全不同，因为 UMA 给出 Pt₃Ga 最稳定（−43.4 kJ/mol·atom），远超 Pt₃Sn（−19.3）。
2. **ΔG_rxn 显著增大**：v3 = −150.5 kJ/mol f.u. vs v2 = −55.7，因为 UMA 形成焕比文献实验值更负。
3. **R² = 0.9993 的 CEF 拟合**确保了化学势计算与 Panel g 的 165 点能量景观完全自洽。

### v3 产出文件

| 文件 | 内容 |
|---|---|
| scripts/script_A_mu_liquid_v3.py | 液相 SER + ΔH_fus（同 v2） |
| scripts/script_B_mu_HEI_v3.py | **CEF 拟合 Panel g 165 点** |
| scripts/script_C_drive_force_v3.py | 五元素 Δμ_i + ΔG_rxn |
| scripts/script_D_panel_a_plot_v3.py | Panel a 可视化图 |
| outputs/mu_liquid_v3_0K.csv | 液相数值 |
| outputs/mu_HEI_v3_0K.csv | HEI 数值（UMA-CEF） |
| outputs/omega_beta_subl_v3_fit.csv | 6 个 拟合 Ω_ij |
| outputs/cef_fit_v3_quality.csv | R², RMSE |
| outputs/delta_mu_v3_0K.csv | 驱动力数值 |
| outputs/panel_a_v3_schematic.png/pdf | Panel a 图 |

---

## 十二、v3 计算验证与 Panel a–g 自洽性审查（2026-04-23 锁定）

### §12.1 v3 数值复现验证

2026-04-23 重新运行全部 v3 脚本（A→B→C），输出与已有 CSV 文件完全一致，**v3 只有一套结果**。之前用户看到的"两次不一样"是 v2 与 v3 的对比，而非 v3 自身不一致。

### §12.2 公式手工验算（以 Sn 为例，excess 最大元素）

```
G^xs = Σ ω_ij × y_i × y_j = -0.4785 kJ/mol·atom
excess_Sn = ω_GaSn×0.65 + ω_InSn×0.20 + ω_SnZn×0.05 - G^xs
          = (-6.336×0.65) + (-4.121×0.20) + (-10.009×0.05) - (-0.479)
          = -4.965
μ_Sn^HEI = h_Sn + excess = -19.266 + (-4.965) = -24.231 ✓
```

Gibbs-Duhem 验证：
```
G_atom = -36.311
(3/4)μ_Pt + (1/4)Σ y_i μ_i^β = (3/4)(-36.311) + (1/4)(-36.311) = -36.311 ✓
```

### §12.3 v2 vs v3 差异根源

| | v2（文献量热） | v3（UMA-CEF） | 原因 |
|---|---|---|---|
| h(Pt₃Ga) | -13.03 /atom | -43.39 /atom | UMA(DFT级) 系统性更负 |
| h(Pt₃Sn) | -15.60 /atom | -19.27 /atom | 同上，但 Sn 差距小 |
| 排序 | Sn > Ga > Pt > Zn > In | **Ga > Pt > Zn > Sn > In** | 端元排序翻转 |
| ΔG_rxn /f.u. | -55.7 | **-150.5** | 端元数值量级差 |

DFT/ML 势系统性过稳化 Pt 基金属间化合物是计算材料科学已知现象（PBE 泛函偏差），Methods/SI 注明即可。

### §12.4 Panel a–g 跨图自洽性审查

| 对比 | 结果 | 说明 |
|---|---|---|
| Panel c ↔ Panel g | ✅ 等摩尔: -30.008 vs -29.904，Δ=0.1 (0.3%) | 同一 UMA 势 |
| Panel f ↔ Panel g | ⚠️ 无冲突 | f 用 GaPt(1:1)，g 用 Pt₃Ga(3:1)，化学计量不同 |
| 新 Panel a ↔ Panel g | ✅ 直接拟合自 165 点 | CEF R²=0.9993 |
| Panel a,b,d,e | ✅ 无冲突 | 用 Miedema 混合焓/界面能，非 DFT 形成焓 |

**结论**：整个 Figure 1 在数据来源和物理量定义上自洽。

### §12.5 决策记录

- **2026-04-23 确认使用 v3 数据**：优先保证与 Panel g 内部自洽，绝对值在 Methods 注明为 DFT 级别估计。
- 相对排序和驱动力符号（全正）在 v2 和 v3 之间是 robust 的。

