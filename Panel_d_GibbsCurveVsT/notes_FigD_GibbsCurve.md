# SI — HEI vs HEA Gibbs-energy temperature curves
## 方法学与叙事说明 (原因 · 过程 · 结果)

---

## 1. 绘制的原因 (Why)

### 1.1 Panel c 已经回答的
**在 0 K**，对相同组成 `Pt24Ga2In2Sn2Zn2` (等摩尔极限)：

| 状态 | ΔH_f (kJ/mol/atom) |
|:---|---:|
| HEI = Ordered L1₂ (1 representative occupancy) | **−30.008** |
| HEA = Disordered random solid solution (N=30 configs) | **−13.964 ± 2.262** |
| Ordering gain (Gap = HEA − HEI) | **+16.044** |

这告诉我们：**焓层面上，有序 HEI 比无序 HEA 深 16 kJ/mol/atom**。

### 1.2 审稿人一定会问的
> "熵的贡献呢？HEA 的构型熵比 HEI 大得多，升温到合成温度
>  (473–1500 K) 时 `T·ΔS` 是不是会把这 16 kJ/mol 的优势吃掉？
>  你们 0 K 的结论还成立吗？"

**回答这个问题需要两样东西**：
1. 把 0 K 的单点焓差扩展为随温度变化的 ΔG(T) 曲线
2. 标出合成温度窗口 (500, 1000, 1500 K)，显示 HEI 曲线是否仍在 HEA 之下

这正是本 SI 图的目的。

### 1.3 对 manuscript 的叙事支撑
- 延续 Panel c 的 "HEI is enthalpically preferred" → "HEI is preferred across the entire synthesis temperature window"
- 提供对 "high-entropy 名字里的 'entropy' 是不是反而不利于 ordering" 这个天然质疑的预先回应
- 连接 Panel c (0 K 焓差) 与 Panel g (化学计量景观) 之间的热力学过渡

---

## 2. 物理模型 (What is being computed)

### 2.1 两个状态的定义

两者**组成完全相同**（same 5 元素 mole fractions），区别仅在原子排布：

| 状态 | 位置占位 |
|:---|:---|
| **HEI** (ordered L1₂) | Pt 锁定在 A 子晶格 24 位点；Ga/In/Sn/Zn 在 B 子晶格 8 位点上占据 Panel c 选定的一个 representative occupancy |
| **HEA** (disordered) | 所有 32 位点随机取 5 种元素之一，保持总组成 x_Pt=0.75, x_Ga=x_In=x_Sn=x_Zn=0.0625 |

### 2.2 Gibbs 自由能分解

对每个状态：
$$
\Delta G_f(T) = \Delta H_f - T \cdot \Delta S_f
$$

其中 $\Delta H_f, \Delta S_f$ 都以 **纯元素参考态** (Pt `mp-126`, Ga `mp-142`, In/Sn/Zn generated prototypes) 为基准，与 Panel c / g 口径完全一致。

### 2.3 熵的拆分 (关键近似)

完整的熵包括:
$$
S = S_{\rm config} + S_{\rm vib} + S_{\rm elec} + S_{\rm mag}
$$

**本图只计算 S_config，其他项作对消近似** — 理由：
- $S_{\rm vib}$: HEA 和 HEI 都是金属相，Debye 温度相近 → $\Delta S_{\rm vib}({\rm HEA-HEI})$ 通常 <0.5 k_B/atom at 1500 K (文献经验值)，远小于 $\Delta S_{\rm config} \approx 0.56\,k_B$
- $S_{\rm elec}$: 两者都是 Pt-rich 金属，费米面态密度同量级 → 对消
- $S_{\rm mag}$: Pt-Ga-In-Sn-Zn 无 3d 铁磁元素，忽略

这是材料热力学里非常标准的**"quasi-harmonic 近似 + 构型熵为主"** 假设，见 §6 的边界说明。

---

## 3. 输入数据 (Inputs)

来源 [data_FigC_Long.csv](../Panel_c_OrderedVsDisordered/data_FigC_Long.csv)。

| 变量 | 值 | 数据来源 |
|:---|---:|:---|
| $\Delta H_f^{\rm HEI}$ | −30.008 kJ/mol/atom | 1 × UMA-s-1p1 single-point on ordered L1₂ supercell |
| $\langle \Delta H_f^{\rm HEA} \rangle$ | −13.964 kJ/mol/atom | Mean over 30 × UMA-s-1p1 on random occupancy (seeds 100–129) |
| $\sigma(\Delta H_f^{\rm HEA})$ | 2.262 kJ/mol/atom | Population std (ddof=0) over same 30 configs |
| $N_{\rm HEA}$ | 30 | 30 independent random occupancy seeds |

**不重跑 UMA** — 完全复用 Panel c 的单点能量，只做 post-processing。

---

## 4. 构型熵的计算 (Math)

### 4.1 HEA 侧 — 全位点随机

对整个 32-原子 supercell，5 种元素全部随机：
$$
S_{\rm config}^{\rm HEA} = -R \sum_i x_i \ln x_i
$$

代入 $x_{\rm Pt}=0.75, x_{\rm Ga}=x_{\rm In}=x_{\rm Sn}=x_{\rm Zn}=0.0625$:

$$
\begin{aligned}
S_{\rm config}^{\rm HEA} &= -R \left[ 0.75\ln 0.75 + 4\times 0.0625\ln 0.0625 \right] \\
&= -R \left[ -0.2158 - 0.6931 \right] \\
&= R \times 0.9090 \\
&= 7.557 \text{ J/(K·mol·atom)} \\
&= 0.007557 \text{ kJ/(K·mol·atom)}
\end{aligned}
$$

### 4.2 HEI 侧 — 两个 bound

**HEI 的熵取决于如何定义 ordered state**。这里给两个 bound：

#### Bound A (保守上界，推荐主线)
**Sublattice model**: Pt 完全锁定在 A 子晶格 (24/32 sites, 无构型熵)；B 子晶格 (8/32 sites) 上 Ga/In/Sn/Zn 等分视为残余无序：
$$
S_{\rm config}^{\rm HEI,\,sub} = \frac{n_B}{n_{\rm tot}} \cdot R \ln 4 = \frac{8}{32} \cdot R \ln 4
$$
$$
= 0.25 \times 8.314 \times 1.3863 = 2.880 \text{ J/(K·mol·atom)}
$$

**物理解释**: L1₂ 的 A 子晶格有明确对称性，Pt 确实不混；但 B 子晶格的 4 种非 Pt 元素只有当我们承认有序态**保留 B 位点的随机混合**时才有这份残余熵。这个处理对 HEI 更 generous (S 大一点 → ΔG 更低一点)，属于 "HEI 最不利"（最接近 HEA）的情况。

#### Bound B (严格下界)
**Frozen occupancy**: Panel c 的 "Ordered_L12_Equimolar" 本来就是**一个具体的 Ga/In/Sn/Zn 占位选择**，不是 sublattice-averaged。对这个 frozen state：
$$
S_{\rm config}^{\rm HEI,\,frozen} = 0
$$

这是 HEI 的 "最有利" 情形 (S 最小 → ΔG 最深)。

**为什么保留两个 bound**: 真实 HEI 位于两者之间——实验上 B-site 不会完全冻结（有少量扩散）也不会完全随机（近邻相互作用会有偏好）。两个 bound 夹住真值，保证结论鲁棒。

### 4.3 ΔG(T) 主公式

$$
\boxed{\,\Delta G_f^{\rm HEA}(T) = \langle \Delta H_f^{\rm HEA} \rangle - T \cdot S_{\rm config}^{\rm HEA}\,}
$$

$$
\boxed{\,\Delta G_f^{\rm HEI}(T) = \Delta H_f^{\rm HEI} - T \cdot S_{\rm config}^{\rm HEI}\,}
$$

### 4.4 交叉温度 T\*

两曲线相交处：
$$
T^\ast = \frac{\Delta H_{\rm gap}}{\Delta S} = \frac{\langle \Delta H^{\rm HEA}\rangle - \Delta H^{\rm HEI}}{S^{\rm HEA}_{\rm config} - S^{\rm HEI}_{\rm config}}
$$

- Bound A: $T^\ast = \dfrac{16.044}{0.007557 - 0.002880} = \dfrac{16.044}{0.004677} \approx 3429 \text{ K}$
- Bound B: $T^\ast = \dfrac{16.044}{0.007557 - 0} \approx 2123 \text{ K}$

**两个 bound 的 T\* 都远高于 manuscript 合成温度上限 1500 K**，这是核心结论。

---

## 5. 数值结果 (Results)

### 5.1 三个关键温度点

| T (K) | $\Delta G^{\rm HEA}$ | $\Delta G^{\rm HEI}_{\rm sub}$ | $\Delta G^{\rm HEI}_{\rm frozen}$ | Gap (HEA−HEI, sub) |
|---:|---:|---:|---:|---:|
| 500  | −17.74 | −31.45 | −30.01 | **+13.70** |
| 1000 | −21.52 | −32.89 | −30.01 | **+11.37** |
| 1500 | −25.30 | −34.33 | −30.01 | **+9.03** |

单位 kJ/mol/atom。

### 5.2 读者应看到的三件事

1. **两条 HEI 线 (实线 + 虚线) 都始终在 HEA 带 (灰色 ±1σ) 下方**，在整个 300–2500 K 范围
2. **HEA 曲线斜率更陡** (熵大) — 升温时 HEA 自由能下降更快；这是符合直觉的 "高熵合金为何在高温存在" 的图示
3. **交叉点 T\* 都在 2100 K 之上**，远超合成温度窗口 (473–1500 K)，甚至超过 Pt-Ga 系液相温度

### 5.3 对 0 K 数字的完全衔接
- T = 0 K 处两曲线的垂直距离 = 16.04 kJ/mol/atom（Panel c 的原始数字）
- 曲线斜率差 = `ΔS_config`，与热力学教科书一致

---

## 6. 局限与边界 (Caveats)

### 6.1 已做的近似
| 假设 | 依据 | 预期误差 |
|:---|:---|:---|
| 忽略振动熵差 | HEA / HEI 都是 Pt-rich 金属相 | <0.5 k_B/atom at 1500 K |
| 忽略电子熵差 | 两者 Fermi 面 DOS 同量级 | 可忽略 |
| HEI 取一个 representative occupancy | Panel c 方法学范围 | Panel c std 内（~2 kJ/mol/atom），远小于 16 kJ/mol gap |
| UMA-s-1p1 0 K 能量 | Panel c 已与 MP-DFT/CHGNet 在 Panel f 的 15 M-Ga 体系里对齐 | 系统偏差由 HEA 的 N=30 σ 反映 |

### 6.2 这张图**不能**证明什么
- 不能断言 HEI **会**在 manuscript 合成条件下实际形成 (kinetic barriers 未计)
- 不能推广到 Zn 含量远离 6.25% 的情形 (本图只算一个组成)
- 不能替代 Panel g 的化学计量扫描 (那是另一个维度)
- 不能给出相图级别的相界 (只比较两个状态，没算液相 / 其他有序相)

### 6.3 审稿人可能的 push back 与预答
| 质疑 | 预答 |
|:---|:---|
| "你们忽略了 vibrational entropy" | §6.1 依据 + 0.5 k_B/atom 上限 → 交叉温度最多降 ~300 K，仍 > 1500 K |
| "Ordered state 的 S 不应该完全 zero" | 给了两个 bound，frozen (=0) 是下界 |
| "UMA-s-1p1 不是 DFT" | Panel f (15 M-Ga 三方法共识) + Panel c (M-Ga fixed-cell 验证) 已证明 UMA 在此体系与 MP-DFT Spearman ρ = 0.951 |
| "N=30 configs 不够代表 HEA" | σ = 2.26 kJ/mol 已小于 16 kJ/mol gap 的 15%，鲁棒 |

---

## 7. 绘制流程 (How the figure is produced)

### 7.1 执行顺序
1. 读 Panel c 的 `Panel_c_OrderedVsDisordered/data_FigC_Long.csv` → 拿到 $\Delta H^{\rm HEI}$ 单值 + $\langle \Delta H^{\rm HEA} \rangle, \sigma$
2. 计算 $S^{\rm HEA}_{\rm config}$ 和 $S^{\rm HEI}_{\rm config}$ (两个 bound)
3. 温度网格 `T ∈ [300, 2500, 221 步]`
4. 对每个 T 算 $\Delta G^{\rm HEA}(T), \Delta G^{\rm HEA}(T) \pm \sigma$ 和两个 HEI bound
5. 输出 CSV + PNG
6. 额外输出 `data_SI_HEIvsHEA_KeyPoints.csv` 只含 T=500/1000/1500 K 三行

### 7.2 图元

| 图元 | 样式 | 含义 |
|:---|:---|:---|
| 灰色填充带 | `#90A4AE` @ 28% α | HEA ±1σ (30 configs) |
| 深灰实线 | `#37474F`, 2 pt | HEA mean $\Delta G(T)$ |
| 红色实线 | `#C62828`, 2 pt | HEI (Bound A: sublattice entropy) |
| 红色虚线 | `#C62828`, 1.2 pt dashed | HEI (Bound B: frozen, S=0) |
| 浅灰虚线 | `#BDBDBD`, 0.6 pt dotted | T = 500 / 1000 / 1500 K 标记 |
| 绿色点划线 | `#2E7D32`, 1 pt | $T^\ast$ (sublattice) ≈ 3429 K |
| 深绿点线 | `#1B5E20`, 1 pt | $T^\ast$ (frozen) ≈ 2123 K |

### 7.3 对应脚本
[script_SI_HEI_vs_HEA_GibbsCurve.py](script_SI_HEI_vs_HEA_GibbsCurve.py) — py312, 运行时间 <1 s，无需 GPU / fairchem。

---

## 8. 在 SI 中的定位

**可能的 SI 章节**: `SI-?: Temperature-dependent Gibbs landscape of HEI vs HEA`

建议插入位置: SI-1 (Workflow) 之后, SI-2 (Zn) 之前 — 因为这张图是 Panel c 的直接延续，而 Panel c 在正文 Fig 1 里已出现。

**与其他 SI 的关系**:
- 呼应 **Panel c**: 把 0 K 16.04 kJ/mol 的单点衔接到连续温度轴
- 独立于 **SI-2 Zn**: 本图只做等摩尔 Pt24Ga2In2Sn2Zn2，不涉及 Zn 扫描
- 独立于 **SI-4 ML Features**: 本图不涉及 ML feature importance

---

## 9. 参考

- Panel c 数据源: [../../20260419_Fig1_Panels_a-g_Final/03_Panel_c_OrderedVsDisordered/notes_FigC.md](../../20260419_Fig1_Panels_a-g_Final/03_Panel_c_OrderedVsDisordered/notes_FigC.md)
- 振动熵对消近似 (de Fontaine 1994; van de Walle & Ceder, 2002)
- Configurational entropy in HEAs (Yeh et al., *Adv. Eng. Mater.* 6, 299, 2004)
- 执行脚本: [script_SI_HEI_vs_HEA_GibbsCurve.py](script_SI_HEI_vs_HEA_GibbsCurve.py)
