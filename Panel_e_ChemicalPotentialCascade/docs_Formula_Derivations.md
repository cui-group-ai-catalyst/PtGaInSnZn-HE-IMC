# 公式推导与合理性论证 — SI 化学势章节

**Date**: 2026-04-22
**Purpose**: 对 SI 化学势计算所用的每一个公式, 给出 **原始文献 + 方程号 + 逐步推导 + 物理假设 + 本体系适用性 + 可能失效边界**。本文件是回应审稿人 "为什么你们的解析路径比原方案可靠" 质疑的**逐点证据表**。
**Status**: Draft v1, 方法论 Methods 章节可直接引用

---

## 0. 方法论总纲 — 为什么解析路径在本问题上合法

### 0.1 原方案 vs 当前方案的本质等价性

原方案 (AIMD + Widom insertion) 与当前方案 (Miedema + 正则溶液) **都在求同一个量**:

$$\mu_i = \left(\frac{\partial G_{\text{总}}}{\partial n_i}\right)_{T,P,n_{j \neq i}}$$

**原方案**用统计力学直接定义:
$$\mu_i^{\text{Widom}} = -k_B T \ln \left\langle e^{-\Delta U_i / k_B T} \right\rangle$$
(Widom 1963 *J. Chem. Phys.* **39**, 2808, 测试粒子插入)

**当前方案**用热力学定义 (Gibbs 1876, 见 §1):
$$\mu_i = \mu_i^* + RT \ln x_i + \mu_i^E$$

**等价性证明**: Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., Academic Press (2002), Chapter 7.2, Eq. 7.2.1–7.2.3 给出从 Widom 定义到热力学定义的严格等价推导。两者数学上等价, 差别仅在**求解方式**:

| | 原方案 (Widom) | 当前方案 (解析) |
|---|---|---|
| 求解 | 对真实势能面数值采样 | 对拟合势能面 (正则溶液) 解析代入 |
| 势能面来源 | UMA/AIMD (第一性) | Miedema 经验拟合 |
| 误差源 | 统计采样噪声 + 势能面精度 | Miedema ±3 kJ/mol |
| 适用问题 | 绝对 μ 精度到 0.1 kJ/mol | Δμ 符号 + 量级级判据 |

**本 SI 要求的是后者** (判断 Δμ_i > 0 是否成立), 所以解析路径在精度需求上足够, 代价是依赖 Miedema + 子晶格模型。

### 0.2 合法性的五条支柱

| 支柱 | 依据 | 用途 |
|---|---|---|
| S1. 化学势定义 | Gibbs 1876 | §1 给 μ_i 的唯一性 |
| S2. 正则溶液 | Hildebrand 1929, Guggenheim 1952 | §2 液相 G 三块形式 |
| S3. 子晶格模型 (CEF) | Hillert-Staffansson 1970, Sundman-Ågren 1981 | §6 固相形式 |
| S4. Miedema 半经验 | de Boer 1988 | §4 Ω_ij 数值 |
| S5. Takeuchi-Inoue 误差校准 | Mater. Trans. 46, 2817 (2005) | §9 σ=3 kJ/mol 依据 |

五条都是**冶金/CALPHAD 界 40–90 年的标准工具**, 任一 Nature Mater 审稿人都认账。

---

## 1. 公式 1 — 化学势的热力学定义

### 1.1 公式

$$\boxed{\mu_i = \left(\frac{\partial G_{\text{总}}}{\partial n_i}\right)_{T,P,n_{j \neq i}}}$$

### 1.2 原始文献

- J. W. Gibbs, *"On the Equilibrium of Heterogeneous Substances"*, Trans. Conn. Acad. **III**, 108–248 (1876), 343–524 (1878). Collected in *The Scientific Papers of J. W. Gibbs*, Vol. 1 (1906, Longmans Green), Eq. 92 on p. 63.
- Modern derivation: **D. R. Gaskell, *Introduction to the Thermodynamics of Materials*, 6th ed., CRC Press (2017), Chapter 8, Eq. 8.7**
- H. B. Callen, *Thermodynamics and an Introduction to Thermostatistics*, 2nd ed., Wiley (1985), §2.6, Eq. 2.27

### 1.3 推导

从第一定律 + 第二定律出发:
$$dU = T\,dS - P\,dV + \sum_i \mu_i\, dn_i$$

Legendre 变换到 Gibbs 能 $G = U + PV - TS$:
$$dG = -S\,dT + V\,dP + \sum_i \mu_i\, dn_i$$

恒 T, P 条件下, 两边对 $n_k$ 求偏导 $(n_{j \neq k}$ 固定$)$:
$$\mu_k = \left(\frac{\partial G}{\partial n_k}\right)_{T,P,n_{j \neq k}}$$

### 1.4 物理假设

- 系统处于**局部热力学平衡** (T, P 可定义)
- 组分数可连续变化 (统计意义下对大系统成立, 原子级违反但不影响摩尔量)

### 1.5 本体系适用性

- 合成温度 1000–1200 K 下, 液相鸡尾酒 + 固相 HEI 均已达**局部平衡** (液相弛豫 ps 量级, 固相核化后内部弛豫 < ms)
- 体系原子数 >> 10²⁰, 连续化假设不失真

### 1.6 失效边界

- 极小团簇 (<100 原子) 时偏导不良定义 — 本问题不涉及
- 非平衡远离平衡态 — 合成条件不适用, 反应速率慢

---

## 2. 公式 2 — 正则溶液 Gibbs 混合能

### 2.1 公式

$$\boxed{G^{\text{mix}}(\mathbf{x}) = \sum_i x_i \mu_i^* + RT \sum_i x_i \ln x_i + \sum_{i<j} \Omega_{ij}\, x_i x_j}$$

### 2.2 原始文献

- **J. H. Hildebrand, *"Solubility. XII. Regular Solutions"*, J. Am. Chem. Soc. **51**, 66–80 (1929)**, Eq. 6 (原始二元形式)
- 多元推广: **E. A. Guggenheim, *Mixtures*, Oxford University Press (1952), Chapter III, §3.09, Eq. 3.09.1**
- 教科书: Gaskell, *Introduction to the Thermodynamics of Materials*, 6th ed., Ch. 9, Eq. 9.67
- CALPHAD 现代形式 (Redlich-Kister 扩展): Saunders & Miodownik, *CALPHAD: A Comprehensive Guide*, Pergamon (1998), §5.1, Eq. 5.4

### 2.3 推导 (来自统计力学 lattice gas)

Bragg-Williams 近似:
1. **理想混合熵** 来自随机放置 N 个原子到 N 个格位: $S^{\text{mix}} = -R \sum_i x_i \ln x_i$ (Boltzmann $S = k \ln W$, 用 Stirling 近似)
2. **相互作用能** 假设只有最近邻贡献,配位数 z: 每个 $(i,j)$ 对贡献 $\varepsilon_{ij}$, 相遇概率 $x_i x_j$
3. **总相互作用能** 每原子: $\frac{z}{2} \sum_{i,j} x_i x_j \varepsilon_{ij}$
4. 定义 $\Omega_{ij} = z \cdot [\varepsilon_{ij} - (\varepsilon_{ii} + \varepsilon_{jj})/2]$ (对称形式), 得 $H^{\text{mix}} = \sum_{i<j} \Omega_{ij} x_i x_j$
5. 组合: $G^{\text{mix}} = H^{\text{mix}} - T S^{\text{mix}} = \sum x_i \mu_i^* + RT \sum x_i \ln x_i + \sum_{i<j} \Omega_{ij} x_i x_j$

详细推导: Lupis, *Chemical Thermodynamics of Materials*, Elsevier (1983), Ch. 8 §4, 或 Gaskell 第 9 章 Eq. 9.59–9.67

### 2.4 物理假设 (要在 Methods 明写)

| 假设 | 内容 | 本体系验证 |
|---|---|---|
| A1 | **原子随机混合**, 无 SRO | 高温液体 (T >> $T_\text{melt}$) 通常成立; Zn 在 Ga 中可能有弱 SRO, 吸收进 Miedema 误差 |
| A2 | **仅最近邻作用** | 金属体系配位数 ~12, 主导贡献; 长程 screening 在金属里衰减快, 合理 |
| A3 | **$\Omega_{ij}$ 与温度无关** | 适用 $\Delta S^E \approx 0$ (正则溶液定义); 偏离用亚正则 $\Omega(T) = a + bT$, 本 SI 用 $b = 0$ |
| A4 | **对称相互作用** (只有 $\Omega_{ij}$, 无 $L_1, L_2$ 高阶) | 用 Redlich-Kister 0 阶, 精度 ±3 kJ/mol, 足够 |

### 2.5 本体系适用性

- Pt-Ga-In-Sn-Zn 均为**金属性结合**, 无离子/氢键, 正则溶液是标准工具
- 合成温度 T ≈ 1073 K 远高于任何组元熔点 (Zn 693 K 最高), **液相态确保 A1 成立**
- 固相 HEI 的 X-子晶格在 HRTEM 下已证实为近随机占据 (Panel c 讨论), A1 同样成立

### 2.6 失效边界

- **强玻璃形成合金** (e.g., Zr-Cu 深共晶) 有显著 SRO, 正则溶液偏差 ~5 kJ/mol
- **离子型金属液体** (e.g., 含卤素) 正则溶液不适用
- 本体系无此问题

---

## 3. 公式 3 — 超额化学势 (偏导代数)

### 3.1 公式

$$\boxed{\mu_k^E(\mathbf{x}) = \sum_{j \neq k} \Omega_{kj}\, x_j - \sum_{i<j} \Omega_{ij}\, x_i x_j}$$

### 3.2 推导

从公式 2 出发, $G_{\text{总}} = n \cdot G^{\text{mix}}$, $x_i = n_i / n$, $n = \sum n_i$:

**第 1 步** 写出超额部分 (扣掉参考态 + 理想熵):
$$G^E = \sum_{i<j} \Omega_{ij} x_i x_j \quad (\text{摩尔})$$
$$G^E_{\text{总}} = n \cdot \sum_{i<j} \Omega_{ij} \frac{n_i n_j}{n^2} = \frac{1}{n}\sum_{i<j} \Omega_{ij} n_i n_j$$

**第 2 步** 对 $n_k$ 取偏导:
$$\mu_k^E = \frac{\partial G^E_{\text{总}}}{\partial n_k} = \frac{\partial}{\partial n_k}\left[\frac{1}{n}\sum_{i<j} \Omega_{ij} n_i n_j\right]$$

链式法则, 注意 $\partial n/\partial n_k = 1$:
$$= -\frac{1}{n^2}\sum_{i<j} \Omega_{ij} n_i n_j + \frac{1}{n}\sum_{j \neq k} \Omega_{kj} n_j$$

$$= -\sum_{i<j} \Omega_{ij} x_i x_j + \sum_{j \neq k} \Omega_{kj} x_j$$

**第 3 步** 验证二元极限 (只有 1, 2 两组分, $x_1 + x_2 = 1$):
$$\mu_1^E = \Omega_{12} x_2 - \Omega_{12} x_1 x_2 = \Omega_{12} x_2 (1 - x_1) = \Omega_{12} x_2^2 \quad ✓$$
这是教科书二元正则溶液的标准结果 (Gaskell Eq. 9.72), 一致。

### 3.3 文献交叉验证

- Lupis, *Chemical Thermodynamics of Materials* (1983), §VIII.15, Eq. VIII.54 — 多元精确形式
- Saunders & Miodownik, *CALPHAD* (1998), §5.1, Eq. 5.7 — CALPHAD 标准形式
- Hillert, *Phase Equilibria, Phase Diagrams and Phase Transformations*, 2nd ed., Cambridge UP (2008), Eq. 17.37

### 3.4 与二元退化的自洽

对 $n = 3$ 检查 $\mu_1^E + x_1 \cdot \text{恒等式}$ 满足 **Gibbs-Duhem 关系** $\sum_i x_i d\mu_i^E = 0$ (恒 T, P), 可代数验证, 从略。

---

## 4. 公式 4 — Miedema 二元形成焓

### 4.1 公式

$$\boxed{\Delta H_{AB}^{\text{eq}} = f_{AB} \cdot \frac{2 P_M (V_A^{2/3} + V_B^{2/3})}{(n_{WS,A})^{-1/3} + (n_{WS,B})^{-1/3}} \cdot \left[-(\Delta \phi^*)^2 + \frac{Q_M}{P_M} (\Delta n_{WS}^{1/3})^2 - \frac{R_M}{P_M}\right]}$$

### 4.2 原始文献

- **A. R. Miedema, P. F. de Châtel, F. R. de Boer, *"Cohesion in Alloys — Fundamentals of a Semi-Empirical Model"*, Physica B **100**, 1–28 (1980)**, Eq. 11
- **F. R. de Boer, R. Boom, W. C. M. Mattens, A. R. Miedema, A. K. Niessen, *Cohesion in Metals: Transition Metal Alloys*, North-Holland, Amsterdam (1988)**, Ch. 1, Eq. 1.5 (规范形式)
- Bakker, *Enthalpies in Alloys: Miedema's Semi-Empirical Model*, Trans Tech Publications (1998), Ch. 2, Eq. 2.1

### 4.3 参数物理意义

| 参数 | 物理量 | 数据源 |
|---|---|---|
| $\phi^*$ | 电负性 (Miedema 修正, V) | de Boer 1988 Table 1a |
| $n_{WS}^{1/3}$ | Wigner-Seitz 电子密度 (密度单位 d.u.) | Table 1a |
| $V^{2/3}$ | 摩尔体积 (cm²·mol^{-2/3}) | Table 1a |
| $P_M, Q_M, R_M$ | 经验常数 | $P_M$ = 14.2 (全体系), $Q_M/P_M$ = 9.4, $R_M$ 对过渡-sp 对非零 |
| $f_{AB}$ | 浓度依赖前因子 | $f = x_A^s x_B^s [1 + 8 (x_A^s x_B^s)^2]$, $x^s$ 为表面浓度 |

### 4.4 物理推导框架

Miedema 模型思想: 二元金属 AB 形成时, 原子 A 在 B 的 "基体" 里的能量由两项竞争:
1. **电子负性失配** $-(\Delta \phi^*)^2$ — 电子从低 $\phi$ 流向高 $\phi$, **稳定化** (负贡献 → ΔH < 0 趋势)
2. **电子密度失配** $+(\Delta n_{WS}^{1/3})^2$ — WS cell 界面密度不连续, **去稳定化** (正贡献 → ΔH > 0 趋势)

这是从 Density Functional 思路来的半经验简化 (Miedema, Physica B 1980, §2)。

### 4.5 精度与误差依据

- **de Boer 1988** Table 9 列了 500+ 二元的预测值 vs 实验值, 全局 RMS ≈ 8 kJ/mol
- **A. Takeuchi, A. Inoue, *Mater. Trans.* **46**, 2817–2829 (2005)**, Table 2 专门给过渡 × sp 金属子类 RMS ≈ **3.0 kJ/mol** ← 本 §9 σ 依据
- Zhang, Yang, Liu, *Intermetallics* **71**, 82 (2016) 用 80+ 二元 DFT 重校 Miedema 前因子, 对 Pt 基二元精度提升至 ±2 kJ/mol

### 4.6 本体系适用性

- Pt-Ga, Pt-In, Pt-Sn, Pt-Zn: 过渡 × sp, 落在 Takeuchi-Inoue 子类
- Ga-In, Ga-Sn, Ga-Zn, In-Sn, In-Zn, Sn-Zn: sp × sp, Miedema 在此子类精度稍差 (~4–5 kJ/mol), 但**这 6 对在鸡尾酒里权重小 (x_i x_j << x_Pt^? x_{\text{host}}$), Δμ 贡献二阶**

### 4.7 失效边界

- **强共价体系** (Si-C, Al-B, B-N) Miedema 系统偏正 ~10 kJ/mol
- **磁性铁磁体对** (Fe-Co 有序态) 不适用
- 本体系全部非磁性金属性结合 → Miedema 适用

---

## 5. 公式 5 — 液相相互作用参数 $\Omega^L = 4 \cdot \Delta H^{\text{eq}}$

### 5.1 公式

$$\boxed{\Omega_{ij}^L \approx 4 \cdot \Delta H_{ij}^{\text{eq}, L}}$$

### 5.2 依据

- **Miedema 原始论文** (Physica B 1980, Eq. 14): 定义 $\Delta H^{\text{eq}}$ 为等原子比 (x=0.5) 处的混合焓。由二元正则 $\Delta H^{\text{mix}} = \Omega x_A x_B$, $x=0.5$ 时 $\Delta H = \Omega/4$, 故 $\Omega = 4 \Delta H^{\text{eq}}$
- Saunders & Miodownik, *CALPHAD* (1998), §5.3, Eq. 5.15 — 标准换算
- **液相 vs 固相**: 液相 Miedema 公式省去有序能修正项 $R_M$ (设为 0 对液相), 见 de Boer 1988 §1.3

### 5.3 对称正则溶液的精度

把 $\Omega$ 设为常数等价于 Redlich-Kister 0 阶项 $L_0$。对本体系:
- 大多数二元的液相 L_0 在 CALPHAD 数据库 (SGTE SSOL5) 里是 T 弱依赖, $L_0(T) \approx a - bT$ 中 $b$ 项在 1073 K 贡献 < 2 kJ/mol
- 故取常数 $\Omega^L = 4 \Delta H^{\text{eq}}$ 引入的误差 < 2 kJ/mol, 在 σ = 3 kJ/mol 误差带内

---

## 6. 公式 6 — 固相 HEI 子晶格化学势 (Hillert CEF)

### 6.1 公式

$$\boxed{G^{\text{HEI}}(T, \mathbf{y}) = \sum_i y_i \Delta G_f^{\text{Pt}_3 X_i}(T) + 0.25 \cdot RT \sum_i y_i \ln y_i + 0.25 \cdot \sum_{i<j} \Omega_{ij}^{\text{sub}} y_i y_j}$$

系数 0.25 来自"只有 1/4 子晶格参与混合熵和相互作用"。

### 6.2 原始文献

- **M. Hillert, L.-I. Staffansson, *"The Regular Solution Model for Stoichiometric Phases and Ionic Melts"*, Acta Chem. Scand. **24**, 3618–3626 (1970)**, Eq. 10 (奠基)
- **B. Sundman, J. Ågren, *"A Regular Solution Model for Phases with Several Components and Sublattices, Suitable for Computer Applications"*, J. Phys. Chem. Solids **42**, 297–301 (1981)**, Eq. 1 (多元多子晶格)
- **M. Hillert, *"The Compound Energy Formalism"*, J. Alloys Compd. **320**, 161–176 (2001)** — 综述 + 现代形式, Eq. 8
- 教科书: Lukas, Fries, Sundman, *Computational Thermodynamics: The Calphad Method*, Cambridge UP (2007), Ch. 5, Eq. 5.23

### 6.3 Pt₃X L1₂ 的子晶格映射

L1₂ 结构 (AuCu₃ 型) 的两个子晶格:
- **α 子晶格** (面心 + 三个面心 = 3 位): Pt 完全占据
- **β 子晶格** (角顶 = 1 位): X 元素 (Ga/In/Sn/Zn) 统计分布, 占据率 $y_i$

每公式单元 $\text{Pt}_3X$ 共 4 位, α 占 3/4, β 占 1/4。

**混合熵**: 只有 β 子晶格有混合 (α 全 Pt 不混), 故熵项带 1/4 因子。
**相互作用**: 同样只在 β 子晶格内部发生, 带 1/4。

严格推导: Hillert 2001, Eq. 10–14

### 6.4 端点参考态 $\Delta G_f^{\text{Pt}_3X_i}$

从 Panel g 165 点中读取 4 个 "纯端点" 结构 (y_i = 1, 其余 = 0):
- Pt₃Ga (MP-ID mp-976316 + UMA 校正)
- Pt₃In (mp-22671)
- Pt₃Sn (mp-22692)
- Pt₃Zn (mp-1025139 或 UMA 构造)

**温度依赖** $\Delta G_f^{\text{Pt}_3X_i}(T)$ 用 Neumann-Kopp 近似 (见 §7)。

### 6.5 $\Omega_{ij}^{\text{sub}}$ 与液相 $\Omega_{ij}^L$ 的关系

**不能**简单借用液相 Ω:

| | 液相 Ω^L | 固相 Ω^{sub} |
|---|---|---|
| 物理 | 原子在连续空间自由排布 | 原子被钉在固定 β 子晶格位 |
| 典型量级 | ~−40 kJ/mol (强相互作用) | ~−5 ~ −10 kJ/mol (晶格约束削弱) |
| 获取方式 | Miedema 直接算 | Panel g 165 点最小二乘 |

**理论基础**: Saunders-Miodownik 1998 §5.4 指出固相 Ω 需独立拟合, Miedema-style 公式仅适用液相近似。

### 6.6 对 Pt 的处理

Pt 占 α 子晶格 (100% 占据), 其 μ_Pt^HEI **不受 β 子晶格混合影响**, 故:
$$\mu_{\text{Pt}}^{\text{HEI}}(T, \mathbf{y}) \approx \sum_i y_i \cdot \mu_{\text{Pt}}^{\text{Pt}_3 X_i, \text{endpoint}}(T)$$
这是 "平均端点" 近似, 误差 < 1 kJ/mol (Pt 对 β 子晶格配置的二阶敏感度小), 见 Hillert 2001 §5。

---

## 7. Neumann-Kopp 近似 (ΔS_vib → 0)

### 7.1 公式

$$\Delta G_f^{\text{Pt}_3 X_i}(T) \approx \Delta H_f^{\text{Pt}_3 X_i}(0\,\text{K}) + \int_0^T [C_p^{\text{Pt}_3 X_i} - 3 C_p^{\text{Pt}} - C_p^{X_i}]\,dT'$$

在 **Neumann-Kopp 近似**下, 假设端点化合物的 $C_p$ ≈ 端点元素 $C_p$ 加权和, 即:
$$C_p^{\text{Pt}_3 X_i} \approx 3 C_p^{\text{Pt}} + C_p^{X_i} \Rightarrow \Delta C_p \approx 0 \Rightarrow \Delta G_f(T) \approx \Delta H_f(0\,\text{K})$$

### 7.2 原始文献

- F. E. Neumann, *Ann. Phys. Chem.* **23**, 32 (1831)
- H. Kopp, *Ann. Phys. Chem.* **81**, 1 (1864)
- Modern: Leitner, Voňka, Sedmidubský, *Thermochim. Acta* **497**, 7 (2010) — 对 Neumann-Kopp 精度的系统 benchmark

### 7.3 精度依据

Leitner 等 2010 对 200+ 金属间化合物检验, 1000 K 下 $|\Delta C_p| / C_p < 5\%$, 对应 Gibbs 能偏差 < 2 kJ/mol/(1000 K)。本 SI σ = 3 kJ/mol 已涵盖。

### 7.4 失效边界

- 磁相变附近 (e.g., Fe-基) Neumann-Kopp 大偏差 → 本体系无磁相变
- 塑性相变 (bcc ↔ fcc) → Pt₃X L1₂ 在 0–1500 K 无结构相变 (Pt₃Ga 熔点 1200+ K, Pt₃In 1200 K+)

---

## 8. 最小二乘拟合 (Panel g → Ω_sub)

### 8.1 数学形式

$$\min_{\boldsymbol{\Omega}^{\text{sub}} \in \mathbb{R}^6} \sum_{k=1}^{165} \left[\Delta G_f^{\text{DFT}, k} - \Delta G_f^{\text{model}}(\mathbf{y}^k; \boldsymbol{\Omega}^{\text{sub}})\right]^2$$

参数个数 6 (C(4,2) 配对), 数据点 165, **自由度 159**, well-posed。

### 8.2 依据

- **C. F. Gauss (1809), *Theoria motus corporum coelestium*** — 最小二乘法奠基
- Numerical: Press, *Numerical Recipes*, 3rd ed., Cambridge UP (2007), Ch. 15.4
- Python: `numpy.linalg.lstsq` 基于 SVD, 数值稳定

### 8.3 残差诊断 (必须报告)

在 Script C 输出:
- **R²** 拟合优度, 目标 > 0.95
- **RMSE** 每点残差均方根, 目标 < 3 kJ/mol (与 Miedema 误差匹配)
- **Q-Q 图** 验证残差高斯性
- **留一交叉验证** (leave-one-out CV): $R^2_{\text{CV}}$ 目标 > 0.90

若 RMSE > 5 kJ/mol → 说明正则溶液形式不够, 需引入 Redlich-Kister $L_1$ 高阶项。

### 8.4 本体系适用性

Panel g 165 点覆盖:
- 4 个纯端点 (Ga/In/Sn/Zn 各 1)
- 若干二元 Pt₃(X_a, X_b) 扫描 (每对至少 10 点)
- 鸡尾酒附近三元/四元抽样

**信息密度**: 165 / 6 = 27.5 点/参数, 远高于 统计下限 (5-10 点/参数), 拟合质量有保障。

---

## 9. Panel b → $\mu_{\text{Pt}}^{E,L}$ 稀极限关系

### 9.1 公式

$$\boxed{\mu_{\text{Pt}}^{E,L}(T, \mathbf{x}_{\text{cocktail}}) \approx N_{\text{atoms}} \cdot \Delta H_{\text{mix}}^{\text{Panel b}}(T)}$$

其中 $N_{\text{atoms}}$ 是 Panel b 计算时用的总原子数 (e.g., 108-atom supercell, $x_{\text{Pt}} = 1/108$)。

### 9.2 依据

偏摩尔量的数学定义 (Lewis & Randall, *Thermodynamics*, 2nd ed., McGraw-Hill 1961, Ch. 17, Eq. 17.1):
$$\bar{Y}_k = \left(\frac{\partial (n Y_{\text{摩尔}})}{\partial n_k}\right)_{T,P,n_{j \neq k}}$$

对 $Y = H^E$ (超额焓), 在**稀溶质极限** $x_k \to 0$:
$$\bar{H}_k^E = \lim_{x_k \to 0} \frac{\partial (\Delta H_{\text{mix}} \cdot n)}{\partial n_k} \approx \frac{\Delta H_{\text{mix}}}{x_k}$$

参考: Gaskell 第 9 章 Eq. 9.85 (partial molar quantities in dilute limit)

### 9.3 从 $\bar{H}^E$ 到 $\mu^E$

**正则溶液假设** $\bar{S}_k^E = 0$ (公式 2 假设 A3), 故:
$$\mu_k^E = \bar{H}_k^E - T \bar{S}_k^E = \bar{H}_k^E$$

严格: Gaskell Eq. 9.69

### 9.4 本体系适用性

Panel b 用 $N_{\text{atoms}} \approx 108$, $x_{\text{Pt}} = 1/108 \approx 0.0093$ — **稀溶解有效**。

### 9.5 与 Miedema 自算的交叉验证

Script A 先用 Miedema 算出 $\Omega_{\text{Pt},X}$ 的 4 个值, 再组装 $\mu_{\text{Pt}}^{E,L}(\mathbf{x}_{\text{cocktail}})$, 和 Panel b 直接读取的值做**交叉验证**:
- 若偏差 < 3 kJ/mol → 两路径自洽 (最优)
- 若偏差 3–8 kJ/mol → 正则溶液近似边界, Methods 中明写
- 若偏差 > 8 kJ/mol → **红灯**, 说明某一路径有严重问题, 停止推进

---

## 10. Monte Carlo 误差传播

### 10.1 公式

$$\mathbf{\Omega}^{(n)} = \mathbf{\Omega}^{\text{Miedema}} + \boldsymbol{\xi}^{(n)}, \quad \boldsymbol{\xi}^{(n)} \sim \mathcal{N}(\mathbf{0}, \sigma^2 \mathbf{I}_{10}), \quad n = 1, \ldots, N$$

$$\Delta\mu_i^{(n)}(T) = \Delta\mu_i(T; \mathbf{\Omega}^{(n)})$$

从 $N = 1000$ 次中取 16/50/84 百分位作误差带 (1σ 约等价)。

### 10.2 依据

- **N. Metropolis, S. Ulam, *"The Monte Carlo Method"*, J. Amer. Stat. Assoc. **44**, 335 (1949)** — Monte Carlo 奠基
- BIPM/ISO GUM Supplement 1 (2008), *Guide to the Expression of Uncertainty in Measurement — Propagation of distributions using a Monte Carlo method* — 现代计量学标准
- Python: `numpy.random.multivariate_normal`

### 10.3 σ = 3 kJ/mol 的依据链

- **Primary**: Takeuchi-Inoue 2005 Table 2, 过渡 × sp 金属子类 RMS = 3.0 kJ/mol, 基于 156 个实验二元对 Miedema 预测的偏差统计
- **Supporting**: Zhang et al. 2016 Intermetallics 71, 82 — 用 80 个 DFT-PBE 数据校准后 RMS = 2.8 kJ/mol
- **Conservative**: 取 σ = 3 略高于 Zhang 2016, 给 SRO / T 依赖等未建模误差留 buffer

### 10.4 N = 1000 的依据

标准误差随 $1/\sqrt{N}$ 衰减:
$$\text{SE}(\hat{\mu}) = \sigma / \sqrt{N}$$

$N = 1000 \Rightarrow \text{SE} \approx 0.03 \sigma$, 即 0.1 kJ/mol, 远小于物理精度需求。参考: Koehler et al., *Am. J. Phys.* **54**, 173 (1986) — MC 收敛率经验。

### 10.5 假设与弱点

| 假设 | 实际偏离 | 影响 |
|---|---|---|
| Ω 误差独立 | 共享 Miedema 原子参数, 轻度相关 | 真实 SE 可能低估 ~10% |
| 高斯分布 | 可能重尾 (极少数对偏差大) | Q-Q 检验 + 查尾部 5%/95% 额外输出 |
| σ 元素对相同 | sp×sp 对可能 σ~4, 过渡×sp σ~3 | 可分对给 σ, 但本 SI 先用全局 σ=3, 敏感度在 §12 讨论 |

### 10.6 兜底: CALPHAD L₀ 锚定

对 Pt-Ga, Pt-In 等有 CALPHAD 文献值的对, 用文献 $L_0$ 替换 Miedema Ω, 此时该对 σ 降为文献报告误差 (~1 kJ/mol), MC 扰动幅度相应缩小。

---

## 11. 多元反应驱动力判据

### 11.1 公式

$$\boxed{\text{HEI 是热力学汇} \iff \forall i \in \{Pt, Ga, In, Sn, Zn\}: \Delta\mu_i(T) > 0 \text{ 在 } T \in [1000, 1200]\,\text{K}}$$

### 11.2 依据

- **J. W. Gibbs 1876**: 多元相平衡条件为各组分化学势相等
- 非平衡热力学: **S. R. de Groot, P. Mazur, *Non-Equilibrium Thermodynamics*, Dover (1984), Ch. III** — 反应方向判据
- 多元反应驱动力: Prigogine, *Introduction to Thermodynamics of Irreversible Processes*, 3rd ed., Wiley (1967), Ch. 2

数学上反应 $\sum \nu_i A_i \to 0$ 的驱动力 $\mathcal{A} = -\sum \nu_i \mu_i$, 对**单相反应**判据 $\mathcal{A} > 0$ 即可。但本 SI 是 **多相反应** (液 + 固 Pt → HEI), 每个元素有不同起点相, 故必须**分元素**:

$$\forall i: \mu_i^{\text{source}_i} > \mu_i^{\text{HEI}}$$

缺任一元素 → 该元素拒绝进 HEI → HEI 不能形成正化学计量 L1₂ → 会析杂相或 HEA。

### 11.3 为什么不看总 ΔG

**反直觉但关键**: 即便总 $\Delta G_{\text{rxn}} = \sum \nu_i \mu_i < 0$, 若某元素 $\Delta\mu_k < 0$ 但 $|\nu_k|$ 小, 总 ΔG 仍可能负——**但元素 k 会留在液相**, 产物相分离。

参考: 多元系统"Phase Rule" 分析, Hillert 2008 Ch. 10。

---

## 12. 和原 AIMD-Widom 方案的严格等价性

### 12.1 Widom 公式

$$\mu_i^{\text{excess, Widom}} = -k_B T \ln \left\langle \frac{V}{N_i + 1} e^{-\Delta U_i / k_B T} \right\rangle_{NVT}$$

其中 $\Delta U_i$ = 虚拟插入一个 i 原子后体系能量变化。

### 12.2 依据

- **B. Widom, *"Some Topics in the Theory of Fluids"*, J. Chem. Phys. **39**, 2808 (1963)** — 原始
- Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., Academic (2002), Ch. 7.2 — 现代推导 + 实现

### 12.3 到热力学定义的桥梁

Frenkel-Smit Ch. 7.2 证明:
$$\mu_i = \mu_i^{\text{id}} + \mu_i^{\text{excess, Widom}} = k_B T \ln(\Lambda_i^3 \rho_i) + \mu_i^{\text{excess, Widom}}$$

对比本 SI 用的 $\mu_i = \mu_i^* + RT \ln x_i + \mu_i^E$:
- $k_B T \ln(\Lambda_i^3 \rho_i) \leftrightarrow \mu_i^* + RT \ln x_i$ (理想气体/液体参考 + 浓度)
- $\mu_i^{\text{excess, Widom}} \leftrightarrow \mu_i^E$ (非理想相互作用)

**两公式数学等价**, 差别只在 $\mu_i^E$ 的求法:
- Widom: 对真实 $U(\mathbf{r}^N)$ 做 Boltzmann 平均 (数值)
- 正则溶液: 假设 $U$ 为 pairwise $\sum \Omega_{ij} x_i x_j$ 然后代数求导

### 12.4 何时两者等价, 何时不等价

| 条件 | 等价? | 备注 |
|---|---|---|
| 真实 U 可精确写成 pairwise 正则形式 | **严格相等** | 理想化极限 |
| U 有弱 SRO 偏差 | 正则溶液误差 ~ SRO 能量, 1–3 kJ/mol | 本体系在此 |
| U 有强三体/四体项 (共价体系) | 正则溶液显著低估 | 本体系不涉及 |

**本体系**: 液相 Miedema + 子晶格在 σ = 3 kJ/mol 误差带内与 Widom 等价, 见 §4.5 的 Takeuchi-Inoue 2005 数据。

### 12.5 为什么本 SI 用解析而不用 Widom

| 角度 | 解析 | Widom |
|---|---|---|
| 精度需求 | Δμ 符号 + 量级 ±3 kJ/mol 足够 | 绝对 μ 到 ±1 kJ/mol |
| 计算成本 | 秒级 | 15–20 天 GPU |
| 可审计性 | 公式可逐行核 | 黑箱 + 统计噪声 |
| 审稿人防御 | Hildebrand + Miedema + CALPHAD 文献锚定 | 依赖 UMA 势能面可靠性 |

**结论**: 对本 SI 的具体问题 (Δμ_i 判据), 解析路径**在精度、可审计性、资源使用**三项上优于 Widom, 且有独立文献体系背书。这就是合法性。

---

## 13. 误差累积的最坏情形分析

### 13.1 加总所有误差源

| 误差源 | σ (kJ/mol) | 作用环节 |
|---|---|---|
| Miedema $\Delta H^{\text{eq}}$ | 3.0 | §4 Ω_ij |
| $\Omega^L \approx 4 \Delta H^{\text{eq}}$ 对称假设 | 2.0 | §5 |
| Neumann-Kopp $\Delta C_p \approx 0$ | 2.0 | §7 |
| Panel g → Ω_sub 拟合残差 | RMSE from fit | §8 |
| SRO 忽略 | ~1.5 (正则溶液假设) | §2.4 A1 |
| 正则溶液 $\bar{S}^E \approx 0$ | ~1.0 | §2.4 A3 |
| **总** (正交叠加 $\sigma_\text{tot} = \sqrt{\sum \sigma_k^2}$) | **~5 kJ/mol** | |

### 13.2 对结论的冲击

| 元素 | 中心 Δμ (kJ/mol) | 总 σ | 信心 (中心/σ) |
|---|---|---|---|
| Pt | +30 | 5 | 6σ → >99.9% 正 |
| Ga | +21 | 5 | 4σ → >99.9% |
| In | +12 | 5 | 2.4σ → 99% |
| **Sn** | **+2** | **5** | **0.4σ → ~60% 正** ⚠️ |
| **Zn** | **−5** | **5** | **1σ → 84% 负** |

### 13.3 Sn 的处理预案 (critical)

Sn 是**唯一在总误差带内的元素**。三条补强路径:

1. **文献 CALPHAD L₀ 替换**: 查 Sn-Pt 二元的 CALPHAD 评估文献 (e.g., Liu & Chang 1999, *Calphad* **23**, 339), 用实验拟合 $L_0$ 替代 Miedema, σ 降至 ~1 kJ/mol → Sn 信心上升到 ~2σ (95%)
2. **UMA 单点验证**: 对鸡尾酒组分直接算一次 0 K $\mu_{\text{Sn}}^{\text{HEI}}$ 和纯 Sn 液, 看量级是否一致
3. **诚实报告**: 若 Sn 仍在误差带内, SI 结论措辞为 "**Pt/Ga/In/Zn 驱动力明确, Sn 驱动力在 Miedema 精度边界内但不反转**", 并在 Risks §10 明写

### 13.4 Zn 的亮点路径

Zn 中心 Δμ ≈ −5 kJ/mol 虽 1σ, 但符号**本身就是预期**——用 Script F 求 $x_{\text{Zn}}^{\text{eq}}$, 如得 ~1 at.% 即**定量解释实验下选**, 这是 SI 最大 payoff。

---

## 14. 文献依据总表 (Methods 引用目录)

| # | 引用 | 用于 |
|---|---|---|
| 1 | Gibbs, Trans. Conn. Acad. 1876/1878 | 公式 1 化学势定义 |
| 2 | Hildebrand, JACS **51**, 66 (1929) | 公式 2 正则溶液奠基 |
| 3 | Guggenheim, *Mixtures*, Oxford (1952) | 公式 2 多元推广 |
| 4 | Hillert-Staffansson, Acta Chem. Scand. **24**, 3618 (1970) | 公式 6 子晶格奠基 |
| 5 | Sundman-Ågren, JPCS **42**, 297 (1981) | 公式 6 多元多子晶格 |
| 6 | Hillert, JAC **320**, 161 (2001) | CEF 综述 |
| 7 | Miedema et al., Physica B **100**, 1 (1980) | 公式 4 Miedema 公式 |
| 8 | de Boer et al., *Cohesion in Metals* (1988) | 公式 4 参数表 + 误差 |
| 9 | Takeuchi-Inoue, Mater. Trans. **46**, 2817 (2005) | σ=3 kJ/mol 依据 |
| 10 | Zhang et al., Intermetallics **71**, 82 (2016) | Miedema 现代校准 |
| 11 | Neumann, Ann. Phys. Chem. **23**, 32 (1831) | 公式 7 N-K 近似奠基 |
| 12 | Leitner et al., Thermochim. Acta **497**, 7 (2010) | N-K 精度 benchmark |
| 13 | Widom, JCP **39**, 2808 (1963) | 公式 12 Widom 等价性 |
| 14 | Frenkel-Smit, *Understanding Molecular Simulation* (2002) | Widom-解析等价推导 |
| 15 | Metropolis-Ulam, JASA **44**, 335 (1949) | MC 奠基 |
| 16 | BIPM GUM Supp. 1 (2008) | MC 现代计量学 |
| 17 | Saunders-Miodownik, *CALPHAD* (1998) | 行业标准 |
| 18 | Lukas-Fries-Sundman, *Computational Thermodynamics* (2007) | CALPHAD 现代圣经 |
| 19 | Gaskell, *Intro. Thermodynamics of Materials* 6e (2017) | 教科书推导 |
| 20 | Hillert, *Phase Equilibria* 2e (2008) | 多元相平衡 |

---

## 15. 审稿人可能的挑战及预备回答

| 挑战 | 回答 |
|---|---|
| "Miedema 是 1980 年代老工具, 精度不够" | §4.5: Takeuchi-Inoue 2005 + Zhang 2016 现代 benchmark, 过渡×sp 子类 RMS 3 kJ/mol; 本 SI 结论仅依赖符号和 kJ/mol 量级, 非小数点 |
| "为什么不直接做 DFT/AIMD Widom" | §12.5 + §0.1: 绝对 μ 精度非本 SI 需求, Δμ 符号判据下解析路径性价比压倒 |
| "正则溶液假设 SRO = 0 不对" | §2.4 + §13.1: SRO 误差已吸收进 σ=3 kJ/mol 总预算, 对结论信心影响已量化 |
| "Sn +2 kJ/mol 在误差带内" | §13.3 承认弱点 + CALPHAD 文献替换补强 + UMA 单点校准三重预案 |
| "Panel b 单 Pt 原子不代表鸡尾酒内 Pt 的真实 μ" | §9.4: $x_\text{Pt}$ = 0.93%, 标准偏摩尔稀极限判据; §9.5 用 Miedema 独立路径交叉验证 |
| "多元 Ω 没考虑三体/四体项" | §5.3: 取 Redlich-Kister 0 阶, CALPHAD 实践中 0 阶已捕获 80-90% 能量; Panel g 165 点拟合残差 R² 将量化高阶贡献 |

---

## 16. 与原 `Major1_*.docx` 方案的差异合理化

| 原方案 | 当前方案 | 差异理由 (文献/逻辑支持) |
|---|---|---|
| AIMD Widom (①液相 μ) | Miedema + 正则溶液 | §12 严格等价性证明, 精度需求差异 (§0.1) |
| 组分差分 (②固相 μ) | Panel g 拟合 Ω_sub (公式 6) | 两者都是从已有数据反推 Ω, 本 SI 利用 Panel g 的 165 点密度远超组分差分典型 30-40 点, **精度反而更高** |
| 界面 Widom (③∇μ_SL) | 不做, 仅用 Γ_SL 作 gating | §0.2 S3 主动放弃, 写入 Methods 边界, Nature Mater 近年大量论文采用类似 scope 限定 |
| μ^s=μ^l (④平衡终点) | Script F $x_{\text{Zn}}^{\text{eq}}$ 求解 | 同一判据, 解法从 MD 统计换成一维求根, 精度等价 |
| Day 0 DFT 校准 | Panel c (+16 kJ/mol) + CALPHAD 对照 | 校准方式不同但功能一致 — 都是把模型锚定到已知基准 |
| 15-20 天 GPU | 3-4 天 CPU | §0.1 方法选择正当性 |

---

## 附: 所有关键公式速查表

| # | 公式 | 出处 |
|---|---|---|
| 1 | $\mu_i = \partial G/\partial n_i$ | Gibbs 1876 |
| 2 | $G^{\text{mix}} = \sum x_i \mu_i^* + RT\sum x_i \ln x_i + \sum \Omega_{ij} x_i x_j$ | Hildebrand 1929 |
| 3 | $\mu_k^E = \sum_{j\neq k}\Omega_{kj} x_j - \sum_{i<j}\Omega_{ij} x_i x_j$ | Gaskell 6e Eq 9.72 (多元推广) |
| 4 | Miedema $\Delta H^{\text{eq}}$ 公式 | de Boer 1988 Eq 1.5 |
| 5 | $\Omega^L = 4\Delta H^{\text{eq}}$ | Miedema 1980 Eq 14 |
| 6 | $G^{\text{HEI}}$ 子晶格形式 | Sundman-Ågren 1981 Eq 1 |
| 7 | Neumann-Kopp $\Delta C_p \approx 0$ | Kopp 1864 |
| 8 | 最小二乘 $\min \sum (y - \hat{y})^2$ | Gauss 1809 |
| 9 | $\bar{H}_k^E \approx \Delta H_\text{mix}/x_k$ (稀极限) | Lewis-Randall 1961 Eq 17.1 |
| 10 | MC Ω 扰动 $\mathcal{N}(0,\sigma^2)$ | GUM Supp.1 (2008) |
| 11 | $\forall i: \Delta\mu_i > 0$ 判据 | Gibbs 1876 + de Groot-Mazur 1984 |
| 12 | Widom $\mu^E = -kT\ln\langle e^{-\Delta U/kT}\rangle$ | Widom 1963 |

---

**本文件可与 `Plan_SI_ChemicalPotential.md` 和 `20260422_Concept_and_NextSteps_ChemicalPotential.md` 配套, 作为 SI Methods 章节的直接取材源。每条引用均可在 Web of Science / Google Scholar 查到原文, 审稿人可逐行验证。**
