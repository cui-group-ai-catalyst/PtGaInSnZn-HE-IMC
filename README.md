# Code Release — Pt₃(Ga,In,Sn,Zn) High-Entropy Intermetallic Assembly

This package accompanies the manuscript *"Data-Driven Liquid-Metal-Mediated
Chemical-Potential-Driven Assembly of High-Entropy Intermetallic Nanocrystals"*
and reproduces the calculations underlying Figure 1 (panels a–f) and the
related Supplementary Information.

## Code availability

The source code in this repository is openly available at
`https://github.com/xiuxiudebobo/PtGaInSnZn-HE-IMC` and archived with a
persistent identifier at `https://doi.org/10.5281/zenodo.20111607`. All in-house code is released under the
**MIT License** (see `LICENSE`). Third-party model weights and datasets are
not redistributed; their sources, versions, and licenses are documented in
`Third_Party_Model_and_Data_Notes.md`. The pretrained UMA-s-1p1 checkpoint
used in this study is provisioned per `UMA_Checkpoint_Setup.md`.

## Layout

```
code_release_v2/
├── README.md                    # this file
├── Figure1_Code_Map.csv         # panel ↔ script ↔ output mapping
├── environment.yml              # Python environment (Python 3.10 + UMA)
├── Code_Availability_Notes.md   # reproducibility scope
├── UMA_Checkpoint_Setup.md      # UMA-s-1p1 checkpoint provisioning
│
├── Panel_a_BinaryHeatmap/                    # Fig 1a — 20-host binary ΔH_mix
├── Panel_b_MultiComponentScatter/            # Fig 1b — 20-host cocktail scatter
├── Panel_c_OrderedVsDisordered/              # Fig 1c — UMA HEI vs HEA (0 K)
├── Panel_d_GibbsCurveVsT/                    # Fig 1d — HEI vs HEA Gibbs vs T
├── Panel_e_ChemicalPotentialCascade/         # Fig 1e — per-element μ cascade
├── Panel_f_Wetting/                          # Fig 1f — γ_SL across 28 hosts
│
├── SI_Figures/
│   ├── SI_Fig02_SizeMismatch/                # Supplementary Fig. 2  — Hume-Rothery δ vs Miedema ΔH_mix, 28 hosts
│   ├── SI_Fig03_TripleConsensus/             # Supplementary Fig. 3  — three-method consensus, 12 M-Ga binaries
│   ├── SI_Fig04_165CompositionLandscape/     # Supplementary Fig. 4  — 165 B-sublattice prototypes, ΔH_f vs Ga at%
│   └── SI_Table01_ZnDownSelection/           # Supplementary Table 1 — Zn eutectic-depressant funnel (deferred figure)
│
├── SI_Note_S1_ChemicalPotential/             # full mathematical derivation
└── shared/
    └── data_periodic_table.py                # Miedema atomic parameters
```

## Reproducibility scope

This package reproduces the **numerical values and reference outputs** of
each Figure 1 panel using the exact scripts that produced the manuscript-
reported results. Two reproducibility levels are documented:

- **Level A — analytical/Miedema panels (a, b, f, SI Fig 2, SI Fig 3 Miedema column)**.
  Self-contained; reproducible from this package alone using a standard
  Python scientific stack (`environment.yml`).

- **Level B — UMA-based panels (c, d, e, SI Fig 4, SI Fig 3 UMA column)**. Requires
  the Fairchem UMA-s-1p1 checkpoint, provisioned per `UMA_Checkpoint_Setup.md`.
  The `fairchem-core==2.13.0` package handles the model loading.

## Quick start

```bash
conda env create -f environment.yml
conda activate ptgainsnzn

# Analytical panels (no UMA required)
python Panel_a_BinaryHeatmap/script_FigA_v2_family_ordered.py
python Panel_b_MultiComponentScatter/script_FigB_v2_family_ordered.py
python Panel_f_Wetting/script_FigF_Wetting.py

# UMA-based panels (set checkpoint per UMA_Checkpoint_Setup.md)
export UMA_CHECKPOINT_PATH=/path/to/uma-s-1p1.pt
python Panel_c_OrderedVsDisordered/script_FigC_OriginReady_ElementRef.py
python Panel_c_OrderedVsDisordered/script_FigC_v2_ElementRef_Preview.py
python Panel_d_GibbsCurveVsT/script_FigD_GibbsCurve.py
```

## Panel summaries

| Panel | What is computed | Method | Reference output |
|---|---|---|---|
| a | Binary ΔH_mix for 20 hosts vs 6 liquid partners | Miedema (de Boer 1988) | `data_FigA_v2_FamilyOrdered_Origin.csv` |
| b | Multi-component cocktail ΔH_mix and ΔG_mix at 500 K | Miedema + regular solution | `data_FigB_v2_Origin_dH.csv`, `data_FigB_v2_Origin_dG.csv` |
| c | Ordered L1₂ vs disordered random ΔH_f at 0 K | UMA-s-1p1 single-point | `data_FigC_Summary.csv` (HEI=−30.01, HEA=−13.96±2.26, gap=16.04) |
| d | HEI vs HEA Gibbs curves vs T (473–2500 K) | UMA c-data + −T·ΔS_config | `data_FigD_GibbsCurve.csv` |
| e | Per-element chemical potential cascade | Hildebrand RS + CEF on 165 UMA points | `outputs/delta_mu_v3_0K.csv` (ΔG_rxn=−150.5 kJ/mol f.u.) |
| f | Solid-liquid γ_SL for 28 hosts | Macroscopic-atom Miedema | `data_FigF_Wetting_Ranked.csv` (Pt=−0.50, 13 wet) |

## Key verified numbers

Each entry below is taken from a specific reference CSV in this package, which
is the canonical numerical source. The Supplementary Information may aggregate
some of these at a higher level (for example, the ordering gap is reported in
the SI as `16.04 ± 0.41 kJ mol⁻¹ atom⁻¹ (mean ± SEM, σ = 2.26, z = 38.9)`,
which combines the per-config HEI and HEA values listed below).

| Quantity | Value | Source |
|---|---|---|
| Pt binary ΔH_mix (Pt-Ga) | −33.76 kJ mol⁻¹ | `Panel_a_BinaryHeatmap/data_FigA_v2_FamilyOrdered_Origin.csv` |
| HEI ΔH_f (0 K, element ref) | −30.01 kJ mol⁻¹ atom⁻¹ | `Panel_c_OrderedVsDisordered/data_FigC_Summary.csv` |
| HEA ΔH_f mean ± σ | −13.96 ± 2.26 kJ mol⁻¹ atom⁻¹ | same |
| Ordering enthalpy gap | 16.04 kJ mol⁻¹ atom⁻¹ | same |
| T* crossover (frozen β) | 2123 K | `Panel_d_GibbsCurveVsT/notes_FigD_GibbsCurve.md` |
| T* crossover (mixed β) | 3429 K | same |
| ΔG_rxn at 0 K | −150.5 kJ mol⁻¹ f.u. | `Panel_e_ChemicalPotentialCascade/outputs/delta_G_rxn_v3_summary.csv` |
| γ_SL (Pt) | −0.50 J m⁻² | `Panel_f_Wetting/data_FigF_Wetting_Ranked.csv` |
| Hosts with γ_SL < 0 | 13 of 28 | same |

## Limitations and known gaps

- `SI_Table01_ZnDownSelection/` contains the script but no precomputed CSV/PNG;
  the script depends on a Python 3.12 environment and CALPHAD lookup tables
  not yet bundled in this release. The Zn down-selection is reported in the
  manuscript as **Supplementary Table 1** (liquidus/solidus of Ga-In-Sn-Zn vs
  Zn loading) rather than a dedicated figure, and is not required to reproduce
  the main-text numerical conclusions.
- Some UMA-derived numbers in the chemical potential pipeline rely on a
  prior precomputed dataset of 165 B-sublattice configurations
  (`Panel_e_ChemicalPotentialCascade/outputs/`). Re-running these from
  scratch requires the UMA checkpoint and ~1 hr of single-GPU computation.

## License and attribution

See `Code_Availability_Notes.md` for licensing and `UMA_Checkpoint_Setup.md`
for third-party model attribution.
