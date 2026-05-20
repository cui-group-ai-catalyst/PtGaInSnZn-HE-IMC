# Panel e — script versions

> **Use only the `_v3` scripts** for paper reproduction. `v1` and `v2`
> are early prototypes whose outputs (e.g. v2 ΔG_rxn = −55.7 kJ/mol f.u.)
> do **not** match the manuscript and Supplementary Table 1.
> Authoritative reference: `script_*_v3.py` → outputs in
> `Panel_e_ChemicalPotentialCascade/outputs/*_v3_*`.

## Workflow (v3)

```
script_A_mu_liquid_v3.py   →  script_B_mu_HEI_v3.py   →  script_C_drive_force_v3.py   →  script_D_panel_a_plot_v3.py
        (μ_liquid)                    (μ_HEI)                    (F_i + ΔG_rxn)                  (cascade plot)
```

| Step | Script | Output (in `../outputs/`) | What it computes |
|---|---|---|---|
| A | `script_A_mu_liquid_v3.py` | `mu_liquid_v3_0K.csv` | SER-referenced chemical potential of each X (Ga,In,Sn,Zn) in the Ga/In/Sn/Zn liquid cocktail (Miedema binary ΔH_mix + Hildebrand–Muggianu + ΔH_fus from CRC) |
| B | `script_B_mu_HEI_v3.py` | `mu_HEI_v3_0K.csv`, `omega_beta_subl_v3_fit.csv`, `cef_fit_v3_quality.csv` | μ_i in the L1₂-Pt₃(GaInSnZn) HEI, from a Compound-Energy-Formalism (CEF) fit to **165 UMA-s-1p1 single-point energies** on Pt₃(Gaₓ In_y Sn_z Zn_w) prototypes (R² = 0.9993, RMSE = 0.14 kJ/mol). |
| C | `script_C_drive_force_v3.py` | `delta_mu_v3_0K.csv`, `delta_G_rxn_v3_summary.csv`, `panel_a_v3_tier_summary.csv` | Per-element driving force **F_i = μ_i_source − μ_i_HEI** (manuscript eq 3, SI Table 1) and total **ΔG_rxn = −150.5 kJ/mol formula unit** |
| D | `script_D_panel_a_plot_v3.py` | `panel_a_v3_schematic.png/pdf` | The three-tier chemical-potential cascade schematic (Fig. 1e) |

## Conventions

- **Composition (β-sublattice):** y_Ga = 0.65, y_In = 0.20, y_Sn = 0.10, y_Zn = 0.05.
- **Reference state:** SER (stable element reference). μ°(X) = 0 for each element in its 0 K bulk ground state.
- **Units:** kJ/mol per atom (so per-formula-unit values are 4× the per-atom values).
- **Temperature:** 0 K enthalpic terms. The dG_rxn at T > 0 K follows from `mu_*_T_table.csv`.

## Reference numbers (must match Supplementary Table 1)

| Element | y_i | μ_source (kJ/mol/atom) | μ_HEI (kJ/mol/atom) | F_i = μ_source − μ_HEI |
|---|---|---|---|---|
| Ga | 0.65 | +5.57 | −43.43 | **+49.00** |
| Zn | 0.05 | +7.30 | −31.20 | **+38.50** |
| Pt | 1.00 |  0.00 | −36.31 | **+36.31** |
| Sn | 0.10 | +6.98 | −24.23 | **+31.21** |
| In | 0.20 | +3.04 | −20.50 | **+23.54** |

Ordering: **F_Ga > F_Zn > F_Pt > F_Sn > F_In**, matching the caption of SI Table 1.

ΔG_rxn = 4·μ_Pt_HEI − Σ y_X · μ_X_source = −150.5 kJ/mol per Pt₃X formula unit (= −Σ ν_i · F_i with ν_Pt = 3, ν_X = y_X).

## Version history (kept for transparency only — do not run for paper reproduction)

- **v1** (2026-04-22, prototype). Pt SER set to 0 with HEI literature values; ordering Pt < everything else; ΔG_rxn omitted. **Do not use.**
- **v2** (2026-04-23, intermediate). Added SER for X elements (CRC ΔH_fus + Miedema Ω). HEI still from literature × 4. Ordering Sn > Ga > Pt > Zn > In; ΔG_rxn = −55.7 kJ/mol f.u. **Do not use.**
- **v3** (2026-04-23, **authoritative**). HEI replaced by CEF fit on 165 UMA single-point energies (`Panel_e/outputs/script_B_v3_meta.json` records the fit quality). Ordering **Ga > Zn > Pt > Sn > In**; ΔG_rxn = **−150.5 kJ/mol f.u.**, consistent with manuscript and SI Table 1.
- **v3 patch (2026-05-20)**. `script_C_drive_force_v3.py`: per-element F_i now computed as μ_source − μ_HEI (per manuscript eq 3 / SI Table 1), no longer as (3 μ_Pt + μ_i)/4 − μ_HEI. ΔG_rxn unchanged (−150.5 kJ/mol f.u.); per-element F_i column in `delta_mu_v3_0K.csv` now matches SI Table 1.

## Source datasets

- Panel g (now SI Fig 4) 165-point UMA scan: `Panel_e_ChemicalPotentialCascade/outputs/` (CEF fit only; raw UMA energies are bundled with SI Fig 4 in `SI_Figures/SI_Fig04_165CompositionLandscape/`).
- Reference-state energies for SER: UMA-s-1p1 single point on the experimental 0 K ground-state geometry of each element (consistent with the SGTE SSUB5 reference state).
