# Code and Data Availability — Reproducibility Notes

## Scope

The scripts in this package reproduce the figures and the numerical results
of Figure 1 (panels a–f) and the Supplementary Information of the
manuscript. Reference outputs (CSV) and preview PNGs are included so that
each panel can be visually compared against the manuscript before any
recomputation is attempted.

## Reproducibility tiers

**Tier A — fully self-contained, analytical**
Panels a, b, f and Supplementary Fig. 2: Miedema-derived. Reproducible from
this package alone with the conda environment `environment.yml`. No external
model weights required.

**Tier B — requires Fairchem UMA-s-1p1**
Panels c, d, e and Supplementary Figs. 3, 4: depend on UMA-s-1p1 single-point
energies. The UMA checkpoint is provisioned per `UMA_Checkpoint_Setup.md`.
Reference UMA-derived CSVs are bundled so reviewers can reproduce the
post-processing without re-running UMA, if desired.

**Tier C — Miedema + Hildebrand-Muggianu (literature-anchored)**
The Zn eutectic-depressant liquidus survey (reported in the manuscript as
**Supplementary Tables 2 and 3**) is reproducible from the bundled
scripts in `SI_Figures/SI_Table01_ZnDownSelection/` (the folder name is
historical; current SI numbering puts this content in Tables 2 & 3):
- `script_SI_LiquidusPredictor.py` extrapolates the six liquid binary
  enthalpies (Witusiewicz CALPHAD assessments) to the Ga-In-Sn-Zn
  quaternary via the Muggianu form, then fits a 2nd-order
  T_l(ΔH_mix) regression to three literature anchors (Daeneke 2018 /
  Wu 2025 / Shentu 2023). Because the fit has 3 parameters and 3
  anchors, in-regime residuals are zero by construction; Bai 2022 is
  documented as the only out-of-regime test point. Predicts a
  eutectic minimum at ~0.8 at% Zn, consistent with the 1 at%
  operational target.
- Pure scientific Python (numpy, pandas, matplotlib, scipy); no CALPHAD
  package required. Outputs `data_SI_Liquidus_Validation.csv`,
  `data_SI_Liquidus_Scan.csv`, `preview_SI_Liquidus.png`,
  `notes_SI_Liquidus_Calibration.txt`.

**Tier D — deferred**
Full CALPHAD treatment of Stage C in `script_SI_DownSelectionFunnel.py`
(pycalphad + COST507 Ga-In-Sn-Zn TDB) is still a `pending-CALPHAD`
placeholder. The four-stage funnel script falls through stages A → B → D
without it; the bundled `script_SI_LiquidusPredictor.py` covers the
liquidus prediction needs for Supplementary Tables 2 and 3 without
requiring CALPHAD.

## What is fixed and what is not

- The Miedema atomic-parameter table (`shared/data_periodic_table.py`)
  is the exact local snapshot used for the manuscript values, derived from
  Takeuchi & Inoue, *Materials Transactions* **46**, 2817 (2005), with
  element-specific molar volumes from de Boer *et al.*, *Cohesion in
  Metals*, North-Holland (1988).
- The UMA checkpoint hash is pinned (`uma-s-1p1.pt`, MD5
  `36a2f071350be0ee4c15e7ebdd16dde1`, snapshot `38529caa`).
- Random-structure generation uses fixed seeds (100–129 in `Panel_c`),
  so re-running yields bit-identical structures.
- The 165-prototype B-sublattice scan is provided as a precomputed CSV
  rather than a recipe; the structure-generation pipeline is documented
  in `Panel_e_ChemicalPotentialCascade/docs_ChemicalPotential_Workflow.md`.

## Numerical conventions

- Energies are reported in **kJ mol⁻¹ atom⁻¹** unless the formula unit
  basis is explicitly stated. The conversion factor used throughout is
  96.485 kJ mol⁻¹ eV⁻¹.
- Formation enthalpies are referenced to the **stable element reference
  (SER)**: the most stable solid phase of each constituent at 0 K is
  assigned μ° = 0. Element reference energies are tabulated in
  `Panel_c_OrderedVsDisordered/data_FigC_Summary.csv` and computed by the
  same UMA-s-1p1 single-point pipeline used for the alloy structures.
- Sample standard deviation (`ddof=1`, the numpy default) is used wherever
  a 1σ variation across the 30 random HEA structures is reported. The
  `mean ± SEM` values quoted in the SI use SEM = σ_sample / √N with N = 30,
  giving the SI's reported `σ = 2.26` and gap-SEM `= 0.41 kJ mol⁻¹ atom⁻¹`.

## Excluded from this release

- Figure-aesthetic Origin templates (`.opju`)
- Manuscript-only figure-caption Word documents
- Internal review memos and draft scripts

These can be supplied on request but are not required for numerical
reproducibility.

## Citation

When reusing any panel-specific script or reference output, please cite the
manuscript and the dependent third-party packages
(`Third_Party_Model_and_Data_Notes.md` not included in this package; see
`UMA_Checkpoint_Setup.md` for the UMA citation).
