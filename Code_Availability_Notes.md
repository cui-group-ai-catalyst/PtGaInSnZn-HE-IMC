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

**Tier C — deferred**
The Zn eutectic-depressant down-selection (reported in the manuscript as
Supplementary Table 1) relies on external CALPHAD tables (Pandat /
Thermo-Calc) plus a Python 3.12 environment that is outside the scope of
this package. The staging script remains in `SI_Table01_ZnDownSelection/`
for transparency.

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
