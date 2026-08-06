# Reviewer reproduction guide

This guide is the shortest path for an editor or reviewer who is not familiar
with this code base. It explains what can be reproduced immediately from the
released files, what requires external access, and what output should be seen.

## 1. What this package proves

The package has two reproducibility levels:

- **Immediate CPU demo**: verifies the self-contained analytical calculations
  behind the Figure 1a Miedema heatmap and Figure 1b central cocktail ranking,
  then checks key bundled reference files and numerical anchors.
- **Panel-level post-processing**: regenerates reviewer-facing `_regen` CSV/PNG
  outputs from the bundled reference data for Figure 1a-f and Supplementary
  Figs. 2-4 / Supplementary Tables 2-3.
- **From-scratch UMA energy calculations**: require the gated UMA-s-1p1
  checkpoint. Panel c and Supplementary Fig. 3 now bundle their deterministic
  structure-generation/CIF inputs and write machine-readable validation
  reports. Other UMA-derived panels may remain table-level reproductions; see
  `UMA_Checkpoint_Setup.md` and `Code_Availability_Notes.md`.

## 2. Install

Install Miniconda or Anaconda, then open a terminal in the repository root.
On Windows, Anaconda Prompt is usually the easiest option.

```bash
conda env create -f environment.yml
conda activate ptgainsnzn
```

Expected install time is approximately 10-25 minutes for the CPU scientific
stack and 20-45 minutes for the full environment including ML packages. Network
speed and conda solver performance dominate the total time.

## 3. Run the quick reviewer demo

```bash
python Panel_a_BinaryHeatmap/regen_FigA_data.py
python Panel_b_MultiComponentScatter/verify_FigB_central.py
python scripts/verify_release.py
```

The demo should complete in less than 5 minutes after installation. On a normal
workstation, the CPU calculations themselves usually complete in less than
1 minute.

## 4. Expected success criteria

The run is successful if:

- `Panel_a_BinaryHeatmap/data_FigA_v2_FamilyOrdered_Origin_regen.csv` is
  written and the Pt-Ga value is `-32.06 kJ/mol`.
- `Panel_b_MultiComponentScatter/data_FigB_central_anchors_regen.csv` is
  written and Pt has `dH_rank = 1`.
- The Panel b Pt values are `dH_mix_kJmol = -24.440` and
  `dG_mix_500K_kJmol = -29.839`.
- `python scripts/verify_release.py` prints
  `All release verification checks passed.`

## 5. Optional panel-level post-processing

The full set of released panel scripts is listed in `Figure1_Code_Map.csv`.
Scripts marked `No` under `UMA_Required` are analytical or post-processing
scripts that do not require the UMA checkpoint. Scripts marked `Yes` or
`Indirect` use bundled UMA-derived reference CSVs in the default release path;
from-scratch UMA energy re-runs require the external checkpoint and may require
additional non-redistributed structure inputs.

The release audit ran the default scripts for Figure 1a-f, Supplementary
Figs. 2-4, and Supplementary Tables 2-3 from the bundled inputs. All default
scripts completed successfully in the tested local environment.

The bundled 30+30 statistics and Panel d entropy model can be regenerated
without loading UMA:

```bash
python Panel_c_OrderedVsDisordered/script_FigC_OrderedEnsemble_Analysis.py
python Panel_c_OrderedVsDisordered/script_FigC_OrderedAll68_Sensitivity.py
python Panel_d_GibbsCurveVsT/script_FigD_GibbsCurve_Ensemble.py
python scripts/verify_release.py
```

To regenerate the ordered energies themselves, obtain the UMA checkpoint and
run `script_FigC_OrderedEnsemble_Rerun_UMA.py --device cuda`. The optional
`--n-ordered 68` path evaluates every symmetry-distinct B-sublattice class.

## 6. Which preview image should I inspect?

Some folders contain both a canonical preview image and a `_regen` preview
image.

- `preview_FigX_*.png` without `_regen` is the manuscript-facing canonical
  preview shipped with the release.
- `preview_FigX_*_regen.png` is produced by re-running the release script. It is
  kept under a separate filename so the canonical preview is not overwritten.

The Nature software checklist does not require two preview images. One demo
output would be enough. This package keeps both because it makes the audit more
transparent: reviewers can inspect the canonical manuscript-facing preview and
also see what the released script regenerates on their machine. If a reviewer
only wants one image, use the non-`_regen` canonical preview listed in
`Figure1_Code_Map.csv`.

## 7. Data and model provenance quick map

Use this map when checking where each result comes from:

- **Figure 1a**: corrected Miedema binary enthalpies from
  `shared/data_periodic_table.py`; no external model.
- **Figure 1b**: Miedema plus Hildebrand-Muggianu regular-solution mixing at
  500 K; no external model.
- **Figure 1c**: historical 1+30 regression plus measured 30+30 ensemble and
  exhaustive 68 ordered-class sensitivity artifacts. UMA reruns require the
  separately obtained checkpoint.
- **Figure 1d**: uses the measured 30 ordered and original 30 disordered
  Figure 1c enthalpies plus ideal configurational-entropy bounds; it does not
  load a new atomistic model.
- **Figure 1e**: liquid term from Hildebrand/Miedema conventions plus a CEF fit
  to the bundled 165-point UMA landscape; default script reads bundled outputs.
- **Figure 1f**: macroscopic-atom/Miedema solid-liquid interfacial-energy
  calculation; no external model.
- **Supplementary Fig. 2**: atomic-size mismatch plus corrected Miedema drive;
  no external model.
- **Supplementary Fig. 3**: bundled Materials Project, UMA and CHGNet values;
  default mode re-plots bundled data. `--rerun-uma` now evaluates the bundled
  15 binary and 16 elemental-reference CIFs, writes separate `_uma_regen`
  outputs, and checks them against the canonical table.
- **Supplementary Fig. 4**: bundled 165-point UMA-derived composition landscape;
  default script post-processes the bundled CSV.
- **Supplementary Tables 2-3**: literature-anchored liquidus predictor using
  Miedema/Hildebrand-Muggianu-style binary enthalpy inputs and calibration
  anchors documented in the script.

Detailed third-party sources, licenses and citations are listed in
`Third_Party_Model_and_Data_Notes.md`.

## 8. Safe use with AI or vibe coding

If using an AI coding assistant, give it these files first:

```text
README.md
REVIEWER_REPRODUCTION_GUIDE.md
Figure1_Code_Map.csv
Third_Party_Model_and_Data_Notes.md
Code_Availability_Notes.md
```

Recommended prompt:

```text
I am reviewing this Nature software package. Do not overwrite canonical files.
Use Figure1_Code_Map.csv to identify each script and reference output. Run the
quick reviewer demo first, then run scripts/verify_release.py. If you regenerate
panel outputs, write or inspect only *_regen files unless a script is explicitly
documented otherwise. Explain which results are self-contained and which depend
on bundled UMA-derived data or external checkpoint access.
```

Safety rules:

- Do not delete or overwrite the non-`_regen` canonical CSV/PNG files.
- Do not commit or upload UMA checkpoint files (`*.pt`, `*.pth`, `*.ckpt`,
  `*.bin`, `*.safetensors`).
- After any edit or regeneration, run `python scripts/verify_release.py`.
- Treat `Third_Party_Model_and_Data_Notes.md` as the source of truth for
  external model/data attribution.

## 9. Bounded validation report

Reviewers who prefer a single entry point can run:

```bash
python experimental_extensions/run_validation.py
```

This validates the external system manifest, reruns the complete-manifold CEF
fit, non-endmember LOOCV, and Ga-count group holdout, and compares Materials
Project DFT, UMA-s-1p1, and CHGNet on the bundled M-Ga binary reference set.
The comparison exports pairwise RMSE, Spearman correlation, top-k overlap, and
every ranking reversal. It then writes:

```text
experimental_extensions/outputs/system_validation/validation_results.json
experimental_extensions/outputs/system_validation/validation_evidence.png
experimental_extensions/outputs/system_validation/validation_evidence.pdf
experimental_extensions/outputs/system_validation/evidence_report.html
experimental_extensions/outputs/system_validation/evidence_report.pdf
```

The report intentionally states that new-host/prototype and N/O/S/P-compound
transferability were not evaluated and that synthesizability is not predicted.
The generic formulas and metric definitions are in
`experimental_extensions/FORMULAS_AND_VALIDATION.md`.

## 10. Troubleshooting

- If `python` is not found, activate the conda environment first:
  `conda activate ptgainsnzn`.
- If Greek letters or dashes display strangely in a Windows terminal, rely on
  the CSV files and the verifier output. This is a terminal encoding display
  issue, not a numerical difference.
- If a UMA checkpoint is requested, follow `UMA_Checkpoint_Setup.md`. The quick
  reviewer demo does not require UMA, Hugging Face access, or GPU hardware.
