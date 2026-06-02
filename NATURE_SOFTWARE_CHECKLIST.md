# Nature software submission checklist mapping

This file maps the Nature Research software submission checklist and Nature
Portfolio reporting standards to concrete files in this repository.

## Repository identifiers

- GitHub repository: `https://github.com/xiuxiudebobo/PtGaInSnZn-HE-IMC`
- Zenodo concept DOI: `https://doi.org/10.5281/zenodo.20111606`
- Exact submitted version DOI: pending until the `v1.3.0` GitHub release is
  archived by Zenodo
- Proposed release tag: `v1.3.0`
- License: MIT

## PDF checklist mapping

| PDF checklist item | Location in this package | Status |
|---|---|---:|
| Compiled standalone software and/or source code | Panel folders, `SI_Figures/`, `shared/` | Yes |
| Small dataset to demo software/code | `demo/README.md`; `REVIEWER_REPRODUCTION_GUIDE.md`; selected bundled CSVs in Panel a and Panel b | Yes |
| README includes system requirements | `README.md`, `environment.yml` | Yes |
| Dependencies and versions | `environment.yml`; README dependency summary | Yes |
| Operating systems and versions | `README.md` system requirements | Yes |
| Tested versions | `README.md`; `environment.yml`; release notes | Yes |
| Required non-standard hardware | `README.md`; `UMA_Checkpoint_Setup.md`; `environment.yml` | Yes |
| Installation guide | `README.md` installation section | Yes |
| Typical install time | `README.md` installation section | Yes |
| Demo instructions | `README.md` demo; `demo/README.md`; `REVIEWER_REPRODUCTION_GUIDE.md` | Yes |
| Expected output | `README.md` demo; `demo/README.md`; `REVIEWER_REPRODUCTION_GUIDE.md`; `scripts/verify_release.py` | Yes |
| Expected demo runtime | `README.md` demo; `demo/README.md`; `REVIEWER_REPRODUCTION_GUIDE.md` | Yes |
| Instructions for use | `README.md` instructions for use; `README.md` section "Using your own data" covers 3 substitution scenarios with worked examples; `scripts/recompute_element_referenced_Hf.py` is the turnkey utility for swapping the 5 element reference energies; `shared/element_reference_structures/*.cif` is bundled for fair-comparison re-runs under a different ML potential; `REVIEWER_REPRODUCTION_GUIDE.md` | Yes |
| Reproduction instructions | `README.md`; `Code_Availability_Notes.md`; `Figure1_Code_Map.csv` | Yes |
| License description | `LICENSE`; README code availability section | Yes |
| Open-source repository link | README and this checklist | Yes |
| Manuscript algorithm/pseudocode location | README section "Where the algorithm is described" | Yes |

## Nature Portfolio reporting standards mapping

| Requirement | Location | Notes |
|---|---|---|
| Minimum dataset transparent to readers | Bundled reference CSVs and paper Source Data | Derived CSVs are included for computational reproduction |
| Public identifiers for code | GitHub URL and Zenodo DOI | Concept DOI present; exact version DOI to add after archive |
| Restrictions disclosed | `UMA_Checkpoint_Setup.md`, `Third_Party_Model_and_Data_Notes.md` | UMA checkpoint is gated and not redistributed |
| Code availability section in manuscript | Manuscript Data/Code availability statement | Existing extraction includes GitHub, Zenodo, and UMA HF URL |
| DOI-minting repository | Zenodo | Use exact version DOI in reference list if required |
| OSI-approved license | MIT license | Copyright holder currently listed as `Bo Xu and co-authors` |
| Version management | Git tags and GitHub releases | Proposed `v1.3.0` tag; not yet committed/tagged |

## Quick reviewer path

```bash
conda env create -f environment.yml
conda activate ptgainsnzn
python Panel_a_BinaryHeatmap/regen_FigA_data.py
python Panel_b_MultiComponentScatter/verify_FigB_central.py
python scripts/verify_release.py
```

Expected total runtime for the quick reviewer path: less than 5 minutes after
environment installation.

For reviewers unfamiliar with the code base, use
`REVIEWER_REPRODUCTION_GUIDE.md` as the primary step-by-step entry point.

## Full reproduction path

1. Run quick reviewer path.
2. For UMA-derived panels, first use the bundled reference CSVs for
   post-processing and comparison.
3. Configure the UMA checkpoint using `UMA_Checkpoint_Setup.md` only if
   re-running UMA energy calculations is required.
4. Compare outputs against reference CSVs and preview PNGs.

Full reproduction time depends on whether the reviewer only post-processes
bundled UMA-derived CSVs or re-runs UMA energy calculations. The precomputed
UMA-derived CSVs are bundled so reviewers can inspect numerical results without
rerunning all UMA calculations.

## Known release-completeness limitations

Disclosed up front so reviewers do not encounter them as surprises:

- **SI Fig 2 — partial reproducer.** The manuscript canonical figure
  `preview_FigE_ThreeWay.png` is an Origin layout built from
  `data_FigE_True_ThreeWay.csv`. The in-release script reproduces only the
  ranked-resistance view (`preview_FigE_Resistance_Plot_regen.png`), not the
  three-way composition. The ThreeWay source CSV is bundled; the ThreeWay
  rendering script is not.
- **`_regen` file convention.** Every panel ships canonical and `_regen`
  parallel files; the `_regen` files are what the in-release scripts produce
  on a fresh run. See README section "Canonical files vs `_regen` files" for
  the convention, the bit-identity guarantee on the CSV side, and the
  cosmetic-only PNG drift for Panel f, SI Fig 2, and SI Fig 3 `_regen`.
- **UMA scope.** The release includes UMA-derived reference CSVs for reviewer
  inspection and post-processing. Re-running UMA energy calculations requires
  the gated UMA-s-1p1 checkpoint with MD5
  `36a2f071350be0ee4c15e7ebdd16dde1`; some archival from-scratch branches also
  require local structure inputs that are not redistributed.

## Items requiring author input before final public release

- Create the `v1.3.0` GitHub release and fill the exact Zenodo version DOI
  after Zenodo archives it.
- Add ORCID/affiliation if the author wants them in `CITATION.cff`.
- Replace `Bo Xu and co-authors` with the final formal author list if required
  by the institution or journal production team.
- Confirm whether the submission uses double-blind peer review. If yes, prepare
  an anonymized review package.
