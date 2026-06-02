# Release notes - v1.3.0 Nature software checklist release

This proposed release upgrades the existing code package for Nature
editor/reviewer review. It does not change the manuscript numerical claims.
It adds the documentation, demo, verification, and metadata clarity requested
by the Nature Research software submission checklist.

## Added

- `NATURE_SOFTWARE_CHECKLIST.md`
  Requirement-by-requirement mapping to Nature/PDF checklist items.
- `demo/README.md`
  CPU-only reviewer demo with commands, expected outputs, and expected runtime.
- `scripts/verify_release.py`
  Lightweight release verifier for key files and numerical anchors.
- Expanded README sections:
  - system requirements;
  - installation;
  - demo;
  - expected output;
  - expected runtime;
  - instructions for use;
  - manuscript algorithm location.

## Changed

- Clarifies GitHub/Zenodo roles:
  - concept DOI for the release family;
  - exact version DOI for the submitted archive.
- Promotes UMA checkpoint restrictions to the main README.
- Recommends `v1.3.0` as the Nature-compliance release target.

## To complete before publication

- Fill the exact Zenodo DOI minted for `v1.3.0`.
- Add the exact version DOI to `CITATION.cff` after Zenodo mints it.
- Update `LICENSE` copyright holder.
- Resolve ORCID/affiliation fields if included.
- Confirm double-blind peer-review handling.
