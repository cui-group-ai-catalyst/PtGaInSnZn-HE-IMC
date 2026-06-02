"""Verify key files and numerical anchors in the PtGaInSnZn code release.

Run from the release root:

    python scripts/verify_release.py

Or run from outside the release root:

    python scripts/verify_release.py --release-root /path/to/code_release_v2

This verifier intentionally uses only the Python standard library so that it
can run before the full scientific environment is installed. It checks bundled
CSV outputs and metadata files. It does not run UMA-dependent calculations.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def require_file(root: Path, rel: str) -> Path:
    path = root / rel
    if not path.exists():
        raise AssertionError(f"Missing required file: {rel}")
    return path


def assert_close(label: str, observed: float, expected: float, tol: float) -> None:
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=tol):
        raise AssertionError(
            f"{label}: observed {observed:.6g}, expected {expected:.6g} +/- {tol}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--release-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Path to the release root. Defaults to the parent of scripts/.",
    )
    args = parser.parse_args()
    root = args.release_root.resolve()

    required = [
        "README.md",
        "LICENSE",
        "CITATION.cff",
        "environment.yml",
        "Figure1_Code_Map.csv",
        "Code_Availability_Notes.md",
        "Third_Party_Model_and_Data_Notes.md",
        "UMA_Checkpoint_Setup.md",
        "NATURE_SOFTWARE_CHECKLIST.md",
        "REVIEWER_REPRODUCTION_GUIDE.md",
        "demo/README.md",
        "RELEASE_NOTES_v1.3.0.md",
        "scripts/recompute_element_referenced_Hf.py",
        "shared/UMA_Element_Reference_Energies.csv",
        "shared/element_reference_structures/README.md",
        "shared/element_reference_structures/element_reference_manifest.csv",
        "shared/element_reference_structures/20260417_Ga_mp142_reference.cif",
        "shared/element_reference_structures/20260417_In_bct_reference.cif",
        "shared/element_reference_structures/20260417_Pt_mp126_reference.cif",
        "shared/element_reference_structures/20260417_Sn_alpha_reference.cif",
        "shared/element_reference_structures/20260417_Zn_hcp_reference.cif",
        "Panel_a_BinaryHeatmap/data_FigA_v2_FamilyOrdered_Origin.csv",
        "Panel_b_MultiComponentScatter/data_FigB_central_anchors.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_Summary.csv",
        "Panel_e_ChemicalPotentialCascade/outputs/delta_G_rxn_v3_summary.csv",
        "Panel_f_Wetting/data_FigF_Wetting_Ranked.csv",
    ]

    for rel in required:
        require_file(root, rel)

    fig_a = read_csv(root / "Panel_a_BinaryHeatmap/data_FigA_v2_FamilyOrdered_Origin.csv")
    ga_row = next(row for row in fig_a if row["Target_Element"] == "Ga")
    assert_close("Fig. 1a Pt-Ga dH_mix", float(ga_row["Pt"]), -32.06, 0.02)
    fig_a_regen = root / "Panel_a_BinaryHeatmap/data_FigA_v2_FamilyOrdered_Origin_regen.csv"
    if fig_a_regen.exists():
        regen_rows = read_csv(fig_a_regen)
        regen_ga = next(row for row in regen_rows if row["Target_Element"] == "Ga")
        assert_close("Fig. 1a regen Pt-Ga dH_mix", float(regen_ga["Pt"]), -32.06, 0.02)

    fig_b = read_csv(root / "Panel_b_MultiComponentScatter/data_FigB_central_anchors.csv")
    pt_b = next(row for row in fig_b if row["Host"] == "Pt")
    if int(pt_b["dH_rank"]) != 1:
        raise AssertionError("Fig. 1b central check: Pt is not ranked #1")
    assert_close("Fig. 1b Pt dH_mix", float(pt_b["dH_mix_kJmol"]), -24.440, 0.005)
    assert_close("Fig. 1b Pt dG_mix 500 K", float(pt_b["dG_mix_500K_kJmol"]), -29.839, 0.005)
    fig_b_regen = root / "Panel_b_MultiComponentScatter/data_FigB_central_anchors_regen.csv"
    if fig_b_regen.exists():
        regen_b = read_csv(fig_b_regen)
        regen_pt = next(row for row in regen_b if row["Host"] == "Pt")
        if int(regen_pt["dH_rank"]) != 1:
            raise AssertionError("Fig. 1b regen central check: Pt is not ranked #1")
        assert_close("Fig. 1b regen Pt dH_mix", float(regen_pt["dH_mix_kJmol"]), -24.440, 0.005)
        assert_close("Fig. 1b regen Pt dG_mix 500 K", float(regen_pt["dG_mix_500K_kJmol"]), -29.839, 0.005)

    fig_c = read_csv(root / "Panel_c_OrderedVsDisordered/data_FigC_Summary.csv")
    metrics = {row["Metric"]: row["Value"] for row in fig_c}
    assert_close("Fig. 1c ordering gap", float(metrics["Ordering_Gap_kJ_mol"]), 16.0436, 0.001)
    if int(float(metrics["Disordered_N"])) != 30:
        raise AssertionError("Fig. 1c disordered ensemble size is not 30")

    fig_e = read_csv(
        root / "Panel_e_ChemicalPotentialCascade/outputs/delta_G_rxn_v3_summary.csv"
    )
    assert_close(
        "Fig. 1e delta_G_rxn per f.u.",
        float(fig_e[0]["delta_G_rxn_per_fu_kJmol"]),
        -150.5354,
        0.001,
    )

    fig_f = read_csv(root / "Panel_f_Wetting/data_FigF_Wetting_Ranked.csv")
    pt_f = next(row for row in fig_f if row["Host"] == "Pt")
    assert_close("Fig. 1f Pt gamma_SL", float(pt_f["Gamma_SL"]), -0.44834, 0.0001)
    if pt_f["Status"] != "Wets":
        raise AssertionError("Fig. 1f Pt wetting status is not Wets")

    readme = (root / "README.md").read_text(encoding="utf-8")
    for text in [
        "https://github.com/cui-group-ai-catalyst/PtGaInSnZn-HE-IMC",
        "10.5281/zenodo.20111606",
    ]:
        if text not in readme:
            raise AssertionError(f"README does not contain required identifier: {text}")

    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    for forbidden in ["TODO", "TO BE FILLED", "<DOI", "<JOURNAL", "<MANUSCRIPT"]:
        if forbidden in citation:
            raise AssertionError(f"CITATION.cff contains unresolved placeholder: {forbidden}")

    ref_rows = read_csv(root / "shared/UMA_Element_Reference_Energies.csv")
    ref_elems = {row["Element"] for row in ref_rows}
    if ref_elems != {"Pt", "Ga", "In", "Sn", "Zn"}:
        raise AssertionError(f"Unexpected element-reference set: {sorted(ref_elems)}")
    manifest = read_csv(root / "shared/element_reference_structures/element_reference_manifest.csv")
    for row in manifest:
        require_file(root, row["Structure_File"])

    print(f"Release root: {root}")
    print("All release verification checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
