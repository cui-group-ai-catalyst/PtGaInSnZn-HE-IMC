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
import hashlib
import json
import math
import statistics
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


def classify_gamma_sl(gamma_sl: float, uncertainty: float = 0.20) -> str:
    if gamma_sl + uncertainty < 0:
        return "favourable"
    if gamma_sl - uncertainty > 0:
        return "unfavourable"
    return "indeterminate"


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
        "demo/run_reviewer_demo.py",
        "demo/tests/test_reviewer_demo.py",
        "demo/examples/three_element/README.md",
        "demo/examples/three_element/manifold_energies.csv",
        "demo/examples/three_element/backend_energies.csv",
        "demo/examples/three_element/manifold_config.json",
        "demo/examples/three_element/system_manifest.json",
        "RELEASE_NOTES_v1.3.2.md",
        "RELEASE_NOTES_v1.4.0.md",
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
        "Panel_b_MultiComponentScatter/data_FigB_composition_sensitivity_regen.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_Summary.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_Raw_UMA_Energies.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_Raw_UMA_Energies_regen.csv",
        "Panel_c_OrderedVsDisordered/script_FigC_Rerun_UMA.py",
        "Panel_c_OrderedVsDisordered/validation_FigC_UMA_Rerun.json",
        "Panel_c_OrderedVsDisordered/script_FigC_OrderedEnsemble_Rerun_UMA.py",
        "Panel_c_OrderedVsDisordered/script_FigC_OrderedEnsemble_Analysis.py",
        "Panel_c_OrderedVsDisordered/data_FigC_OrderedEnsemble_Raw_UMA_Energies_regen.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_OrderedEnsemble_StructureManifest_regen.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_OrderedVsDisordered_Ensembles_regen.csv",
        "Panel_c_OrderedVsDisordered/validation_FigC_OrderedEnsemble_UMA_regen.json",
        "Panel_c_OrderedVsDisordered/validation_FigC_OrderedVsDisordered_Ensembles_regen.json",
        "Panel_c_OrderedVsDisordered/script_FigC_OrderedAll68_Sensitivity.py",
        "Panel_c_OrderedVsDisordered/data_FigC_OrderedAll68_Raw_UMA_Energies_sensitivity.csv",
        "Panel_c_OrderedVsDisordered/data_FigC_OrderedAll68_StructureManifest_sensitivity.csv",
        "Panel_c_OrderedVsDisordered/validation_FigC_OrderedAll68_UMA_sensitivity.json",
        "Panel_c_OrderedVsDisordered/validation_FigC_OrderedAll68_Sensitivity.json",
        "Panel_d_GibbsCurveVsT/script_FigD_GibbsCurve_Ensemble.py",
        "Panel_d_GibbsCurveVsT/data_FigD_GibbsCurve_Ensemble_regen.csv",
        "Panel_d_GibbsCurveVsT/data_FigD_KeyPoints_Ensemble_regen.csv",
        "Panel_d_GibbsCurveVsT/validation_FigD_GibbsCurve_Ensemble_regen.json",
        "Panel_d_GibbsCurveVsT/preview_FigD_GibbsCurve_Ensemble_regen.pdf",
        "Panel_d_GibbsCurveVsT/preview_FigD_GibbsCurve_Ensemble_regen.png",
        "Panel_e_ChemicalPotentialCascade/outputs/delta_G_rxn_v3_summary.csv",
        "Panel_e_ChemicalPotentialCascade/scripts/validate_liquid_thermodynamics_v3.py",
        "Panel_e_ChemicalPotentialCascade/outputs/validation_liquid_thermodynamics_v3.json",
        "Panel_e_ChemicalPotentialCascade/outputs/cef_validation_summary.csv",
        "Panel_e_ChemicalPotentialCascade/outputs/beta_diffusion_potentials_v3_0K.csv",
        "Panel_f_Wetting/data_FigF_Wetting_Ranked.csv",
        "Panel_f_Wetting/data_FigF_Wetting_Ranked_regen.csv",
        "SI_Figures/SI_Fig02_SizeMismatch/data_FigE_Resistance_Ranked_regen.csv",
        "SI_Figures/SI_Fig03_TripleConsensus/data_FigF_CHGNet_References.csv",
        "SI_Figures/SI_Fig03_TripleConsensus/data_FigF_TripleConsensus_Data_uma_regen.csv",
        "SI_Figures/SI_Fig03_TripleConsensus/script_FigF_Rerun_UMA.py",
        "SI_Figures/SI_Fig03_TripleConsensus/validation_FigF_UMA_Rerun.json",
        "SI_Figures/SI_Fig03_TripleConsensus/inputs/structures/binary_manifest.csv",
        "SI_Figures/SI_Fig03_TripleConsensus/inputs/structures/element_manifest.csv",
        "SI_Figures/SI_Fig03_TripleConsensus/inputs/structures/CHECKSUMS.sha256",
        "experimental_extensions/README.md",
        "experimental_extensions/FORMULAS_AND_VALIDATION.md",
        "experimental_extensions/system_manifest.json",
        "experimental_extensions/system_manifest.schema.json",
        "experimental_extensions/module_contracts.json",
        "experimental_extensions/contracts.py",
        "experimental_extensions/compare_energy_backends.py",
        "experimental_extensions/pt3_gainsnzn_regression.json",
        "experimental_extensions/run_manifold.py",
        "experimental_extensions/run_validation.py",
        "experimental_extensions/reporting.py",
        "experimental_extensions/visualization.py",
        "experimental_extensions/tests/test_regression.py",
        "experimental_extensions/tests/test_manifest_workflow.py",
        "experimental_extensions/outputs/pt3_gainsnzn_regression/summary.json",
        "experimental_extensions/outputs/system_validation/validation_results.json",
        "experimental_extensions/outputs/system_validation/validation_evidence.png",
        "experimental_extensions/outputs/system_validation/validation_evidence.pdf",
        "experimental_extensions/outputs/system_validation/evidence_report.html",
        "experimental_extensions/outputs/system_validation/evidence_report.pdf",
        "experimental_extensions/outputs/system_validation/uma_mp_dft_binary_rank/comparison_metrics.csv",
        "experimental_extensions/outputs/system_validation/uma_mp_dft_binary_rank/ranking_reversals.csv",
        "experimental_extensions/outputs/system_validation/uma_mp_dft_binary_rank/top_k_members.csv",
        # Figure 3
        "Figure3/README.md",
        "Figure3/Figure3_Code_Map.csv",
        "Figure3/fig3_paths.py",
        "Figure3/analysis/analyze_h_l8_projected_rows.py",
        "Figure3/analysis/audit_layer_coordinates.py",
        "Figure3/analysis/complete_lattice_all_layers.py",
        "Figure3/analysis/fit_per_layer_circles.py",
        "Figure3/analysis/particle_geometry.py",
        "Figure3/analysis/pipeline_atom_depth.py",
        "Figure3/build/render_compact_h.py",
        "Figure3/build/render_white_layered_3d_readable_v2.py",
        "Figure3/build/build_layer9_de_panels.py",
        "Figure3/build/build_fig1_eh_radial_model.py",
        "Figure3/build/build_c1_c2_true_numeric_axes_v8.py",
        "Figure3/build/build_c3_true_numeric_axes_v8.py",
        "Figure3/source_data/h_compact_source_data.csv",
        "Figure3/source_data/h_compact_row_pairs.csv",
        "Figure3/source_data/h_row_summary.csv",
        "Figure3/source_data/j_source_data.csv",
        "Figure3/source_data/full_audit_source_80_cells.csv",
        "Figure3/source_data/heatmap_matrix_16x5.csv",
        "Figure3/source_data/layer09_periodic_contrast_input.csv",
        "Figure3/source_data/3d_A_displayed_detected_peaks.csv",
        "Figure3/source_data/3d_B_sampled_intensity_points.csv",
        "Figure3/data/atoms_db.npz",
        "Figure3/data/columns_db.npz",
        "Figure3/data/complete_lattice_coordinates_all_layers.csv",
        "Figure3/data/complete_lattice_model.npz",
        "Figure3/data/per_layer_circles.csv",
        "Figure3/data/gray8_plain_layer_08_of_16_scale2nm.tif",
        "Figure3/data/gray8_plain_layer_09_of_16_scale2nm.tif",
        "Figure3/data/registered_stack.npz",
        "Figure3/interactive_3d/A_white_layer_resolved_readable_interactive.html",
    ]

    for rel in required:
        require_file(root, rel)

    demo_source = (root / "demo/run_reviewer_demo.py").read_text(encoding="utf-8")
    for required_mode in ['"quick"', '"full-uma"']:
        if required_mode not in demo_source:
            raise AssertionError(f"Reviewer demo is missing mode: {required_mode}")
    for unsafe_claim in [
        '"new_host_or_prototype_transferability": "demonstrated"',
        '"nonmetal_compound_transferability": "demonstrated"',
    ]:
        if unsafe_claim in demo_source:
            raise AssertionError(f"Reviewer demo contains an unsafe claim: {unsafe_claim}")

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

    fig_b_sensitivity = read_csv(
        root / "Panel_b_MultiComponentScatter/data_FigB_composition_sensitivity_regen.csv"
    )
    for c_host in (0.25, 0.75):
        rows = [row for row in fig_b_sensitivity if float(row["c_host"]) == c_host]
        pt = next(row for row in rows if row["Host"] == "Pt")
        if int(pt["dH_rank_within_composition"]) != 1:
            raise AssertionError(f"Fig. 1b sensitivity: Pt is not rank #1 at c_host={c_host}")

    fig_c = read_csv(root / "Panel_c_OrderedVsDisordered/data_FigC_Summary.csv")
    metrics = {row["Metric"]: row["Value"] for row in fig_c}
    assert_close("Fig. 1c ordering gap", float(metrics["Ordering_Gap_kJ_mol"]), 16.0436, 0.001)
    if int(float(metrics["Disordered_N"])) != 30:
        raise AssertionError("Fig. 1c disordered ensemble size is not 30")

    fig_c_raw = read_csv(
        root / "Panel_c_OrderedVsDisordered/data_FigC_Raw_UMA_Energies.csv"
    )
    if len(fig_c_raw) != 31:
        raise AssertionError("Fig. 1c raw UMA table must contain 31 structures")
    ordered_raw = [row for row in fig_c_raw if row["Type"] == "Ordered_L12"]
    disordered_raw = [
        row for row in fig_c_raw if row["Type"] == "Disordered_Random"
    ]
    if len(ordered_raw) != 1 or len(disordered_raw) != 30:
        raise AssertionError("Fig. 1c raw UMA table has an unexpected type split")
    observed_seeds = sorted(int(row["Seed"]) for row in disordered_raw)
    if observed_seeds != list(range(100, 130)):
        raise AssertionError("Fig. 1c raw UMA seeds are not exactly 100-129")
    ordered_energy = float(ordered_raw[0]["Energy_eV_atom"])
    disordered_mean_energy = sum(
        float(row["Energy_eV_atom"]) for row in disordered_raw
    ) / len(disordered_raw)
    raw_gap = (disordered_mean_energy - ordered_energy) * 96.485
    assert_close("Fig. 1c raw UMA ordering gap", raw_gap, 16.0436, 0.001)
    panel_c_report = json.loads(
        (root / "Panel_c_OrderedVsDisordered/validation_FigC_UMA_Rerun.json").read_text(
            encoding="utf-8"
        )
    )
    if panel_c_report.get("status") != "passed":
        raise AssertionError("Fig. 1c UMA rerun report did not pass")
    if float(panel_c_report["max_abs_energy_difference_eV_atom"]) > float(
        panel_c_report["tolerance_eV_atom"]
    ):
        raise AssertionError("Fig. 1c UMA rerun exceeds its reported tolerance")

    ordered_ensemble = read_csv(
        root
        / "Panel_c_OrderedVsDisordered/data_FigC_OrderedEnsemble_Raw_UMA_Energies_regen.csv"
    )
    if len(ordered_ensemble) != 30:
        raise AssertionError("Fig. 1c ordered ensemble must contain 30 structures")
    if len({row["Canonical_Class"] for row in ordered_ensemble}) != 30:
        raise AssertionError("Fig. 1c ordered ensemble symmetry classes are not unique")
    if len({row["Assignment_SHA256"] for row in ordered_ensemble}) != 30:
        raise AssertionError("Fig. 1c ordered ensemble occupancy fingerprints are not unique")
    historical_anchor = next(
        row for row in ordered_ensemble if row["Is_Historical_Anchor"] == "True"
    )
    assert_close(
        "Fig. 1c ordered ensemble historical anchor",
        float(historical_anchor["Energy_eV_atom"]),
        -4.945261,
        1.0e-6,
    )
    ensemble_report = json.loads(
        (
            root
            / "Panel_c_OrderedVsDisordered/validation_FigC_OrderedVsDisordered_Ensembles_regen.json"
        ).read_text(encoding="utf-8")
    )
    if ensemble_report.get("status") != "passed":
        raise AssertionError("Fig. 1c 30+30 ensemble report did not pass")
    if ensemble_report["ordered"]["n"] != 30 or ensemble_report["disordered"]["n"] != 30:
        raise AssertionError("Fig. 1c ensemble report does not contain 30 structures per group")
    assert_close(
        "Fig. 1c 30+30 mean ordering gap",
        float(ensemble_report["mean_gap_disordered_minus_ordered_kJ_mol_atom"]),
        15.90498,
        1.0e-4,
    )
    if float(ensemble_report["nonoverlap_margin_kJ_mol_atom"]) <= 0:
        raise AssertionError("Fig. 1c ordered and disordered ensemble ranges overlap")

    all68_report = json.loads(
        (
            root
            / "Panel_c_OrderedVsDisordered/validation_FigC_OrderedAll68_Sensitivity.json"
        ).read_text(encoding="utf-8")
    )
    if all68_report.get("status") != "passed":
        raise AssertionError("Fig. 1c all-68-class sensitivity report did not pass")
    if all68_report["all_68_equal_class_weight"]["n"] != 68:
        raise AssertionError("Fig. 1c all-class sensitivity does not contain 68 classes")
    if all68_report["all_68_degeneracy_weighted"]["raw_assignment_count"] != 2520:
        raise AssertionError("Fig. 1c class degeneracies do not recover 2520 assignments")
    if abs(float(all68_report["weighted_all68_minus_sample30_mean_kJ_mol_atom"])) > 0.01:
        raise AssertionError("Fig. 1c 30-class sample does not represent the all-class mean")

    fig_d_report = json.loads(
        (
            root
            / "Panel_d_GibbsCurveVsT/validation_FigD_GibbsCurve_Ensemble_regen.json"
        ).read_text(encoding="utf-8")
    )
    if fig_d_report.get("status") != "passed":
        raise AssertionError("Fig. 1d 30+30 ensemble report did not pass")
    if fig_d_report["ordered"]["n"] != 30 or fig_d_report["disordered"]["n"] != 30:
        raise AssertionError("Fig. 1d must use 30 ordered and 30 disordered structures")
    assert_close(
        "Fig. 1d ensemble 0 K gap",
        float(fig_d_report["gap_0K_kJ_mol_atom"]),
        15.90498,
        1.0e-4,
    )
    fig_d_key = read_csv(
        root / "Panel_d_GibbsCurveVsT/data_FigD_KeyPoints_Ensemble_regen.csv"
    )
    expected_gaps = {500.0: 13.567223, 1000.0: 11.229466, 1500.0: 8.891708}
    for temperature, expected in expected_gaps.items():
        row = next(item for item in fig_d_key if float(item["T_K"]) == temperature)
        assert_close(
            f"Fig. 1d gap at {temperature:.0f} K",
            float(row["Gap_disordered_minus_ordered_Bmix"]),
            expected,
            1.0e-5,
        )

    fig_e = read_csv(
        root / "Panel_e_ChemicalPotentialCascade/outputs/delta_G_rxn_v3_summary.csv"
    )
    assert_close(
        "Fig. 1e delta_G_rxn per f.u.",
        float(fig_e[0]["delta_G_rxn_per_fu_kJmol"]),
        -150.5481,
        0.001,
    )
    liquid_validation = json.loads(
        (
            root
            / "Panel_e_ChemicalPotentialCascade/outputs/validation_liquid_thermodynamics_v3.json"
        ).read_text(encoding="utf-8")
    )
    if liquid_validation.get("status") != "passed":
        raise AssertionError("Liquid thermodynamics validation report did not pass")
    if abs(float(liquid_validation["gibbs_duhem_euler_residual_kJmol"])) > 1.0e-10:
        raise AssertionError("Liquid partial molar terms violate the Euler identity")
    if abs(float(liquid_validation["delta_G_identity_residual_kJmol"])) > 1.0e-3:
        raise AssertionError("Panel e reaction-energy identity residual is too large")

    cef_validation = {
        row["metric"]: float(row["value"])
        for row in read_csv(
            root / "Panel_e_ChemicalPotentialCascade/outputs/cef_validation_summary.csv"
        )
    }
    assert_close("CEF training R2", cef_validation["training_R2"], 0.9992773421, 1e-9)
    assert_close(
        "CEF non-endmember LOOCV RMSE",
        cef_validation["nonendmember_LOOCV_RMSE"],
        0.1422504672,
        1e-9,
    )

    diffusion = read_csv(
        root / "Panel_e_ChemicalPotentialCascade/outputs/beta_diffusion_potentials_v3_0K.csv"
    )
    expected_diffusion = {"Ga": 0.0, "In": 91.7376, "Sn": 76.7948, "Zn": 48.9312}
    for element, expected in expected_diffusion.items():
        row = next(item for item in diffusion if item["element"] == element)
        assert_close(
            f"beta diffusion potential {element}-Ga",
            float(row["diffusion_potential_kJmol_element"]),
            expected,
            0.001,
        )

    fig_f = read_csv(root / "Panel_f_Wetting/data_FigF_Wetting_Ranked.csv")
    wetting_by_host = {row["Host"]: row for row in fig_f}
    pt_f = wetting_by_host["Pt"]
    assert_close("Fig. 1f Pt gamma_SL", float(pt_f["Gamma_SL"]), -0.44834, 0.0001)
    expected_interpretations = {
        "Pt": "favourable",
        "Co": "indeterminate",
        "Cr": "unfavourable",
    }
    for host, expected in expected_interpretations.items():
        observed = classify_gamma_sl(float(wetting_by_host[host]["Gamma_SL"]))
        if observed != expected:
            raise AssertionError(
                f"Fig. 1f {host}: interpretation {observed!r}, expected {expected!r}"
            )
    wetting_regen = read_csv(root / "Panel_f_Wetting/data_FigF_Wetting_Ranked_regen.csv")
    regen_by_host = {row["Host"]: row for row in wetting_regen}
    for host, expected in expected_interpretations.items():
        if regen_by_host[host]["Interpretation"] != expected:
            raise AssertionError(f"Fig. 1f regen {host} interpretation mismatch")

    si_fig2 = read_csv(
        root / "SI_Figures/SI_Fig02_SizeMismatch/data_FigE_Resistance_Ranked_regen.csv"
    )
    si_fig2_by_host = {row["Host"]: row for row in si_fig2}
    assert_close(
        "SI Fig. 2 corrected Pt enthalpy",
        float(si_fig2_by_host["Pt"]["Enthalpy_Drive"]),
        -33.2755,
        0.001,
    )

    si_fig3_root = root / "SI_Figures/SI_Fig03_TripleConsensus/inputs/structures"
    binary_manifest = read_csv(si_fig3_root / "binary_manifest.csv")
    element_manifest = read_csv(si_fig3_root / "element_manifest.csv")
    if len(binary_manifest) != 15 or len(element_manifest) != 16:
        raise AssertionError("SI Fig. 3 structure manifests must contain 15 binaries and 16 elements")
    binary_hosts = {row["Host"] for row in binary_manifest}
    element_names = {row["Element"] for row in element_manifest}
    if element_names != binary_hosts | {"Ga"}:
        raise AssertionError("SI Fig. 3 element references do not cover all binary hosts plus Ga")
    for row in binary_manifest:
        require_file(
            si_fig3_root,
            row["CIF_File"].replace("\\", "/"),
        )
    for row in element_manifest:
        require_file(
            si_fig3_root,
            row["CIF_File"].replace("\\", "/"),
        )
    checksum_lines = (si_fig3_root / "CHECKSUMS.sha256").read_text(
        encoding="utf-8"
    ).splitlines()
    if len(checksum_lines) != 33:
        raise AssertionError("SI Fig. 3 checksum manifest must cover all 33 source files")
    for line in checksum_lines:
        expected_hash, relative_path = line.split(maxsplit=1)
        source_path = require_file(si_fig3_root, relative_path)
        observed_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
        if observed_hash != expected_hash:
            raise AssertionError(f"SI Fig. 3 checksum mismatch: {relative_path}")
    si_fig3_report = json.loads(
        (
            root
            / "SI_Figures/SI_Fig03_TripleConsensus/validation_FigF_UMA_Rerun.json"
        ).read_text(encoding="utf-8")
    )
    if si_fig3_report.get("status") != "passed":
        raise AssertionError("SI Fig. 3 UMA rerun report did not pass")

    extension_summary = json.loads(
        (
            root
            / "experimental_extensions/outputs/pt3_gainsnzn_regression/summary.json"
        ).read_text(encoding="utf-8")
    )
    if extension_summary["n_compositions"] != 165:
        raise AssertionError("Experimental extension did not validate 165 compositions")
    assert_close(
        "Experimental extension training R2",
        float(extension_summary["training_R2"]),
        0.9992773421,
        1e-9,
    )
    if extension_summary.get("transferability_claim") != "none":
        raise AssertionError("Experimental extension must not claim scientific transferability")

    bounded = json.loads(
        (
            root
            / "experimental_extensions/outputs/system_validation/validation_results.json"
        ).read_text(encoding="utf-8")
    )
    if bounded.get("overall_status") != "passed":
        raise AssertionError("Bounded P1 validation did not pass")
    if bounded.get("scientific_transferability_claim") != "none":
        raise AssertionError("P1 report must not claim scientific transferability")
    expected_boundaries = {
        "new_host_or_prototype_transferability": "not_evaluated",
        "nonmetal_compound_transferability": "not_evaluated",
        "synthesizability": "not_predicted",
    }
    for key, expected in expected_boundaries.items():
        observed = bounded["claim_boundary"].get(key)
        if observed != expected:
            raise AssertionError(f"Unsafe P1 claim boundary {key}: {observed!r}")

    manifold_result = bounded["modules"]["fixed_l12_manifold"]["result"]
    assert_close(
        "P1 group-holdout RMSE minimum",
        float(manifold_result["group_holdout"]["RMSE_min"]),
        0.0718060334,
        1e-9,
    )
    assert_close(
        "P1 group-holdout RMSE maximum",
        float(manifold_result["group_holdout"]["RMSE_max"]),
        0.2319621909,
        1e-9,
    )
    pairwise = bounded["modules"]["uma_mp_dft_binary_rank"]["result"][
        "pairwise_results"
    ]
    pair_map = {
        (row["subset_id"], row["pair_id"]): row for row in pairwise
    }
    assert_close(
        "P1 MP-DFT vs UMA all-host Spearman",
        float(pair_map[("all_hosts", "mp_dft__uma_s_1p1")]["spearman_rho"]),
        0.8678571429,
        1e-9,
    )
    assert_close(
        "P1 MP-DFT vs CHGNet all-host Spearman",
        float(pair_map[("all_hosts", "mp_dft__chgnet")]["spearman_rho"]),
        0.95,
        1e-9,
    )
    if pair_map[("all_hosts", "mp_dft__uma_s_1p1")]["ranking_reversals"] != 15:
        raise AssertionError("Unexpected MP-DFT/UMA ranking-reversal count")
    report_html = (
        root / "experimental_extensions/outputs/system_validation/evidence_report.html"
    ).read_text(encoding="utf-8")
    for required_boundary in ["not_evaluated", "does not predict synthesizability"]:
        if required_boundary not in report_html:
            raise AssertionError(f"P1 HTML report is missing boundary: {required_boundary}")
    report_pdf = (
        root / "experimental_extensions/outputs/system_validation/evidence_report.pdf"
    ).read_bytes()
    if not report_pdf.startswith(b"%PDF"):
        raise AssertionError("P1 PDF report is not a valid PDF file")

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
    if 'version: "1.4.0-dev"' not in citation:
        raise AssertionError("CITATION.cff does not identify the current development version")

    ref_rows = read_csv(root / "shared/UMA_Element_Reference_Energies.csv")
    ref_elems = {row["Element"] for row in ref_rows}
    if ref_elems != {"Pt", "Ga", "In", "Sn", "Zn"}:
        raise AssertionError(f"Unexpected element-reference set: {sorted(ref_elems)}")
    manifest = read_csv(root / "shared/element_reference_structures/element_reference_manifest.csv")
    for row in manifest:
        require_file(root, row["Structure_File"])

    # ---- Figure 3 -----------------------------------------------------------

    # h: row-resolved projected-column modulation.
    h_pairs = read_csv(root / "Figure3/source_data/h_compact_row_pairs.csv")
    if len(h_pairs) != 8:
        raise AssertionError("Fig. 3h row-pair table must contain 8 adjacent pairs")
    contrasts = [float(row["signed_fractional_contrast"]) for row in h_pairs]
    assert_close("Fig. 3h median relative contrast",
                 statistics.median(contrasts), 0.147395, 1e-5)
    h_compact = read_csv(root / "Figure3/source_data/h_compact_source_data.csv")
    if len(h_compact) != 16:
        raise AssertionError("Fig. 3h compact profile must contain 16 complete rows")

    # i2: Layer 9 periodic-contrast score matrix is 14 x 14.
    i2_lines = (
        root / "Figure3/source_data/layer09_periodic_contrast_input.csv"
    ).read_text(encoding="utf-8-sig").splitlines()
    if len(i2_lines) != 15:
        raise AssertionError("Fig. 3i2 periodic-contrast matrix must have 15 lines (header + 14 rows)")
    for line in i2_lines[1:]:
        if len([value for value in line.split(",") if value.strip()]) != 15:
            raise AssertionError("Fig. 3i2 matrix row does not contain 14 score values")

    # i3 / j: lattice-support cells and the Layer 9 centre-edge difference.
    def clean_bom(rows):
        return [{key.lstrip("﻿"): value for key, value in row.items()} for row in rows]

    i3_cells = clean_bom(read_csv(root / "Figure3/source_data/full_audit_source_80_cells.csv"))
    if len(i3_cells) != 80:
        raise AssertionError("Fig. 3i3 audit table must contain 80 layer-radial cells")
    l9_center_cand = sum(
        int(row["candidate_locked_sites"]) for row in i3_cells
        if int(row["layer"]) == 9 and int(row["radial_bin"]) in (0, 1)
    )
    l9_center_match = sum(
        int(row["matched_measured_sites"]) for row in i3_cells
        if int(row["layer"]) == 9 and int(row["radial_bin"]) in (0, 1)
    )
    l9_edge_cand = sum(
        int(row["candidate_locked_sites"]) for row in i3_cells
        if int(row["layer"]) == 9 and int(row["radial_bin"]) == 4
    )
    l9_edge_match = sum(
        int(row["matched_measured_sites"]) for row in i3_cells
        if int(row["layer"]) == 9 and int(row["radial_bin"]) == 4
    )
    assert_close("Fig. 3j Layer 9 centre fraction",
                 l9_center_match / l9_center_cand, 127 / 131, 1e-9)
    assert_close("Fig. 3j Layer 9 edge fraction",
                 l9_edge_match / l9_edge_cand, 201 / 289, 1e-9)
    j_rows = clean_bom(read_csv(root / "Figure3/source_data/j_source_data.csv"))
    j9 = next(row for row in j_rows if int(row["layer"]) == 9)
    assert_close("Fig. 3j Layer 9 delta",
                 float(j9["delta_center_minus_outer"]), 0.273964, 1e-5)

    print(f"Release root: {root}")
    print("All release verification checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
