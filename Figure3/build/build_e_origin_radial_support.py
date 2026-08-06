"""Build two editable Origin panels for particle-wide radial lattice support.

E1 is a 16 x 5 layer/radius heatmap. E2 is a three-curve weighted radial
summary. The source values are copied into each OPJU and are never smoothed,
fitted, normalized again, or interpreted as phase fractions.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import originpro as op


PROJECT = Path(r"C:\Users\13102\Documents\Project Review")
SOURCE_DIR = PROJECT / "outputs" / "FIG1_EH_RADIAL_MODEL_20260803_v3"
SOURCE = SOURCE_DIR / "panel_e_radial_depth_source_data.csv"
OUTPUT = PROJECT / "outputs" / "FIG1_E_ORIGIN_RADIAL_SUPPORT_20260803_v5"
COMPAT_SRC = PROJECT / "editaplot_runtime_origin2026b_compat_v1" / "src"
sys.path.insert(0, str(COMPAT_SRC))

from origin_sciplot.origin_backend.export_utils import export_graph  # noqa: E402
from origin_sciplot.origin_backend.base_style_contract import pt_to_origin_width_units  # noqa: E402


FONT_CODE = 73  # Arial in the locally verified Origin 2026b environment.
COLORS = {"Center": "#3C5488", "Middle": "#00A087", "Outer edge": "#D55E5E"}
MARKERS = {"Center": 2, "Middle": 3, "Outer edge": 4}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric_hash(frame: pd.DataFrame) -> str:
    values = np.ascontiguousarray(frame.to_numpy(dtype=np.float64))
    return hashlib.sha256(values.tobytes()).hexdigest()


def tick_text(values) -> str:
    return " ".join(f"{value:g}" for value in values)


def style_axis(layer, axis: str, decimals: int, label_size: float = 15.0) -> None:
    layer.set_int(f"{axis}.label.type", 1)
    layer.set_int(f"{axis}.label.numFormat", 1)
    layer.set_int(f"{axis}.label.decimalPlaces", decimals)
    layer.set_int(f"{axis}.ticks", 10)  # major and minor ticks point outward
    layer.set_int(f"{axis}.minorTicks", 0)
    layer.set_float(f"{axis}.ticklength", 5.5)
    layer.set_float(f"{axis}.mticklength", 3.0)
    layer.set_float(f"{axis}.tickthickness", 1.7)
    layer.set_float(f"{axis}.mtickthickness", 1.1)
    layer.set_int(f"{axis}.showLabels", 1)
    layer.set_int(f"{axis}.label.show", 1)
    layer.set_int(f"{axis}.label.font", FONT_CODE)
    layer.set_float(f"{axis}.label.pt", label_size)


def exact_ticks(layer, x_major, x_minor, y_major, y_minor) -> None:
    command = (
        f'layer.x.ticksbydata$="{tick_text(x_major)}";'
        f'layer.x.mticksbydata$="{tick_text(x_minor)}";'
        f'layer.y.ticksbydata$="{tick_text(y_major)}";'
        f'layer.y.mticksbydata$="{tick_text(y_minor)}";doc -uw;'
    )
    if not layer.obj.LT_execute(command):
        raise RuntimeError("Origin rejected exact major/minor tick positions")


def set_titles(layer, x_title: str, y_title: str, size: float = 18.0) -> None:
    command = (
        f'xb.text$="{x_title}";yl.text$="{y_title}";'
        f"xb.show=1;yl.show=1;xb.font={FONT_CODE};yl.font={FONT_CODE};"
        f"xb.fsize={size:g};yl.fsize={size:g};doc -uw;"
    )
    if not layer.obj.LT_execute(command):
        raise RuntimeError("Origin rejected axis-title settings")


def axis_readback(layer, axis: str) -> dict[str, object]:
    major_var = f"__{axis}MajorE"
    minor_var = f"__{axis}MinorE"
    if not layer.obj.LT_execute(
        f"string {major_var}$=layer.{axis}.ticksbydata$;"
        f"string {minor_var}$=layer.{axis}.mticksbydata$;"
    ):
        raise RuntimeError(f"Could not read back {axis} ticks")
    return {
        "from": float(layer.get_float(f"{axis}.from")),
        "to": float(layer.get_float(f"{axis}.to")),
        "increment": float(layer.get_float(f"{axis}.inc")),
        "tick_direction_code": int(layer.get_int(f"{axis}.ticks")),
        "major_length": float(layer.get_float(f"{axis}.ticklength")),
        "minor_length": float(layer.get_float(f"{axis}.mticklength")),
        "label_font_code": int(layer.get_int(f"{axis}.label.font")),
        "label_size_pt": float(layer.get_float(f"{axis}.label.pt")),
        "major_ticks_by_data": op.get_lt_str(f"{major_var}$"),
        "minor_ticks_by_data": op.get_lt_str(f"{minor_var}$"),
    }


def add_page_label(layer, name: str, text: str, left: float, top: float, size: float = 12.0):
    label = layer.add_label(text)
    if label is None:
        raise RuntimeError(f"Could not create label {name}")
    label.name = name
    label.set_int("attach", 0)
    label.set_float("left", left)
    label.set_float("top", top)
    label.set_int("font", FONT_CODE)
    label.set_float("fsize", size)
    label.set_int("frame", 0)
    label.set_int("showframe", 0)
    return label


def export_project(graph, out: Path, report: dict[str, object]) -> dict[str, object]:
    graph.activate()
    op.lt_exec("doc -uw;")
    opju = out / "result.opju"
    if not op.save(str(opju)):
        raise RuntimeError("Origin failed to save the editable OPJU")
    export_graph(op, graph, out / "result.png", out / "result.pdf", out / "result.tif")
    report["artifacts"] = {
        name: {"path": str(out / name), "sha256": sha256_file(out / name)}
        for name in ("result.opju", "result.png", "result.pdf", "result.tif")
    }
    (out / "origin_readback.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return report


def prepare_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source = pd.read_csv(SOURCE)
    if len(source) != 80 or set(source["layer"]) != set(range(1, 17)):
        raise RuntimeError("Expected exactly 80 radial cells across Layers 1-16")
    expected = source["matched_measured_sites"] / source["candidate_locked_sites"]
    if not np.allclose(expected, source["lattice_support_fraction"], rtol=0, atol=1e-14):
        raise RuntimeError("Source support fractions fail matched/candidate recomputation")

    matrix = (
        source.pivot(index="layer", columns="radial_bin", values="lattice_support_fraction")
        .reindex(index=range(1, 17), columns=range(5))
    )
    matrix.columns = ["r_0.0_0.2", "r_0.2_0.4", "r_0.4_0.6", "r_0.6_0.8", "r_0.8_1.0"]
    matrix.insert(0, "Layer index", matrix.index.astype(int))
    matrix = matrix.reset_index(drop=True)

    group_map = {0: "Center", 1: "Center", 2: "Middle", 3: "Middle", 4: "Outer edge"}
    grouped = source.assign(radial_group=source["radial_bin"].map(group_map))
    weighted = (
        grouped.groupby(["layer", "radial_group"], sort=False)[
            ["matched_measured_sites", "candidate_locked_sites"]
        ]
        .sum()
        .reset_index()
    )
    weighted["weighted_support"] = (
        weighted["matched_measured_sites"] / weighted["candidate_locked_sites"]
    )
    trends = weighted.pivot(index="layer", columns="radial_group", values="weighted_support")
    trends = trends.reindex(index=range(1, 17), columns=["Center", "Middle", "Outer edge"])
    trends.insert(0, "Layer index", trends.index.astype(int))
    trends = trends.reset_index(drop=True)
    return source, matrix, trends


def build_heatmap(source: pd.DataFrame, matrix_frame: pd.DataFrame) -> dict[str, object]:
    out = OUTPUT / "E1_particle_wide_support_heatmap"
    out.mkdir(parents=True)
    source.to_csv(out / "full_audit_source_80_cells.csv", index=False, encoding="utf-8-sig")
    matrix_frame.to_csv(out / "heatmap_matrix_16x5.csv", index=False, encoding="utf-8-sig")

    values = matrix_frame.iloc[:, 1:].to_numpy(float)
    expected_hash = numeric_hash(matrix_frame)
    op.set_show(False)
    try:
        op.new(asksave=False)
        audit_sheet = op.new_sheet("w", "E1 Full Audit Source")
        table_sheet = op.new_sheet("w", "E1 Heatmap Table")
        matrix = op.new_sheet("m", "E1 Support Matrix")
        if audit_sheet is None or table_sheet is None or matrix is None:
            raise RuntimeError("Could not create E1 source objects")
        audit_sheet.from_df(source)
        table_sheet.from_df(matrix_frame)
        matrix.from_np(values)
        matrix.xymap = (0.1, 0.9, 1.0, 16.0)

        graph = op.new_graph("E1 particle-wide radial support", template="heatmap")
        if graph is None:
            raise RuntimeError("Could not create E1 graph")
        graph.set_int("background", op.ocolor("#FFFFFF"))
        layer = graph[0]
        plot = layer.add_plot(matrix, colz=0)
        if plot is None:
            raise RuntimeError("Could not create E1 matrix plot")
        plot.colormap = "Viridis.PAL"
        layer.rescale("z")
        if not layer.obj.LT_execute(
            "layer.cmap.zmin=0;layer.cmap.zmax=1;layer.cmap.flippal=1;"
            "layer.cmap.updateScale();doc -uw;"
        ):
            raise RuntimeError("Could not apply E1 0-1 sequential color scale")

        layer.axis("x").set_limits(0.0, 1.0, 0.2)
        layer.axis("y").set_limits(0.5, 16.5, 1.0)
        style_axis(layer, "x", 1)
        style_axis(layer, "y", 0)
        x_major = np.arange(0.0, 1.01, 0.2)
        x_minor = np.arange(0.1, 1.0, 0.2)
        y_major = range(1, 17)
        y_minor = [value + 0.5 for value in range(1, 16)]
        exact_ticks(layer, x_major, x_minor, y_major, y_minor)
        set_titles(layer, "Normalized radius, r/R", "Reconstructed layer index")

        layer.set_int("unit", 1)
        layer.set_float("left", 13.0)
        layer.set_float("top", 7.0)
        layer.set_float("width", 68.0)
        layer.set_float("height", 78.0)
        graph.activate()
        graph.obj.SetWidth(8.0)
        graph.obj.SetHeight(7.2)

        legend = layer.label("Legend") or layer.label("legend")
        if legend is not None:
            legend.remove()

        color_major = np.arange(0.0, 1.01, 0.2)
        color_minor = np.arange(0.1, 1.0, 0.2)
        graph.activate()
        page_width = float(op.lt_float("page.width"))
        page_height = float(op.lt_float("page.height"))
        if not layer.obj.LT_execute(
            "spectrum1.title=0;spectrum1.show=1;spectrum1.attach=0;"
            "spectrum1.labels.autodisp=0;spectrum1.labels.fsize=12;"
            f"spectrum1.labels.font={FONT_CODE};spectrum1.labels.bold=0;"
            "spectrum1.labels.numFormat=1;spectrum1.labels.decimalPlaces=1;"
            "spectrum1.top=page.height*0.11;spectrum1.left=page.width*0.855;"
            "spectrum1.width=page.width*0.11;spectrum1.height=page.height*0.68;"
            "spectrum1.barthick=100;doc -uw;"
        ):
            raise RuntimeError("Could not style the native E1 colorbar")
        add_page_label(
            layer, "E1ColorTitle", "Measured lattice support fraction",
            0.68 * page_width, 0.025 * page_height, 11.5,
        )

        readback = table_sheet.to_df()[matrix_frame.columns.tolist()]
        if numeric_hash(readback) != expected_hash:
            raise RuntimeError("E1 Origin worksheet differs from input matrix")
        report = {
            "panel": "E1",
            "scientific_metric": "matched measured lattice sites / locked candidate lattice sites",
            "not_interpretable_as": ["crystallinity percentage", "phase fraction", "composition"],
            "source_rows": 80,
            "locked_coordinate_baseline": 8748,
            "candidate_sites_in_radial_domain": int(source["candidate_locked_sites"].sum()),
            "outside_r_over_R_0_to_1": 8748 - int(source["candidate_locked_sites"].sum()),
            "matrix_shape": list(values.shape),
            "input_numeric_hash": expected_hash,
            "origin_numeric_hash": numeric_hash(readback),
            "source_values_modified": False,
            "smoothing_or_fit": False,
            "x_axis": axis_readback(layer, "x"),
            "y_axis": axis_readback(layer, "y"),
            "palette": plot.colormap,
            "palette_flipped": bool(round(layer.get_float("cmap.flippal"))),
            "color_semantics": "low support = dark purple/blue; high support = bright yellow",
            "color_range": [0.0, 1.0],
            "colorbar_major_ticks": list(map(float, color_major)),
            "colorbar_minor_ticks": list(map(float, color_minor)),
        }
        return export_project(graph, out, report)
    finally:
        op.exit()


def build_trends(source: pd.DataFrame, trends: pd.DataFrame) -> dict[str, object]:
    out = OUTPUT / "E2_radial_support_trends"
    out.mkdir(parents=True)
    source.to_csv(out / "full_audit_source_80_cells.csv", index=False, encoding="utf-8-sig")
    trends.to_csv(out / "weighted_radial_trends.csv", index=False, encoding="utf-8-sig")
    expected_hash = numeric_hash(trends)

    op.set_show(False)
    try:
        op.new(asksave=False)
        audit_sheet = op.new_sheet("w", "E2 Full Audit Source")
        sheet = op.new_sheet("w", "E2 Weighted Trends")
        if audit_sheet is None or sheet is None:
            raise RuntimeError("Could not create E2 worksheets")
        audit_sheet.from_df(source)
        sheet.from_df(trends)
        sheet.cols_axis("xyyy")

        graph = op.new_graph("E2 radial support trends", template="Line")
        if graph is None:
            raise RuntimeError("Could not create E2 graph")
        graph.set_int("background", op.ocolor("#FFFFFF"))
        layer = graph[0]
        plots = []
        for label in ("Center", "Middle", "Outer edge"):
            plot = layer.add_plot(sheet, label, "Layer index", type="y")
            if plot is None:
                raise RuntimeError(f"Could not create E2 series: {label}")
            plot.color = COLORS[label]
            plot.set_cmd(f"-c color({COLORS[label]})", f"-w {pt_to_origin_width_units(2.2)}", "-d 0")
            plot.symbol_kind = MARKERS[label]
            plot.symbol_interior = 1
            plot.symbol_size = 8.0
            plot.set_cmd("-kh 35")
            plots.append(plot)

        layer.axis("x").set_limits(0.5, 16.5, 1.0)
        layer.axis("y").set_limits(0.0, 1.0, 0.2)
        style_axis(layer, "x", 0)
        style_axis(layer, "y", 1)
        x_major = range(1, 17)
        x_minor = [value + 0.5 for value in range(1, 16)]
        y_major = np.arange(0.0, 1.01, 0.2)
        y_minor = np.arange(0.1, 1.0, 0.2)
        exact_ticks(layer, x_major, x_minor, y_major, y_minor)
        set_titles(layer, "Reconstructed layer index", "Measured lattice support fraction")
        layer.set_int("unit", 1)
        layer.set_float("left", 13.0)
        layer.set_float("top", 7.0)
        layer.set_float("width", 78.0)
        layer.set_float("height", 76.0)
        graph.activate()
        graph.obj.SetWidth(8.0)
        graph.obj.SetHeight(6.0)

        layer.obj.LT_execute("legend -r;doc -uw;")
        legend = layer.label("legend") or layer.label("Legend")
        if legend is not None:
            legend.text = "\\l(1) Center (r/R = 0-0.4)\n\\l(2) Middle (r/R = 0.4-0.8)\n\\l(3) Outer edge (r/R = 0.8-1.0)"
            legend.set_int("attach", 0)
            legend.set_int("font", FONT_CODE)
            legend.set_float("fsize", 11.5)
            legend.set_int("frame", 0)
            legend.set_int("showframe", 0)
            graph.activate()
            page_width = float(op.lt_float("page.width"))
            page_height = float(op.lt_float("page.height"))
            legend.set_float("left", 0.55 * page_width)
            legend.set_float("top", 0.08 * page_height)

        readback = sheet.to_df()[trends.columns.tolist()]
        if numeric_hash(readback) != expected_hash:
            raise RuntimeError("E2 Origin worksheet differs from weighted input")
        report = {
            "panel": "E2",
            "scientific_metric": "candidate-count-weighted measured lattice support fraction",
            "locked_coordinate_baseline": 8748,
            "candidate_sites_in_radial_domain": int(source["candidate_locked_sites"].sum()),
            "outside_r_over_R_0_to_1": 8748 - int(source["candidate_locked_sites"].sum()),
            "aggregation": {
                "Center": "sum matched in r/R 0-0.4 / sum candidates in r/R 0-0.4",
                "Middle": "sum matched in r/R 0.4-0.8 / sum candidates in r/R 0.4-0.8",
                "Outer edge": "sum matched in r/R 0.8-1.0 / sum candidates in r/R 0.8-1.0",
            },
            "unweighted_bin_average_used": False,
            "source_values_modified": False,
            "smoothing_or_fit": False,
            "input_numeric_hash": expected_hash,
            "origin_numeric_hash": numeric_hash(readback),
            "series_colors": COLORS,
            "series_markers": MARKERS,
            "x_axis": axis_readback(layer, "x"),
            "y_axis": axis_readback(layer, "y"),
            "legend_labels": [
                "Center (r/R = 0-0.4)",
                "Middle (r/R = 0.4-0.8)",
                "Outer edge (r/R = 0.8-1.0)",
            ],
        }
        return export_project(graph, out, report)
    finally:
        op.exit()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    source, matrix, trends = prepare_inputs()
    source_hash_before = sha256_file(SOURCE)

    contract = {
        "status": "confirmed from the user's latest request",
        "source": str(SOURCE),
        "source_sha256": source_hash_before,
        "locked_coordinate_baseline": 8748,
        "candidate_sites_in_radial_domain": int(source["candidate_locked_sites"].sum()),
        "outside_r_over_R_0_to_1": 8748 - int(source["candidate_locked_sites"].sum()),
        "E1_role": "particle-wide layer-by-radius distribution of measured lattice support",
        "E2_role": "center/middle/outer weighted depth trend",
        "visible_primary_columns": ["layer", "radial_bin", "lattice_support_fraction"],
        "support_only_columns": ["radius_inner", "radius_outer", "candidate_locked_sites", "matched_measured_sites"],
        "retain_not_render": [
            "parity0_count", "parity1_count", "auc_local_order_s_parity0",
            "alternating_row_order_score", "order_cell_valid", "low_sample_reason",
        ],
        "uncertain_columns": [],
        "forbidden_claims": ["phase fraction", "composition", "crystallinity percentage", "time evolution"],
        "display_transforms": ["pivot to 16x5 matrix", "candidate-count-weighted radial aggregation"],
        "not_performed": ["smoothing", "fitting", "normalization", "outlier removal", "phase assignment"],
    }
    (OUTPUT / "figure_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    reports = [build_heatmap(source, matrix), build_trends(source, trends)]
    if sha256_file(SOURCE) != source_hash_before:
        raise RuntimeError("The immutable source CSV changed during Origin rendering")
    summary = {
        "status": "origin_exports_complete_pending_visual_qa",
        "origin_version": "2026b",
        "source_sha256_unchanged": True,
        "locked_coordinate_baseline": 8748,
        "candidate_sites_in_radial_domain": int(source["candidate_locked_sites"].sum()),
        "outside_r_over_R_0_to_1": 8748 - int(source["candidate_locked_sites"].sum()),
        "reports": reports,
    }
    (OUTPUT / "origin_verify_report.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    readme = "# Fig. 1E Origin radial-support figures\n\n"
    readme += "- E1: 8,512 locked candidate sites inside r/R = 0-1 summarized as a 16-layer x 5-radius-bin heatmap.\n"
    readme += "- The immutable baseline contains 8,748 locked coordinates; 236 outside the fitted particle radius are retained but not binned.\n"
    readme += "- E2: candidate-count-weighted Center/Middle/Outer-edge support trends.\n"
    readme += "- Color/Y values are measured lattice support fractions, not phase fraction, composition, or crystallinity percentage.\n"
    readme += "- No smoothing, fitting, outlier removal, or coordinate modification was performed.\n"
    (OUTPUT / "README.md").write_text(readme, encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "panels": ["E1", "E2"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
