"""Build editable Origin2025b h and compact j panels for Figure 1.

The OPJU retains the full audited source tables. Panel h is a fixed-lattice-
parity ECDF. Panel j is a compact two-curve summary derived from the same
16 x 5 layer/radius audit table used by i3. No smoothing or fitting is used.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(r"C:\Users\13102\Documents\Project Review")
VENDOR = PROJECT / "tools" / "originpro_vendor"
RUNTIME_SRC = PROJECT / "editaplot_runtime_origin2026b_compat_v1" / "src"
sys.path.insert(0, str(VENDOR))
sys.path.insert(0, str(RUNTIME_SRC))

import originpro as op  # noqa: E402

from origin_sciplot.origin_backend.base_style_contract import (  # noqa: E402
    pt_to_origin_width_units,
)
from origin_sciplot.origin_backend.export_utils import export_graph  # noqa: E402


PYTHON_BUNDLE = PROJECT / "outputs" / "FIG1_HJ_SCIENCE_LOGIC_PYTHON_20260804_v1"
H_SOURCE = PYTHON_BUNDLE / "h_column_source_data.csv"
J_SOURCE = PYTHON_BUNDLE / "j_source_data.csv"
J_AUDIT_SOURCE = (
    PROJECT
    / "outputs"
    / "FIG1_E_ORIGIN_RADIAL_SUPPORT_20260803_v5"
    / "E1_particle_wide_support_heatmap"
    / "full_audit_source_80_cells.csv"
)
OUTPUT = PROJECT / "outputs" / "FIG1_HJ_ORIGIN2025B_20260804_v4"

FONT_CODE = 73
BLUE = "#356D9C"
ORANGE = "#D58A3A"
GREY = "#777777"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric_hash(frame: pd.DataFrame) -> str:
    values = np.ascontiguousarray(frame.to_numpy(dtype=np.float64))
    return hashlib.sha256(values.tobytes()).hexdigest()


def assert_numeric_equal(expected: pd.DataFrame, actual: pd.DataFrame, name: str) -> None:
    if expected.shape != actual.shape:
        raise RuntimeError(f"{name} shape mismatch: {expected.shape} != {actual.shape}")
    if not np.allclose(
        expected.to_numpy(dtype=float),
        actual.to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    ):
        raise RuntimeError(f"{name} Origin worksheet differs from input")


def tick_text(values) -> str:
    return " ".join(f"{value:g}" for value in values)


def style_axis(layer, axis: str, decimals: int, label_size: float = 9.0) -> None:
    layer.set_int(f"{axis}.label.type", 1)
    layer.set_int(f"{axis}.label.numFormat", 1)
    layer.set_int(f"{axis}.label.decimalPlaces", decimals)
    layer.set_float(f"{axis}.label.rotate", 0.0)
    layer.set_int(f"{axis}.ticks", 10)
    layer.set_int(f"{axis}.minorTicks", 0)
    layer.set_float(f"{axis}.ticklength", 4.0)
    layer.set_float(f"{axis}.mticklength", 2.2)
    layer.set_float(f"{axis}.tickthickness", 1.2)
    layer.set_float(f"{axis}.mtickthickness", 0.8)
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
        raise RuntimeError("Origin rejected exact tick positions")


def set_titles(layer, x_title: str, y_title: str, size: float = 10.5) -> None:
    command = (
        f'xb.text$="{x_title}";yl.text$="{y_title}";'
        f"xb.show=1;yl.show=1;xb.font={FONT_CODE};yl.font={FONT_CODE};"
        f"xb.fsize={size:g};yl.fsize={size:g};doc -uw;"
    )
    if not layer.obj.LT_execute(command):
        raise RuntimeError("Origin rejected axis titles")


def add_page_label(layer, name: str, text: str, left: float, top: float, size: float):
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


def axis_readback(layer, axis: str, suffix: str) -> dict[str, object]:
    major_var = f"__{axis}Major{suffix}"
    minor_var = f"__{axis}Minor{suffix}"
    if not layer.obj.LT_execute(
        f"string {major_var}$=layer.{axis}.ticksbydata$;"
        f"string {minor_var}$=layer.{axis}.mticksbydata$;"
    ):
        raise RuntimeError(f"Could not read back {axis} ticks")
    return {
        "from": float(layer.get_float(f"{axis}.from")),
        "to": float(layer.get_float(f"{axis}.to")),
        "tick_direction_code": int(layer.get_int(f"{axis}.ticks")),
        "major_length": float(layer.get_float(f"{axis}.ticklength")),
        "minor_length": float(layer.get_float(f"{axis}.mticklength")),
        "label_font_code": int(layer.get_int(f"{axis}.label.font")),
        "label_size_pt": float(layer.get_float(f"{axis}.label.pt")),
        "major_ticks_by_data": op.get_lt_str(f"{major_var}$"),
        "minor_ticks_by_data": op.get_lt_str(f"{minor_var}$"),
    }


def step_ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=float))
    y = np.arange(1, x.size + 1, dtype=float) / x.size
    return np.r_[x[0], np.repeat(x, 2)[1:]], np.r_[0.0, np.repeat(y, 2)[:-1]]


def prepare_h() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    source = pd.read_csv(H_SOURCE)
    required = {"parity_i_plus_j", "local_normalized_intensity"}
    if not required.issubset(source.columns):
        raise RuntimeError("H source columns are incomplete")
    if set(source["parity_i_plus_j"].unique()) != {0, 1}:
        raise RuntimeError("H source must contain fixed parity groups 0 and 1")
    v0 = source.loc[source["parity_i_plus_j"] == 0, "local_normalized_intensity"].to_numpy()
    v1 = source.loc[source["parity_i_plus_j"] == 1, "local_normalized_intensity"].to_numpy()
    x0, y0 = step_ecdf(v0)
    x1, y1 = step_ecdf(v1)
    length = max(len(x0), len(x1))
    plot = pd.DataFrame(index=range(length))
    plot["Parity 0 intensity"] = pd.Series(x0)
    plot["Parity 0 ECDF"] = pd.Series(y0)
    plot["Parity 1 intensity"] = pd.Series(x1)
    plot["Parity 1 ECDF"] = pd.Series(y1)
    stats = {
        "n_total": int(len(source)),
        "n_parity_0": int(len(v0)),
        "n_parity_1": int(len(v1)),
        "median_parity_0": float(np.median(v0)),
        "median_parity_1": float(np.median(v1)),
        "median_ratio": float(np.median(v0) / np.median(v1)),
        "descriptive_auc": 0.9590402230187176,
    }
    return source, plot, stats


def prepare_j() -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(J_SOURCE)
    audit = pd.read_csv(J_AUDIT_SOURCE)
    required = {
        "layer",
        "center_support_fraction",
        "outer_edge_support_fraction",
        "delta_center_minus_outer",
        "center_candidate_sites",
        "outer_edge_candidate_sites",
        "center_matched_sites",
        "outer_edge_matched_sites",
    }
    if not required.issubset(summary.columns):
        raise RuntimeError("J summary columns are incomplete")
    if len(summary) != 16 or len(audit) != 80:
        raise RuntimeError("J requires 16 summary rows and 80 audit cells")
    recomputed = summary["center_support_fraction"] - summary["outer_edge_support_fraction"]
    if not np.allclose(recomputed, summary["delta_center_minus_outer"], atol=1e-12):
        raise RuntimeError("J signed differences fail recomputation")
    return summary, audit


def create_definitions_sheet(h_stats: dict[str, float]) -> None:
    definitions = pd.DataFrame(
        {
            "Item": [
                "h role",
                "h grouping",
                "h normalization",
                "h inference",
                "j role",
                "j center",
                "j outer edge",
                "j display reduction",
                "forbidden claims",
            ],
            "Definition": [
                "Real-space validation of ordered-region reciprocal-space evidence",
                "Fixed fitted-lattice parity: (i+j) mod 2; not intensity-defined",
                "Background-subtracted integrated column intensity divided by median of 8 nearest neighbors",
                f"Descriptive only; n={int(h_stats['n_parity_0'])}/{int(h_stats['n_parity_1'])}; no p-value",
                "Compact quantitative summary derived from the same source cells as i3",
                "sum matched / sum candidates for r/R = 0-0.4",
                "sum matched / sum candidates for r/R = 0.8-1.0",
                "Candidate counts retained in worksheets, omitted from compact graph",
                "No phase fraction, crystallinity, composition, chemical identification, or time evolution",
            ],
        }
    )
    sheet = op.new_sheet("w", "HJ Definitions")
    if sheet is None:
        raise RuntimeError("Could not create definitions worksheet")
    sheet.from_df(definitions)


def build_h(source: pd.DataFrame, plot_frame: pd.DataFrame, stats: dict[str, float]):
    raw_sheet = op.new_sheet("w", "H Column Source")
    plot_sheet = op.new_sheet("w", "H ECDF Plot Data")
    if raw_sheet is None or plot_sheet is None:
        raise RuntimeError("Could not create H worksheets")
    raw_sheet.from_df(source)
    plot_sheet.from_df(plot_frame)
    plot_sheet.cols_axis("xyxy")

    graph = op.new_graph("Fig1 h sublattice intensity ECDF", template="Line")
    if graph is None:
        raise RuntimeError("Could not create H graph")
    graph.set_int("background", op.ocolor("#FFFFFF"))
    layer = graph[0]

    p0 = layer.add_plot(plot_sheet, "Parity 0 ECDF", "Parity 0 intensity", type="y")
    p1 = layer.add_plot(plot_sheet, "Parity 1 ECDF", "Parity 1 intensity", type="y")
    if p0 is None or p1 is None:
        raise RuntimeError("Could not create H ECDF curves")
    for plot, color in ((p0, BLUE), (p1, ORANGE)):
        plot.color = op.ocolor(color)
        plot.symbol_kind = 0
        plot.set_cmd(
            f"-c color({color})",
            f"-w {pt_to_origin_width_units(1.45)}",
            "-k 0",
            "-d 0",
        )

    median_frame = pd.DataFrame(
        {
            "Median x 0": [stats["median_parity_0"]],
            "Median y 0": [0.5],
            "Median x 1": [stats["median_parity_1"]],
            "Median y 1": [0.5],
        }
    )
    median_sheet = op.new_sheet("w", "H Median Markers")
    if median_sheet is None:
        raise RuntimeError("Could not create H median worksheet")
    median_sheet.from_df(median_frame)
    median_sheet.cols_axis("xyxy")
    for y_col, x_col, color in (
        ("Median y 0", "Median x 0", BLUE),
        ("Median y 1", "Median x 1", ORANGE),
    ):
        marker = layer.add_plot(median_sheet, y_col, x_col, type="s")
        if marker is None:
            raise RuntimeError("Could not create H median marker")
        marker.color = op.ocolor(color)
        marker.symbol_kind = 2
        marker.symbol_interior = 1
        marker.symbol_size = 5.8
        marker.set_cmd(f"-c color({color})", "-kh 35")

    layer.axis("x").set_limits(0.5, 1.32, 0.1)
    layer.axis("y").set_limits(0.0, 1.0, 0.25)
    style_axis(layer, "x", 1, 8.8)
    style_axis(layer, "y", 2, 8.8)
    exact_ticks(
        layer,
        np.arange(0.6, 1.31, 0.1),
        np.arange(0.55, 1.31, 0.1),
        np.arange(0, 1.01, 0.25),
        np.arange(0.125, 1.0, 0.25),
    )
    set_titles(layer, "Column intensity / local 8-neighbor median", "Cumulative fraction")
    layer.set_int("unit", 1)
    layer.set_float("left", 17.0)
    layer.set_float("top", 12.0)
    layer.set_float("width", 79.0)
    layer.set_float("height", 76.0)
    graph.activate()
    graph.obj.SetWidth(7.8)
    graph.obj.SetHeight(5.7)

    layer.obj.LT_execute("legend -r;doc -uw;")
    legend = layer.label("legend") or layer.label("Legend")
    if legend is None:
        raise RuntimeError("H legend is missing")
    legend.text = (
        f"\\l(1) Parity 0 (n={int(stats['n_parity_0'])})\n"
        f"\\l(2) Parity 1 (n={int(stats['n_parity_1'])})"
    )
    legend.set_int("attach", 0)
    legend.set_int("font", FONT_CODE)
    legend.set_float("fsize", 8.2)
    legend.set_int("frame", 0)
    legend.set_int("showframe", 0)
    graph.activate()
    page_width = float(op.lt_float("page.width"))
    page_height = float(op.lt_float("page.height"))
    legend.set_float("left", 0.19 * page_width)
    legend.set_float("top", 0.10 * page_height)
    add_page_label(layer, "panel_h", "h", 8.0, 5.5, 13.0)
    add_page_label(
        layer,
        "HStats",
        f"Median ratio = {stats['median_ratio']:.2f}\nDescriptive AUC = {stats['descriptive_auc']:.3f}",
        0.19 * page_width,
        0.24 * page_height,
        7.2,
    )

    readback = plot_sheet.to_df()[plot_frame.columns.tolist()]
    assert_numeric_equal(plot_frame, readback, "H plot data")
    return graph, {
        "panel": "h",
        "source_rows": int(len(source)),
        "statistics": stats,
        "normalization": "column intensity / median intensity of 8 nearest neighboring columns",
        "input_numeric_hash": numeric_hash(plot_frame),
        "origin_numeric_hash": numeric_hash(readback),
        "source_values_modified": False,
        "smoothing_or_fit": False,
        "x_axis": axis_readback(layer, "x", "H"),
        "y_axis": axis_readback(layer, "y", "H"),
    }


def build_j(summary: pd.DataFrame, audit: pd.DataFrame):
    audit_sheet = op.new_sheet("w", "J Full Audit 80 Cells")
    summary_sheet = op.new_sheet("w", "J Compact Summary")
    if audit_sheet is None or summary_sheet is None:
        raise RuntimeError("Could not create J worksheets")
    audit_sheet.from_df(audit)
    summary_sheet.from_df(summary)

    plot_frame = summary[
        ["layer", "center_support_fraction", "outer_edge_support_fraction"]
    ].copy()
    plot_sheet = op.new_sheet("w", "J Plot Data")
    if plot_sheet is None:
        raise RuntimeError("Could not create J plot worksheet")
    plot_sheet.from_df(plot_frame)
    plot_sheet.cols_axis("xyy")

    graph = op.new_graph("Fig1 j compact center outer support", template="Line")
    if graph is None:
        raise RuntimeError("Could not create J graph")
    graph.set_int("background", op.ocolor("#FFFFFF"))
    layer = graph[0]
    series = [
        ("center_support_fraction", BLUE, 2),
        ("outer_edge_support_fraction", ORANGE, 3),
    ]
    for column, color, marker_kind in series:
        plot = layer.add_plot(plot_sheet, column, "layer", type="y")
        if plot is None:
            raise RuntimeError(f"Could not create J series {column}")
        plot.color = op.ocolor(color)
        plot.set_cmd(f"-c color({color})", f"-w {pt_to_origin_width_units(1.45)}", "-d 0")
        plot.symbol_kind = marker_kind
        plot.symbol_interior = 1
        plot.symbol_size = 5.8
        plot.set_cmd("-kh 35")

    layer.axis("x").set_limits(0.5, 16.5, 1.0)
    layer.axis("y").set_limits(0.0, 1.0, 0.2)
    style_axis(layer, "x", 0, 8.3)
    style_axis(layer, "y", 1, 8.8)
    exact_ticks(
        layer,
        range(1, 17),
        [value + 0.5 for value in range(1, 16)],
        np.arange(0, 1.01, 0.2),
        np.arange(0.1, 1.0, 0.2),
    )
    set_titles(layer, "Reconstructed layer index", "Measured lattice support fraction")
    layer.set_int("unit", 1)
    layer.set_float("left", 15.5)
    layer.set_float("top", 17.0)
    layer.set_float("width", 81.0)
    layer.set_float("height", 72.0)
    graph.activate()
    graph.obj.SetWidth(8.8)
    graph.obj.SetHeight(5.5)

    layer.obj.LT_execute("legend -r;doc -uw;")
    legend = layer.label("legend") or layer.label("Legend")
    if legend is None:
        raise RuntimeError("J legend is missing")
    legend.text = "\\l(1) Center (r/R = 0-0.4)   \\l(2) Outer edge (r/R = 0.8-1.0)"
    legend.set_int("attach", 0)
    legend.set_int("font", FONT_CODE)
    legend.set_float("fsize", 7.8)
    legend.set_int("frame", 0)
    legend.set_int("showframe", 0)
    graph.activate()
    page_width = float(op.lt_float("page.width"))
    page_height = float(op.lt_float("page.height"))
    legend.set_float("left", 0.18 * page_width)
    legend.set_float("top", 0.025 * page_height)
    add_page_label(layer, "panel_j", "j", 7.5, 5.0, 13.0)
    layer9 = summary.loc[summary["layer"] == 9].iloc[0]
    add_page_label(
        layer,
        "JDeltaL9",
        f"Delta = {layer9['delta_center_minus_outer']:+.3f}",
        0.59 * page_width,
        0.23 * page_height,
        7.5,
    )

    readback = plot_sheet.to_df()[plot_frame.columns.tolist()]
    assert_numeric_equal(plot_frame, readback, "J plot data")
    return graph, {
        "panel": "j",
        "source_rows": int(len(summary)),
        "full_audit_rows": int(len(audit)),
        "scientific_metric": "candidate-count-weighted measured lattice support fraction",
        "display": "compact Center and Outer-edge raw curves; candidate denominators retained in workbook",
        "layer9_delta": float(layer9["delta_center_minus_outer"]),
        "layer12_delta": float(
            summary.loc[summary["layer"] == 12, "delta_center_minus_outer"].iloc[0]
        ),
        "candidate_counts_visible_in_graph": False,
        "candidate_counts_retained_in_origin_workbook": True,
        "input_numeric_hash": numeric_hash(plot_frame),
        "origin_numeric_hash": numeric_hash(readback),
        "source_values_modified": False,
        "smoothing_or_fit": False,
        "x_axis": axis_readback(layer, "x", "J"),
        "y_axis": axis_readback(layer, "y", "J"),
    }


def write_records(h_stats: dict[str, float], reports: list[dict[str, object]]) -> None:
    contract = {
        "core_conclusion": (
            "Ordered-region reciprocal-space evidence is supported by alternating projected-column contrast, "
            "and measured lattice support is heterogeneous across reconstructed layer and radius."
        ),
        "origin_target": "Origin2025b 10.25.212",
        "panel_h_role": "real-space quantitative validation after g1-g3",
        "panel_j_role": "compact quantitative summary derived from i3",
        "j_compaction": {
            "kept": ["Center curve", "Outer-edge curve", "Layer 9 signed difference"],
            "removed_from_visible_panel": ["32 candidate-site labels", "Layer 12 annotation"],
            "retained_in_opju": ["80-cell audit table", "all per-layer denominators", "all signed differences"],
        },
        "forbidden_interpretations": [
            "formal chemical order parameter",
            "atomic chemical identification",
            "phase fraction",
            "crystallinity percentage",
            "composition",
            "time evolution",
        ],
        "source_sha256": {
            "h": sha256_file(H_SOURCE),
            "j": sha256_file(J_SOURCE),
            "j_80_cell_audit": sha256_file(J_AUDIT_SOURCE),
        },
        "reports": reports,
    }
    (OUTPUT / "figure_contract_and_origin_readback.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    readme = f"""# Figure 1 h/j Origin2025b 审计包

## 图的作用

- h：紧接 g1-g3，用固定晶格奇偶分组的实空间柱强度 ECDF 验证 ordered 区域的倒易空间证据。
- j：由 i3 同一批 16 x 5 审计单元得到的 Center/Outer-edge 定量总结，不是独立数据集。

## h 审计结果

- 固定 parity 0/1 样本数：{int(h_stats['n_parity_0'])}/{int(h_stats['n_parity_1'])}。
- 中位数：{h_stats['median_parity_0']:.6f}/{h_stats['median_parity_1']:.6f}；比值 {h_stats['median_ratio']:.3f}。
- 描述性 AUC：{h_stats['descriptive_auc']:.6f}。单颗粒区域内技术样本具有空间相关性，因此不报告 p 值。
- 强度归一化：每个柱的背景扣除积分强度除以其 8 个最近邻柱积分强度的中位数。

## j 紧凑化决定

- 图内仅保留 Center、Outer edge 和 Layer 9 的差值标注。
- 32 个逐层候选位点数字及 Layer 12 的差值不在小面板中显示。
- 所有候选位点数、匹配位点数、逐层差值和 80 单元源表均保留在 OPJU 工作表和 CSV 中。
- 原始点不平滑、不拟合、不做时间解释。

## 文件

- `result.opju`：同时包含 h/j 图页、绘图数据、完整源数据和定义工作表。
- `Fig1_h_origin.*` 与 `Fig1_j_origin_compact.*`：Origin 导出的 PNG/PDF/TIF。
- `figure_contract_and_origin_readback.json`：数值回读、轴设置和科学边界。
"""
    (OUTPUT / "README_CN.md").write_text(readme, encoding="utf-8")


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    for source in (H_SOURCE, J_SOURCE, J_AUDIT_SOURCE):
        shutil.copy2(source, OUTPUT / source.name)

    h_source, h_plot, h_stats = prepare_h()
    j_summary, j_audit = prepare_j()
    source_hashes_before = {
        path: sha256_file(path) for path in (H_SOURCE, J_SOURCE, J_AUDIT_SOURCE)
    }

    op.set_show(False)
    try:
        op.new(asksave=False)
        create_definitions_sheet(h_stats)
        graph_h, report_h = build_h(h_source, h_plot, h_stats)
        graph_j, report_j = build_j(j_summary, j_audit)

        opju = OUTPUT / "result.opju"
        if not op.save(str(opju)):
            raise RuntimeError("Origin failed to save the editable OPJU")
        export_graph(
            op,
            graph_h,
            OUTPUT / "Fig1_h_origin.png",
            OUTPUT / "Fig1_h_origin.pdf",
            OUTPUT / "Fig1_h_origin.tif",
            raster_width=2100,
        )
        export_graph(
            op,
            graph_j,
            OUTPUT / "Fig1_j_origin_compact.png",
            OUTPUT / "Fig1_j_origin_compact.pdf",
            OUTPUT / "Fig1_j_origin_compact.tif",
            raster_width=2300,
        )
    finally:
        op.exit()

    if source_hashes_before != {
        path: sha256_file(path) for path in (H_SOURCE, J_SOURCE, J_AUDIT_SOURCE)
    }:
        raise RuntimeError("An immutable source changed during Origin rendering")
    reports = [report_h, report_j]
    write_records(h_stats, reports)
    artifacts = {}
    for name in (
        "result.opju",
        "Fig1_h_origin.png",
        "Fig1_h_origin.pdf",
        "Fig1_h_origin.tif",
        "Fig1_j_origin_compact.png",
        "Fig1_j_origin_compact.pdf",
        "Fig1_j_origin_compact.tif",
    ):
        path = OUTPUT / name
        artifacts[name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (OUTPUT / "artifact_manifest.json").write_text(
        json.dumps(
            {
                "status": "origin_exports_complete_pending_visual_qa",
                "artifacts": artifacts,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(json.dumps({"output": str(OUTPUT), "panels": ["h", "j"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
