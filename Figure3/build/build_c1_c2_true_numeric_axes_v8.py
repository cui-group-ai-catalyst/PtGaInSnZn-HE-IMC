from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import originpro as op

from origin_sciplot.origin_backend.export_utils import export_graph


ROOT = Path(
    r"C:\Users\13102\ns_review_4dstem\DEEPSEEK_FINAL_v2_1\FIG1C_FFT_RAW_ORIGIN_v1"
)
PANELS = {
    "C1": {
        "source": ROOT / "origin_exports_magma_axes_FINAL_v7" / "C1",
        "output": ROOT / "C1_true_numeric_axes_FINAL_v8b",
    },
    "C2": {
        "source": ROOT / "origin_exports_magma_axes_FINAL_v7" / "C2",
        "output": ROOT / "C2_true_numeric_axes_FINAL_v8b",
    },
}
MAJOR_CANDIDATES = (-10.0, -5.0, 0.0, 5.0, 10.0)
MINOR_CANDIDATES = (-7.5, -2.5, 2.5, 7.5)


def array_hash(values: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(values, dtype=np.float64)
    return hashlib.sha256(contiguous.tobytes()).hexdigest()


def tick_text(values: tuple[float, ...]) -> str:
    return " ".join(f"{value:g}" for value in values)


def axis_state(layer, axis: str) -> dict[str, float | int | str]:
    return {
        "from": float(layer.get_float(f"{axis}.from")),
        "to": float(layer.get_float(f"{axis}.to")),
        "increment": float(layer.get_float(f"{axis}.inc")),
        "first_tick": float(layer.get_float(f"{axis}.firstTick")),
        "ticks": int(layer.get_int(f"{axis}.ticks")),
        "minor_ticks": int(layer.get_int(f"{axis}.minorTicks")),
        "major_length": float(layer.get_float(f"{axis}.ticklength")),
        "minor_length": float(layer.get_float(f"{axis}.mticklength")),
        "major_thickness": float(layer.get_float(f"{axis}.tickthickness")),
        "minor_thickness": float(layer.get_float(f"{axis}.mtickthickness")),
        "label_type": int(layer.get_int(f"{axis}.label.type")),
        "label_size_pt": float(layer.get_float(f"{axis}.label.pt")),
        "label_font_code": int(layer.get_int(f"{axis}.label.font")),
        "ticks_by_data": op.get_lt_str(f"layer.{axis}.ticksbydata$"),
        "minor_ticks_by_data": op.get_lt_str(f"layer.{axis}.mticksbydata$"),
    }


def build_panel(panel: str, source_dir: Path, output_dir: Path) -> dict[str, object]:
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing directory: {output_dir}")
    output_dir.mkdir(parents=True)

    base_copy = output_dir / "base_v7_copy.opju"
    shutil.copy2(source_dir / "result.opju", base_copy)
    shutil.copy2(source_dir / "input_copy.csv", output_dir / "input_copy.csv")
    shutil.copy2(source_dir / "render-plan.json", output_dir / "source_render-plan.json")
    shutil.copy2(
        source_dir / "origin_verify_report.json",
        output_dir / "source_v7_origin_verify_report.json",
    )

    frame = pd.read_csv(source_dir / "input_copy.csv")
    qx = np.asarray([float(value) for value in frame.columns[1:]], dtype=float)
    qy = frame.iloc[:, 0].to_numpy(dtype=float)
    if not np.all(np.diff(qx) > 0):
        raise ValueError(f"{panel}: qx coordinates must be strictly increasing")
    if not np.all(np.diff(qy) < 0):
        raise ValueError(f"{panel}: qy coordinates must be strictly decreasing")

    lower = max(float(qx.min()), float(qy.min()))
    upper = min(float(qx.max()), float(qy.max()))
    major_values = tuple(value for value in MAJOR_CANDIDATES if lower < value < upper)
    minor_values = tuple(value for value in MINOR_CANDIDATES if lower < value < upper)
    if not major_values or not minor_values:
        raise ValueError(f"{panel}: requested ticks do not fit the reciprocal-space range")
    major_text = tick_text(major_values)
    minor_text = tick_text(minor_values)

    op.set_show(False)
    try:
        if not op.open(str(base_copy), readonly=False, asksave=False):
            raise RuntimeError(f"Origin could not open the {panel} v7 copy")
        graph = op.find_graph(0)
        matrix = op.find_sheet("m", 0)
        if graph is None or matrix is None:
            raise RuntimeError(f"{panel}: graph or matrix sheet is missing")
        layer = graph[0]
        plot = layer.plot_list()[0]

        values_before = matrix.to_np2d().copy()
        hash_before = array_hash(values_before)
        xymap_before = tuple(float(value) for value in matrix.xymap)

        matrix.xymap = (float(qx[0]), float(qx[-1]), float(qy[0]), float(qy[-1]))
        layer.axis("x").set_limits(float(qx[0]), float(qx[-1]), 5.0)
        layer.axis("y").set_limits(float(qy[-1]), float(qy[0]), 5.0)

        for axis in ("x", "y"):
            layer.set_int(f"{axis}.label.type", 1)
            layer.set_int(f"{axis}.label.numFormat", 1)
            layer.set_int(f"{axis}.label.decimalPlaces", 0)
            layer.set_float(f"{axis}.label.rotate", 0.0)
            layer.set_int(f"{axis}.ticks", 10)  # major out (2) + minor out (8)
            layer.set_int(f"{axis}.minorTicks", 0)
            layer.set_float(f"{axis}.ticklength", 5.5)
            layer.set_float(f"{axis}.mticklength", 3.0)
            layer.set_float(f"{axis}.tickthickness", 1.8)
            layer.set_float(f"{axis}.mtickthickness", 1.2)
            layer.set_float(f"{axis}.firstTick", major_values[0])
            layer.set_int(f"{axis}.showLabels", 1)
            layer.set_int(f"{axis}.label.show", 1)

        if not layer.obj.LT_execute(
            f'layer.x.ticksbydata$="{major_text}";'
            f'layer.y.ticksbydata$="{major_text}";'
            f'layer.x.mticksbydata$="{minor_text}";'
            f'layer.y.mticksbydata$="{minor_text}";'
            "doc -uw;"
        ):
            raise RuntimeError(f"{panel}: Origin could not apply exact tick positions")

        values_after = matrix.to_np2d().copy()
        hash_after = array_hash(values_after)
        if hash_after != hash_before or not np.array_equal(
            values_after, values_before, equal_nan=True
        ):
            raise RuntimeError(f"{panel}: matrix values changed while setting coordinate metadata")

        graph.activate()
        op.set_show(True)
        op.lt_exec("doc -uw;")
        palette = plot.colormap
        output_opju = output_dir / "result.opju"
        if not op.save(str(output_opju)):
            raise RuntimeError(f"Origin could not save {panel} v8")
        exports = export_graph(
            op,
            graph,
            output_dir / "result.png",
            output_dir / "result.pdf",
            output_dir / "result.tif",
        )

        report = {
            "panel": panel,
            "source_project": str(source_dir / "result.opju"),
            "output_project": str(output_opju),
            "matrix_shape": list(values_after.shape),
            "matrix_hash_before": hash_before,
            "matrix_hash_after": hash_after,
            "matrix_values_modified": False,
            "xymap_before": list(xymap_before),
            "xymap_after": [float(value) for value in matrix.xymap],
            "major_tick_values": list(major_values),
            "minor_tick_values": list(minor_values),
            "tick_direction": "out",
            "palette": palette,
            "x_axis": axis_state(layer, "x"),
            "y_axis": axis_state(layer, "y"),
            "exports": exports,
            "source_data_modified": False,
        }
        report_text = json.dumps(report, ensure_ascii=False, indent=2)
        (output_dir / "numeric_axis_readback.json").write_text(report_text, encoding="utf-8")
        return report
    finally:
        op.exit()


reports = []
for panel_name, paths in PANELS.items():
    reports.append(build_panel(panel_name, paths["source"], paths["output"]))
print(json.dumps(reports, ensure_ascii=False, indent=2))
