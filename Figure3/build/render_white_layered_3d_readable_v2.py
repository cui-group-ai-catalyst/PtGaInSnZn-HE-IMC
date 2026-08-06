"""Readable white-background 3-D rendering for the approved layer stack.

This is a non-destructive presentation upgrade. It reads the locked detected-peak
CSV and layer geometry, keeps every coordinate and default camera setting unchanged,
and writes a separate output package with larger, darker axis text and editable
style controls.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import plotly.graph_objects as go

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",  # keep SVG labels as editable text
    "pdf.fonttype": 42,       # keep PDF labels as TrueType text
})


import fig3_paths

SOURCE_DIR = fig3_paths.SOURCE_DATA
CIRCLES_FILE = fig3_paths.DATA / "per_layer_circles.csv"
SCRIPT_DIR = Path(__file__).resolve().parent
OUT = fig3_paths.OUTPUT
STATIC = OUT / "static"
INTERACTIVE = OUT / "interactive"
ORIGIN = OUT / "origin_ready"
RECORDS = OUT / "records"
COPIED_SOURCE = OUT / "source_data"

A_FILE = SOURCE_DIR / "3d_A_displayed_detected_peaks.csv"
B_FILE = SOURCE_DIR / "3d_B_sampled_intensity_points.csv"
PX_NM = 0.01138848395
COMPACT_DZ = 0.50
EXPLODED_DZ = 0.90
EXPECTED_A_SHA256 = "c158e6b6421db8f30d8142b29362a2fbdfae9516e617277ca103ab29453be7f8"
EXPECTED_B_SHA256 = "9d96fd2a4f48fb9fd5c4a26ee9740223fa798ba4478608553b6ec708c9f613b4"

# Edit these values to change the static figure without touching the data.
STYLE = {
    "text": "#1F2933",
    "axis": "#16212B",
    "grid": "#B7C0C9",
    "pane": "#F8FAFC",
    "pane_edge": "#66727D",
    "ring": "#6B747D",
    "axis_label_size": 11.5,
    "tick_size": 10.0,
    "title_size": 13.0,
    "colorbar_label_size": 10.0,
    "colorbar_tick_size": 9.0,
    "axis_line_width": 1.35,
    "tick_width": 1.05,
    "grid_width": 0.65,
    "projection": "ortho",
    "elev": 17,
    "azim": 40,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_all(fig, stem: Path, dpi: int = 400) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    export_kwargs = {"bbox_inches": "tight", "pad_inches": 0.24, "facecolor": "white"}
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, **export_kwargs)
    fig.savefig(stem.with_suffix(".pdf"), **export_kwargs)
    fig.savefig(stem.with_suffix(".svg"), **export_kwargs)
    plt.close(fig)


def style_white_3d(ax, zmax: float) -> None:
    ax.set_proj_type(STYLE["projection"])
    ax.set_facecolor("white")
    ax.view_init(elev=STYLE["elev"], azim=STYLE["azim"])
    ax.set_xlim(4.0, 11.5)
    ax.set_ylim(4.0, 11.5)
    ax.set_zlim(0.4, zmax)
    ax.set_xticks(np.arange(4, 12, 1))
    ax.set_yticks(np.arange(4, 12, 1))
    ax.set_zticks(np.arange(1, int(np.floor(zmax)) + 1, 1))
    ax.set_xlabel("x (nm)", labelpad=10, fontsize=STYLE["axis_label_size"],
                  color=STYLE["text"], fontweight="medium")
    ax.set_ylabel("y (nm)", labelpad=10, fontsize=STYLE["axis_label_size"],
                  color=STYLE["text"], fontweight="medium")
    ax.set_zlabel("display depth (a.u.)", labelpad=3,
                  fontsize=STYLE["axis_label_size"], color=STYLE["text"],
                  fontweight="medium")
    ax.tick_params(colors=STYLE["text"], labelsize=STYLE["tick_size"],
                   width=STYLE["tick_width"], pad=2)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = True
        axis.pane.set_facecolor(STYLE["pane"])
        axis.pane.set_edgecolor(STYLE["pane_edge"])
        axis._axinfo["grid"].update({
            "color": STYLE["grid"],
            "linewidth": STYLE["grid_width"],
            "linestyle": "-",
        })
        axis.line.set_color(STYLE["axis"])
        axis.line.set_linewidth(STYLE["axis_line_width"])
    ax.set_box_aspect((1, 1, 1.18))


def add_rings(ax, circles: pd.DataFrame, dz: float) -> None:
    theta = np.linspace(0, 2 * np.pi, 360)
    for layer in (1, 3, 5, 8, 10, 12, 14, 16):
        row = circles.loc[layer]
        radius = float(row.display_radius_px) * PX_NM
        cx = float(row.center_x_px) * PX_NM
        cy = float(row.center_y_px) * PX_NM
        ax.plot(cx + radius * np.cos(theta), cy + radius * np.sin(theta),
                np.full_like(theta, (17 - layer) * dz), color=STYLE["ring"],
                linewidth=0.75, alpha=0.72)


def render_static(a: pd.DataFrame, circles: pd.DataFrame, lo: float, hi: float,
                  dz: float, name: str, title: str, alpha: float) -> None:
    fig = plt.figure(figsize=(8.1, 7.7), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")
    z = (17 - a.layer.to_numpy()) * dz
    scatter = ax.scatter(a.x_nm, a.y_nm, z, c=a.brightness, cmap="viridis",
                         norm=Normalize(lo, hi), s=7.5, marker="o", alpha=alpha,
                         linewidths=0, depthshade=True, rasterized=True)
    add_rings(ax, circles, dz)
    style_white_3d(ax, 16.6 * dz)
    ax.set_title(title, pad=12, color=STYLE["text"], fontsize=STYLE["title_size"],
                 fontweight="medium")
    colorbar = fig.colorbar(scatter, ax=ax, shrink=0.60, pad=0.08, aspect=24)
    colorbar.set_label("fitted peak brightness (a.u.)",
                       fontsize=STYLE["colorbar_label_size"], color=STYLE["text"])
    colorbar.ax.tick_params(labelsize=STYLE["colorbar_tick_size"],
                            width=STYLE["tick_width"], colors=STYLE["text"])
    colorbar.outline.set_linewidth(0.75)
    colorbar.outline.set_edgecolor(STYLE["axis"])
    save_all(fig, STATIC / name)


def render_layer_montage(a: pd.DataFrame, circles: pd.DataFrame, lo: float, hi: float) -> None:
    fig, axes = plt.subplots(4, 4, figsize=(7.8, 7.8), facecolor="white")
    last = None
    for layer, ax in enumerate(axes.flat, start=1):
        one = a[a.layer == layer]
        last = ax.scatter(one.x_nm, one.y_nm, c=one.brightness, cmap="viridis",
                          norm=Normalize(lo, hi), s=3.2, linewidths=0)
        row = circles.loc[layer]
        rr = float(row.display_radius_px) * PX_NM
        cx = float(row.center_x_px) * PX_NM
        cy = float(row.center_y_px) * PX_NM
        ax.add_patch(plt.Circle((cx, cy), rr, fill=False, color=STYLE["ring"], linewidth=0.55))
        ax.set_aspect("equal")
        ax.set_xlim(4.0, 11.5)
        ax.set_ylim(4.0, 11.5)
        ax.set_title(f"Layer {layer}", pad=2, color=STYLE["text"], fontsize=8.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(STYLE["pane_edge"])
            spine.set_linewidth(0.55)
    fig.subplots_adjust(left=0.035, right=0.90, bottom=0.035, top=0.965,
                        wspace=0.10, hspace=0.18)
    cax = fig.add_axes([0.925, 0.25, 0.018, 0.50])
    cb = fig.colorbar(last, cax=cax)
    cb.set_label("fitted peak brightness (a.u.)", fontsize=8.5, color=STYLE["text"])
    cb.ax.tick_params(labelsize=7.5, width=STYLE["tick_width"], colors=STYLE["text"])
    cb.outline.set_linewidth(0.6)
    save_all(fig, STATIC / "A_white_each_layer_4x4_readable")


def make_interactive(a: pd.DataFrame, circles: pd.DataFrame, lo: float, hi: float) -> None:
    fig = go.Figure()
    point_counts = []
    for layer in range(1, 17):
        one = a[a.layer == layer]
        point_counts.append(int(len(one)))
        fig.add_trace(go.Scatter3d(
            x=one.x_nm, y=one.y_nm,
            z=np.full(len(one), (17 - layer) * EXPLODED_DZ),
            mode="markers", name=f"Layer {layer}", showlegend=False,
            customdata=np.column_stack([one.layer, one.brightness]),
            hovertemplate=("Layer %{customdata[0]}<br>brightness=%{customdata[1]:.1f}"
                           "<br>x=%{x:.2f} nm<br>y=%{y:.2f} nm<extra></extra>"),
            marker=dict(size=4.2, symbol="circle", color=one.brightness,
                        colorscale="Viridis", cmin=lo, cmax=hi, opacity=0.82,
                        line=dict(width=0), showscale=(layer == 1),
                        colorbar=dict(
                            title=dict(text="Fitted peak<br>brightness (a.u.)",
                                        font=dict(size=14, color=STYLE["text"])),
                            tickfont=dict(size=12, color=STYLE["text"]),
                            thickness=16, len=0.64, x=1.01, y=0.50,
                        )),
        ))

    theta = np.linspace(0, 2 * np.pi, 300)
    for layer in range(1, 17):
        row = circles.loc[layer]
        rr = float(row.display_radius_px) * PX_NM
        cx = float(row.center_x_px) * PX_NM
        cy = float(row.center_y_px) * PX_NM
        fig.add_trace(go.Scatter3d(
            x=cx + rr * np.cos(theta), y=cy + rr * np.sin(theta),
            z=np.full_like(theta, (17 - layer) * EXPLODED_DZ),
            mode="lines", hoverinfo="skip", showlegend=False,
            line=dict(color="rgba(107,116,125,0.72)", width=2),
        ))

    axis_common = dict(
        color=STYLE["axis"], gridcolor=STYLE["grid"], linecolor=STYLE["axis"],
        linewidth=3, zeroline=False, showline=True, showbackground=True,
        backgroundcolor=STYLE["pane"], showspikes=False,
        tickfont=dict(size=13, color=STYLE["text"]),
    )
    axis_title = dict(font=dict(size=16, color=STYLE["text"]))
    fig.update_layout(
        paper_bgcolor="white", plot_bgcolor="white", autosize=True,
        margin=dict(l=4, r=46, t=52, b=4),
        title=dict(text="A | Layer-resolved detected peaks", x=0.5,
                   font=dict(size=20, color=STYLE["text"])),
        scene=dict(
            xaxis=dict(axis_common, title=dict(text="x (nm)", **axis_title), range=[4, 11.5]),
            yaxis=dict(axis_common, title=dict(text="y (nm)", **axis_title), range=[4, 11.5]),
            zaxis=dict(axis_common, title=dict(text="display depth (a.u.)", **axis_title),
                       range=[0.7, 15.0]),
            aspectmode="manual", aspectratio=dict(x=1, y=1, z=1.35),
            camera=dict(projection=dict(type="orthographic"),
                        eye=dict(x=1.65, y=1.35, z=0.78)),
        ),
        showlegend=False,
    )
    html = fig.to_html(
        include_plotlyjs=True, full_html=True, div_id="layered3d_readable",
        config=dict(displaylogo=False, responsive=True, scrollZoom=True,
                    # The custom Export PNG button preserves the current on-screen sizes.
                    modeBarButtonsToRemove=["lasso3d", "select2d", "toImage"]),
    )
    controls = """
<div id="controls" role="region" aria-label="3D display controls">
  <label>Layer spacing <input id="spacing" type="range" min="0.35" max="1.40" step="0.05" value="0.90"><output id="spacingValue">0.90</output></label>
  <label>Projection <select id="projection"><option value="orthographic">Orthographic</option><option value="perspective">Perspective</option></select></label>
  <label>Visible layer <select id="layer"><option value="all">All layers</option></select></label>
  <label>Opacity <input id="opacity" type="range" min="0.20" max="1.00" step="0.05" value="0.82"></label>
  <label>Point size <input id="pointSize" type="range" min="1.5" max="8" step="0.25" value="4.25"></label>
  <label>Axis titles <input id="axisTitleSize" type="range" min="10" max="24" step="1" value="16"></label>
  <label>Axis ticks <input id="axisTickSize" type="range" min="8" max="22" step="1" value="13"></label>
  <label>Axis colour <input id="axisColor" type="color" value="#16212B"></label>
  <label>Grid colour <input id="gridColor" type="color" value="#B7C0C9"></label>
  <label><input id="axes" type="checkbox" checked> Axes</label>
  <label><input id="grid" type="checkbox" checked> Grid</label>
  <label><input id="rings" type="checkbox" checked> Layer guides</label>
  <label><input id="bar" type="checkbox" checked> Colour bar</label>
  <label>Bar x <input id="barX" type="range" min="0.80" max="1.10" step="0.01" value="1.01"></label>
  <label>Bar y <input id="barY" type="range" min="0.20" max="0.80" step="0.02" value="0.50"></label>
  <label>Bar title <input id="barTitleSize" type="range" min="9" max="24" step="1" value="14"></label>
  <label>Bar ticks <input id="barTickSize" type="range" min="8" max="20" step="1" value="12"></label>
  <button id="reset" type="button">Reset view</button>
  <button id="export" type="button">Export PNG</button>
</div>
"""
    styles = """
<style>
:root { color-scheme:light; --bg:#FFFFFF; --panel:#F3F5F7; --line:#C9CFD5; --text:#20262C; --muted:#46515B; --focus:#1261A0; }
html,body { width:100%; height:100%; margin:0; overflow:hidden; background:var(--bg); color:var(--text); font-family:"Segoe UI",Arial,sans-serif; }
body { display:flex; flex-direction:column; }
#controls { box-sizing:border-box; flex:0 0 auto; min-height:76px; padding:10px 16px; display:flex; flex-wrap:wrap; align-items:center; gap:7px 15px; background:var(--panel); border-bottom:1px solid var(--line); }
#controls label { display:flex; align-items:center; gap:7px; min-height:34px; font-size:13px; white-space:nowrap; color:var(--muted); }
#controls input[type="range"] { width:106px; accent-color:#1261A0; }
#controls input[type="checkbox"] { width:18px; height:18px; accent-color:#1261A0; }
#controls input[type="color"] { width:28px; height:24px; padding:1px; border:1px solid #9EA8B1; background:#FFFFFF; }
#controls select,#controls button { min-height:34px; padding:0 11px; color:var(--text); background:#FFFFFF; border:1px solid #9EA8B1; border-radius:4px; }
#controls button { cursor:pointer; }
#controls button:hover { background:#E5E9ED; }
#controls :focus-visible { outline:2px solid var(--focus); outline-offset:2px; }
#spacingValue { min-width:34px; color:var(--text); }
#layered3d_readable { box-sizing:border-box; flex:1 1 auto; width:100% !important; max-width:100vw; min-width:0; min-height:300px; height:auto !important; }
@media (max-width:760px) { #controls { max-height:330px; overflow-y:auto; overflow-x:hidden; gap:4px 11px; padding:8px 12px; } #controls label { min-height:31px; font-size:12px; } #controls input[type="range"] { width:102px; } }
</style>
"""
    js = f"""
<script>
(() => {{
  const gd = document.getElementById('layered3d_readable');
  const controls = document.getElementById('controls');
  const pointTraces = [...Array(16).keys()];
  const pointLayers = {json.dumps(list(range(1, 17)))};
  const pointCounts = {json.dumps(point_counts)};
  const ringStart = 16;
  const ringLayers = {json.dumps(list(range(1, 17)))};
  const overviewRings = [1,3,5,8,10,12,14,16];
  const spacing = document.getElementById('spacing');
  const spacingValue = document.getElementById('spacingValue');
  const layerSelect = document.getElementById('layer');
  const projection = document.getElementById('projection');
  const axes = document.getElementById('axes');
  const grid = document.getElementById('grid');
  const rings = document.getElementById('rings');
  const bar = document.getElementById('bar');
  pointLayers.forEach(v => layerSelect.add(new Option('Layer ' + v, String(v))));
  const axisNames = ['xaxis','yaxis','zaxis'];
  const applyAxis = (key, value) => {{ const u = {{}}; axisNames.forEach(n => u['scene.'+n+'.'+key] = value); Plotly.relayout(gd, u); }};
  function updateSpacing() {{
    const dz = Number(spacing.value); spacingValue.value = dz.toFixed(2);
    pointLayers.forEach((layer, i) => Plotly.restyle(gd, {{z:[Array(pointCounts[i]).fill((17-layer)*dz)]}}, [i]));
    ringLayers.forEach((layer, i) => Plotly.restyle(gd, {{z:[Array(300).fill((17-layer)*dz)]}}, [ringStart+i]));
    Plotly.relayout(gd, {{'scene.zaxis.range':[Math.max(0.2,dz*0.75),16.6*dz], 'scene.aspectratio.z':Math.min(1.55,Math.max(0.72,1.35*dz/0.9))}});
  }}
  function updateLayer() {{
    const selected = layerSelect.value;
    pointLayers.forEach((layer, i) => Plotly.restyle(gd, {{visible:selected === 'all' || Number(selected) === layer}}, [i]));
    ringLayers.forEach((layer, i) => Plotly.restyle(gd, {{visible:rings.checked && (selected === 'all' ? overviewRings.includes(layer) : Number(selected) === layer)}}, [ringStart+i]));
    pointLayers.forEach((layer, i) => Plotly.restyle(gd, {{'marker.showscale':bar.checked && (selected === 'all' ? i === 0 : Number(selected) === layer)}}, [i]));
  }}
  function resizePlot() {{
    const mobile = window.innerWidth <= 760;
    gd.style.setProperty('height', Math.max(300, window.innerHeight-controls.offsetHeight)+'px', 'important');
    Plotly.relayout(gd, {{'title.text':mobile ? 'A | Layer-resolved peaks' : 'A | Layer-resolved detected peaks', 'title.font.size':mobile ? 15 : 20, 'margin.r':mobile ? 4 : 46, 'margin.t':mobile ? 36 : 52}}).then(() => Plotly.Plots.resize(gd));
  }}
  spacing.addEventListener('input', updateSpacing);
  projection.addEventListener('change', () => Plotly.relayout(gd, {{'scene.camera.projection.type':projection.value}}));
  layerSelect.addEventListener('change', updateLayer);
  document.getElementById('opacity').addEventListener('input', e => Plotly.restyle(gd, {{'marker.opacity':Number(e.target.value)}}, pointTraces));
  document.getElementById('pointSize').addEventListener('input', e => Plotly.restyle(gd, {{'marker.size':Number(e.target.value)}}, pointTraces));
  document.getElementById('axisTitleSize').addEventListener('input', e => applyAxis('title.font.size', Number(e.target.value)));
  document.getElementById('axisTickSize').addEventListener('input', e => applyAxis('tickfont.size', Number(e.target.value)));
  document.getElementById('axisColor').addEventListener('input', e => {{ applyAxis('color', e.target.value); applyAxis('linecolor', e.target.value); applyAxis('title.font.color', e.target.value); applyAxis('tickfont.color', e.target.value); }});
  document.getElementById('gridColor').addEventListener('input', e => applyAxis('gridcolor', e.target.value));
  axes.addEventListener('change', () => Plotly.relayout(gd, {{'scene.xaxis.visible':axes.checked,'scene.yaxis.visible':axes.checked,'scene.zaxis.visible':axes.checked}}));
  grid.addEventListener('change', () => Plotly.relayout(gd, {{'scene.xaxis.showgrid':grid.checked,'scene.yaxis.showgrid':grid.checked,'scene.zaxis.showgrid':grid.checked}}));
  rings.addEventListener('change', updateLayer); bar.addEventListener('change', updateLayer);
  document.getElementById('barX').addEventListener('input', e => Plotly.restyle(gd, {{'marker.colorbar.x':Number(e.target.value)}}, pointTraces));
  document.getElementById('barY').addEventListener('input', e => Plotly.restyle(gd, {{'marker.colorbar.y':Number(e.target.value)}}, pointTraces));
  document.getElementById('barTitleSize').addEventListener('input', e => Plotly.restyle(gd, {{'marker.colorbar.title.font.size':Number(e.target.value)}}, pointTraces));
  document.getElementById('barTickSize').addEventListener('input', e => Plotly.restyle(gd, {{'marker.colorbar.tickfont.size':Number(e.target.value)}}, pointTraces));
  document.getElementById('reset').addEventListener('click', () => Plotly.relayout(gd, {{'scene.camera':{{projection:{{type:projection.value}},eye:{{x:1.65,y:1.35,z:0.78}}}}}}));
  document.getElementById('export').addEventListener('click', async () => {{
    const button = document.getElementById('export');
    button.disabled = true;
    try {{
      // Keep the exported proportions tied to the adjusted on-screen plot.
      await Plotly.Plots.resize(gd);
      await new Promise(resolve => requestAnimationFrame(resolve));
      const box = gd.getBoundingClientRect();
      const width = Math.max(900, Math.round(box.width || 1000));
      const height = Math.max(650, Math.round(box.height || 800));
      const imageData = await Plotly.toImage(gd, {{format:'png', filename:'A_white_layer_resolved_readable', width, height, scale:2}});
      const link = document.createElement('a');
      link.href = imageData;
      link.download = 'A_white_layer_resolved_readable.png';
      document.body.appendChild(link);
      link.click();
      link.remove();
    }} finally {{
      button.disabled = false;
    }}
  }});
  window.addEventListener('resize', resizePlot); new ResizeObserver(resizePlot).observe(controls);
  updateSpacing(); updateLayer(); resizePlot();
}})();
</script>
"""
    html = html.replace("</head>", styles + "</head>", 1)
    html = html.replace("<body>", "<body>" + controls, 1)
    html = html.replace("</body>", js + "</body>", 1)
    INTERACTIVE.mkdir(parents=True, exist_ok=True)
    (INTERACTIVE / "A_white_layer_resolved_readable_interactive.html").write_text(html, encoding="utf-8")


def write_origin_files(a: pd.DataFrame, b: pd.DataFrame) -> None:
    ORIGIN.mkdir(parents=True, exist_ok=True)
    a_out = a.copy()
    a_out["z_compact"] = (17 - a_out.layer) * COMPACT_DZ
    a_out["z_exploded"] = (17 - a_out.layer) * EXPLODED_DZ
    a_out = a_out[["x_nm", "y_nm", "z_compact", "z_exploded", "layer", "brightness"]]
    a_out.columns = ["X_nm", "Y_nm", "Z_compact", "Z_exploded", "Layer", "Fitted_Peak_Brightness_au"]
    a_out.to_csv(ORIGIN / "Origin_A_detected_peaks_readable.csv", index=False, encoding="utf-8-sig")
    b_out = b[["x_nm", "y_nm", "display_z", "layer", "sampled_intensity", "source"]].copy()
    b_out.columns = ["X_nm", "Y_nm", "Z_display", "Layer", "Sampled_Reconstruction_Intensity_au", "Source"]
    b_out.to_csv(ORIGIN / "Origin_B_intensity_cloud.csv", index=False, encoding="utf-8-sig")
    guide = f"""# Origin drawing guide: readable white 3-D version

Import `Origin_A_detected_peaks_readable.csv` and plot `X_nm`, `Y_nm`, and either
`Z_exploded` (layer-resolved display) or `Z_compact` (compact volume). Map symbol
colour to `Fitted_Peak_Brightness_au`.

Recommended presentation settings: white page; orthographic camera; axis and tick
colour `{STYLE['axis']}`; grid `{STYLE['grid']}`; axis-title size {STYLE['axis_label_size']} pt;
tick size {STYLE['tick_size']} pt; axis line width {STYLE['axis_line_width']} pt.
These are display settings only. The colour scale is fitted peak brightness (a.u.),
not phase, composition, or atom-by-atom elemental identity.
"""
    (ORIGIN / "ORIGIN_DRAWING_GUIDE_READABLE.md").write_text(guide, encoding="utf-8")


def main() -> None:
    for folder in (STATIC, INTERACTIVE, ORIGIN, RECORDS, COPIED_SOURCE):
        folder.mkdir(parents=True, exist_ok=True)
    before = {"A": sha256(A_FILE), "B": sha256(B_FILE)}
    if before != {"A": EXPECTED_A_SHA256, "B": EXPECTED_B_SHA256}:
        raise RuntimeError(f"Source hash mismatch: {before}")
    a = pd.read_csv(A_FILE)
    b = pd.read_csv(B_FILE)
    circles = pd.read_csv(CIRCLES_FILE).set_index("layer")
    if len(a) != 10332 or len(b) != 20000:
        raise AssertionError(f"Unexpected source sizes: A={len(a)}, B={len(b)}")
    lo, hi = np.percentile(a.brightness, [2, 98])
    render_static(a, circles, float(lo), float(hi), EXPLODED_DZ,
                  "A_white_exploded_orthographic_readable",
                  "A  layer-resolved detected peaks", 0.78)
    render_static(a, circles, float(lo), float(hi), COMPACT_DZ,
                  "A_white_compact_volume_readable",
                  "A  compact detected-peak volume", 0.86)
    render_layer_montage(a, circles, float(lo), float(hi))
    make_interactive(a, circles, float(lo), float(hi))
    write_origin_files(a, b)
    shutil.copy2(A_FILE, COPIED_SOURCE / A_FILE.name)
    shutil.copy2(B_FILE, COPIED_SOURCE / B_FILE.name)
    after = {"A": sha256(A_FILE), "B": sha256(B_FILE)}
    if before != after:
        raise AssertionError("Source file changed during rendering")
    params = {
        "style": STYLE,
        "data": {"A_points": int(len(a)), "B_points": int(len(b)),
                 "A_sha256": before["A"], "B_sha256": before["B"],
                 "brightness_percentile_limits": [float(lo), float(hi)]},
        "defaults": {"projection": "orthographic", "layer_spacing": EXPLODED_DZ,
                     "camera_eye": {"x": 1.65, "y": 1.35, "z": 0.78}},
        "scientific_note": "A contains detected atom-column candidates; colour is fitted peak brightness (a.u.), not elemental identity or a formal order parameter.",
    }
    (RECORDS / "figure_parameters.json").write_text(json.dumps(params, ensure_ascii=False, indent=2), encoding="utf-8")
    (RECORDS / "validation_summary.json").write_text(json.dumps({
        "source_files_unchanged": True, "coordinates_changed": False,
        "A_points": int(len(a)), "B_points": int(len(b)), "source_sha256": before,
        "rendered_outputs": [
            "A_white_compact_volume_readable.png/.svg/.pdf",
            "A_white_exploded_orthographic_readable.png/.svg/.pdf",
            "A_white_each_layer_4x4_readable.png/.svg/.pdf",
            "A_white_layer_resolved_readable_interactive.html",
        ],
        "axis_upgrade": {"axis_color": STYLE["axis"], "grid_color": STYLE["grid"],
                         "axis_label_size_pt": STYLE["axis_label_size"],
                         "tick_size_pt": STYLE["tick_size"]},
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    (RECORDS / "README.md").write_text(
        "# Readable white-background 3-D package\n\n"
        "This package is a display-only upgrade of the approved white-background 3-D view.\n"
        "The point coordinates, brightness values, layer indices, default camera and layer\n"
        "spacing are unchanged. Edit `render_white_layered_3d_readable_v2.py` and the\n"
        "`STYLE` dictionary to regenerate static outputs. Open the interactive HTML to\n"
        "adjust axis title/tick size, axis/grid colours, opacity, point size, layer spacing,\n"
        "projection, guide rings and colour-bar placement. The HTML export uses the current\n"
        "plot dimensions at 2x resolution, so adjusted sizes are preserved. SVG/PDF text\n"
        "remains editable.\n",
        encoding="utf-8",
    )
    script_target = OUT / Path(__file__).name
    if Path(__file__).resolve() != script_target.resolve():
        shutil.copy2(__file__, script_target)
    print(OUT)


if __name__ == "__main__":
    main()
