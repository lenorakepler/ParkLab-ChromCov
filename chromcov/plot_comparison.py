"""Render the per-chromosome coverage comparison to an interactive HTML.

Panels:
  1. Mean depth per chromosome, grouped by tool (samtools / pandepth / mosdepth)
  2. Approx. copy number (2 x depth / autosomal baseline), colored by gain/loss
  3. Breadth at depth: % of each chromosome covered at >=1x/10x/20x/30x
     (from the mosdepth global distribution)
  4. The full comparison table

Inputs:
  out/coverage_comparison.csv            (from chromcov/__coverage.py build_comparison)
  data/testmos.mosdepth.global.dist.txt  (mosdepth cumulative depth distribution)
"""
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DATA_DIR = Path().resolve() / "data"
OUT_DIR = Path().resolve() / "out"

COMPARISON_CSV = OUT_DIR / "coverage_comparison.csv"
MOSDEPTH_DIST = DATA_DIR / "testmos.mosdepth.global.dist.txt"
HTML_OUT = DATA_DIR / "coverage_comparison.html"

TOOLS = [("samtools", "#3b6ea5"), ("pandepth", "#27ae60"), ("mosdepth", "#e67e22"),
         ("idx_depth", "#9b59b6")]
THRESHOLDS = [1, 10, 20, 30]


def breadth_at_depth(chroms, dist_file=MOSDEPTH_DIST, thresholds=THRESHOLDS):
	"""For each chromosome, % of bases covered at >= each threshold, from the
	mosdepth global.dist (chrom, depth, proportion-at->=-depth)."""
	d = pd.read_csv(dist_file, sep="\t", names=["chrom", "depth", "prop"])
	out = {}
	for c in chroms:
		lut = d[d.chrom == c].set_index("depth")["prop"]
		out[c] = [100 * float(lut.get(t, 0.0)) for t in thresholds]
	return out


def build_figure():
	df = pd.read_csv(COMPARISON_CSV)
	chroms = df["chrom"].tolist()
	cn_color = ["#c0392b" if r > 1.25 else "#2c7fb8" if r < 0.75 else "#7f8c8d"
	            for r in df["ratio"]]
	breadth = breadth_at_depth(chroms)

	fig = make_subplots(
		rows=4, cols=1, vertical_spacing=0.055,
		row_heights=[0.26, 0.24, 0.24, 0.26],
		specs=[[{"type": "xy"}], [{"type": "xy"}], [{"type": "xy"}], [{"type": "table"}]],
		subplot_titles=("Mean depth per chromosome by tool",
		                "Approx. copy number (2 x depth / autosomal baseline)",
		                "Breadth at depth: % of chromosome covered at >= Nx",
		                "Comparison table"))

	# Panel 1: grouped mean-depth bars by tool (only those present)
	for tool, color in TOOLS:
		if tool in df.columns:
			fig.add_bar(x=chroms, y=df[tool], name=tool, marker_color=color, row=1, col=1)

	# Panel 2: approx CN colored by gain/loss, with indexcov overlaid as an
	# independent (index-only) CN estimate (scaled coverage x2 ~= copy number).
	fig.add_bar(x=chroms, y=df["approx_cn"], marker_color=cn_color, name="approx_cn (depth)",
	            showlegend=True, hovertemplate="%{x}<br>CN %{y:.2f}<extra></extra>", row=2, col=1)
	if "indexcov" in df.columns:
		fig.add_scatter(x=chroms, y=2 * df["indexcov"], mode="markers", name="indexcov x2",
		                marker=dict(color="black", size=7, symbol="diamond"),
		                hovertemplate="%{x}<br>indexcov CN %{y:.2f}<extra></extra>", row=2, col=1)
	fig.add_hline(y=2, line=dict(color="black", dash="dash", width=1), row=2, col=1,
	              annotation_text="CN = 2 (diploid)")

	# Panel 3: breadth at depth, grouped bars per threshold
	shades = {1: "#bdd7e7", 10: "#6baed6", 20: "#3182bd", 30: "#08519c"}
	for i, t in enumerate(THRESHOLDS):
		fig.add_bar(x=chroms, y=[breadth[c][i] for c in chroms], name=f">={t}x",
		            marker_color=shades[t], legendgroup="breadth", row=3, col=1)

	# Panel 4: table — every comparison column that's present, in a sensible order
	preferred = ["samtools", "pandepth", "mosdepth", "idx_depth", "indexcov",
	             "idx_ratio", "ratio", "approx_cn"]
	table_cols = [c for c in preferred if c in df.columns]
	cells = [df["chrom"]] + [df[c].round(3) for c in table_cols]
	fig.add_trace(go.Table(
		header=dict(values=["chrom"] + table_cols,
		            fill_color="#2c3e50", font=dict(color="white"), align="left"),
		cells=dict(values=cells,
		           fill_color=[["#f7f9fb", "#eef2f6"] * len(df)], align="left")),
		row=4, col=1)

	fig.update_yaxes(title_text="mean depth (x)", row=1, col=1)
	fig.update_yaxes(title_text="approx. CN", row=2, col=1)
	fig.update_yaxes(title_text="% covered", range=[0, 100], row=3, col=1)
	fig.update_layout(title="COLO829T per-chromosome coverage — tool comparison",
	                  template="plotly_white", height=1250, width=1250, barmode="group",
	                  legend=dict(orientation="h", y=1.05))
	return fig


if __name__ == "__main__":
	fig = build_figure()
	fig.write_html(HTML_OUT, include_plotlyjs="cdn")
	print(f"wrote {HTML_OUT}")
