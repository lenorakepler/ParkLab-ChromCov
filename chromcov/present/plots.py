"""
Plots, read from the saved windowed means (never from a live CRAM recompute) so
styling is decoupled from the expensive coverage pass.

  bar_by_chromosome  -> mean depth per primary chromosome (the headline table).
  scatter_windows    -> windowed copy-ratio along the genome with a cytoband
                        ideogram beneath; intrachromosomal CNV breakpoints show
                        up as step changes between windows.

bar_by_chromosome uses matplotlib (Agg backend, headless). scatter_windows uses
plotly so it can emit both a static PNG (via kaleido) and a self-contained
interactive HTML from one figure.
"""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..categories import STRATUM_ORDER
from .cytobands import BAND_OUTLINE, STAIN_COLORS, chrom_length, load_cytobands

# Primary assembly in karyotypic order; decoys/unplaced omitted from the headline
# plots (their per-base means are multi-mapping artifacts).
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Callability tier -> point color (best -> worst callability), + the fallback
# color when no strata were supplied.
STRATA_COLORS = {"easy": "#2ca02c", "difficult": "#ff7f0e", "extreme": "#d62728"}
UNSTRATIFIED_COLOR = "#4C72B0"


def bar_by_chromosome(chrom_means: dict[str, float], out_path: Path, baseline: float | None = None):
    chroms = [c for c in PRIMARY if c in chrom_means]
    vals = [chrom_means[c] for c in chroms]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(chroms, vals, color="#4C72B0")
    if baseline:
        ax.axhline(baseline, color="crimson", ls="--", lw=1, label=f"autosomal median ({baseline:.1f}x)")
        ax.legend()
    ax.set_ylabel("mean depth (x)")
    ax.set_title("Mean coverage per chromosome")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


def _genome_layout(windows: list[dict], min_easy_frac: float):
    """Lay the primary chromosomes end-to-end on one x-axis.

    Returns (layout, rows_by_chrom) where layout is a list of
    (chrom, offset, span, chrom_len): `offset` is the left edge on the shared
    axis, `span` the width used by this chromosome's windows, and `chrom_len` the
    assembled length (for scaling the ideogram onto that same span).
    """
    layout, rows_by_chrom = [], {}
    offset = 0
    for chrom in PRIMARY:
        rows = [w for w in windows
                if w["chrom"] == chrom and w.get("easy_frac", 1.0) >= min_easy_frac]
        if not rows:
            continue
        span = max(w["start"] for w in rows) or 1
        clen = chrom_length(chrom) or span
        layout.append((chrom, offset, span, clen))
        rows_by_chrom[chrom] = rows
        offset += span + 1
    return layout, rows_by_chrom


def _add_scatter(fig, layout, rows_by_chrom, baseline, ploidy, cap_cn):
    """Copy-number points, grouped into one trace per callability tier (so the
    legend and colors match the strata) on row 1. Returns whether strata were
    present (drives the title/legend)."""
    by_tier: dict[str, dict[str, list]] = {}   # tier -> {x, y, cd}
    for chrom, offset, span, _clen in layout:
        for w in rows_by_chrom[chrom]:
            x = offset + w["start"]
            y = min(ploidy * w["mean"] / baseline, cap_cn) if baseline else 0
            tier = w.get("stratum", "") or ""
            d = by_tier.setdefault(tier, {"x": [], "y": [], "cd": []})
            d["x"].append(x)
            d["y"].append(y)
            d["cd"].append((chrom.replace("chr", ""), w["start"], w["mean"]))

    stratified = any(tier for tier in by_tier)
    if stratified:
        order = [t for t in STRATUM_ORDER if t in by_tier] + \
                [t for t in by_tier if t and t not in STRATUM_ORDER]
    else:
        order = [""]

    for tier in order:
        d = by_tier.get(tier, {"x": [], "y": [], "cd": []})
        color = STRATA_COLORS.get(tier, "#888888") if stratified else UNSTRATIFIED_COLOR
        fig.add_trace(go.Scattergl(
            x=d["x"], y=d["y"], mode="markers",
            name=tier or "windows",
            marker=dict(size=3, opacity=0.5, color=color),
            customdata=d["cd"],
            hovertemplate="chr%{customdata[0]}:%{customdata[1]:,}<br>"
                          "depth %{customdata[2]:.1f}x<br>CN≈%{y:.2f}<extra>%{fullData.name}</extra>",
            showlegend=stratified,
        ), row=1, col=1)
    return stratified


def _add_ideogram(fig, layout):
    """Cytoband ideogram on row 2: one horizontal bar per band, each chromosome
    scaled onto the same x-span its windows occupy on row 1 so bands sit under
    the matching copy-number points."""
    bands = load_cytobands()
    base, widths, colors, cd = [], [], [], []
    for chrom, offset, span, clen in layout:
        scale = span / clen if clen else 1.0
        for b in bands.get(chrom, []):
            base.append(offset + b.start * scale)
            widths.append((b.end - b.start) * scale)
            colors.append(STAIN_COLORS.get(b.stain, "#cccccc"))
            cd.append((chrom.replace("chr", ""), b.name, b.stain))

    fig.add_trace(go.Bar(
        x=widths, base=base, y=[0] * len(base), orientation="h",
        marker=dict(color=colors, line=dict(color=BAND_OUTLINE, width=0.3)),
        customdata=cd, width=0.8,
        hovertemplate="chr%{customdata[0]}%{customdata[1]} (%{customdata[2]})<extra></extra>",
        showlegend=False,
    ), row=2, col=1)


# Injected into the HTML export: builds a checkbox per copy-number (scattergl)
# trace that shows/hides that callability tier via Plotly.restyle. HTML-only (the
# PNG is static); references the figure by its known div id.
_SCATTER_DIV_ID = "chromcov-scatter"
_CHECKBOX_SCRIPT = """
<script>
(function () {
  var gd = document.getElementById('%s');
  if (!gd) return;
  function build() {
    if (!gd.data) { return setTimeout(build, 50); }
    var panel = document.createElement('div');
    panel.style.cssText = 'font-family:sans-serif;font-size:13px;margin:6px 0 4px 72px;';
    var lab = document.createElement('span');
    lab.textContent = 'show strata: ';
    lab.style.marginRight = '6px';
    panel.appendChild(lab);
    var any = false;
    gd.data.forEach(function (tr, i) {
      if (tr.type !== 'scattergl') return;
      any = true;
      var w = document.createElement('label');
      w.style.cssText = 'margin-right:14px;cursor:pointer;';
      var cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.checked = (tr.visible === undefined || tr.visible === true);
      cb.addEventListener('change', function () {
        Plotly.restyle(gd, {visible: cb.checked ? true : 'legendonly'}, [i]);
      });
      w.appendChild(cb);
      w.appendChild(document.createTextNode(' ' + (tr.name || ('trace ' + i))));
      panel.appendChild(w);
    });
    if (any) gd.parentNode.insertBefore(panel, gd);
  }
  build();
})();
</script>
""" % _SCATTER_DIV_ID


def scatter_windows(windows: list[dict], out_path: Path, baseline: float, ploidy: int = 2,
                    cap_cn: float = 6.0, min_easy_frac: float = 0.0) -> dict[str, Path]:
    """windows: rows with chrom/start/mean(/easy_frac/stratum), restricted to
    PRIMARY. Y is copy ratio (ploidy*mean/baseline), capped so pileup spikes don't
    flatten the axis. A GRCh38 cytoband ideogram is drawn beneath, aligned to the
    same genomic x-axis.

    When strata were supplied, every window is shown and colored by its dominant
    callability tier (easy/difficult/extreme), so segmental CN steps and the
    repeat/centromere pileups are visible *and* distinguishable rather than the
    latter being silently dropped. `min_easy_frac > 0` instead restricts to the
    callable ('easy') windows (the old callable-only view).

    `out_path` is the PNG destination; a self-contained interactive HTML is
    written alongside with the same stem. Returns {"png": ..., "html": ...}.
    """
    out_path = Path(out_path)
    fig = build_scatter_figure(windows, baseline, ploidy=ploidy, cap_cn=cap_cn,
                               min_easy_frac=min_easy_frac)
    html_path = out_path.with_suffix(".html")
    fig.write_image(str(out_path), scale=2)
    html = fig.to_html(include_plotlyjs=True, full_html=True, div_id=_SCATTER_DIV_ID)
    html = html.replace("</body>", _CHECKBOX_SCRIPT + "</body>")
    html_path.write_text(html)
    return {"png": out_path, "html": html_path}


def build_scatter_figure(windows: list[dict], baseline: float, ploidy: int = 2,
                         cap_cn: float = 6.0, min_easy_frac: float = 0.0) -> go.Figure:
    """Assemble the copy-number scatter + ideogram figure (see scatter_windows).
    Returned unwritten so callers can embed it (e.g. a Dash dcc.Graph)."""
    layout, rows_by_chrom = _genome_layout(windows, min_easy_frac)

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        row_heights=[0.82, 0.18], vertical_spacing=0.03)

    stratified = _add_scatter(fig, layout, rows_by_chrom, baseline, ploidy, cap_cn)
    _add_ideogram(fig, layout)

    # Chromosome dividers (through both the scatter and the ideogram) + CN=2
    # diploid reference on the scatter row.
    for _chrom, offset, _span, _clen in layout[1:]:
        for r in (1, 2):
            fig.add_vline(x=offset, line=dict(color="#000075", width=2), row=r, col=1)
    fig.add_hline(y=ploidy, line=dict(color="grey", width=1, dash="dash"), row=1, col=1)

    if min_easy_frac:
        note = f" — callable windows only (easy≥{min_easy_frac:g})"
    elif stratified:
        note = " — colored by callability tier"
    else:
        note = ""

    xticks = [offset + span / 2 for _c, offset, span, _cl in layout]
    xlabels = [chrom.replace("chr", "") for chrom, *_ in layout]
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(tickvals=xticks, ticktext=xlabels, tickfont=dict(size=9),
                     showgrid=False, row=2, col=1)
    fig.update_yaxes(title_text=f"approx copy number (cap {cap_cn:g})",
                     range=[0, cap_cn], row=1, col=1)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False,
                     fixedrange=True, row=2, col=1)
    fig.update_layout(
        title=f"Windowed copy number across the genome{note}",
        template="simple_white", barmode="overlay", bargap=0,
        legend=dict(title="callability", itemsizing="constant"),
        width=1600, height=560, margin=dict(l=70, r=30, t=60, b=40),
    )
    return fig
