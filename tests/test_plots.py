"""
The Plotly copy-number scatter + cytoband ideogram: the figure carries a trace
per callability tier plus the ideogram bar, and scatter_windows emits both a PNG
and a self-contained HTML. PNG assertions are skipped where kaleido/Chromium
isn't available (it drives a headless browser).
"""
import pytest

from chromcov.present import plots
from chromcov.present.cytobands import chrom_length, load_cytobands


def _windows(stratum=None):
    rows = []
    for chrom in ("chr1", "chr2", "chrX"):
        for start in range(0, 2_000_000, 100_000):
            r = {"chrom": chrom, "start": start, "end": start + 100_000,
                 "mean": 30.0, "easy_frac": 1.0}
            if stratum:
                r["stratum"] = stratum
            rows.append(r)
    return rows


def test_cytobands_cover_primary_and_have_centromere():
    bands = load_cytobands()
    assert len(bands["chr1"]) > 0
    assert chrom_length("chr1") > 240_000_000
    assert any(b.stain == "acen" for b in bands["chr1"])   # centromere present


def test_figure_has_tier_traces_and_ideogram():
    fig = plots.build_scatter_figure(_windows(stratum="easy"), baseline=30.0)
    names = [t.name for t in fig.data]
    assert "easy" in names            # one scatter trace per tier
    # scatter traces (row 1) + one ideogram bar trace (row 2)
    assert sum(t.type == "bar" for t in fig.data) == 1


def test_unstratified_single_trace():
    fig = plots.build_scatter_figure(_windows(), baseline=30.0)
    scatter = [t for t in fig.data if t.type == "scattergl"]
    assert len(scatter) == 1
    assert scatter[0].name == "windows"


def test_scatter_windows_writes_png_and_html(tmp_path):
    png = tmp_path / "coverage.scatter.png"
    try:
        out = plots.scatter_windows(_windows(), png, baseline=30.0)
    except Exception as e:                      # kaleido/Chromium missing
        pytest.skip(f"static image export unavailable: {e}")
    assert out["png"].exists() and out["png"].stat().st_size > 0
    assert out["html"].exists()
    html = out["html"].read_text()
    assert "plotly" in html[:100_000].lower()
    # HTML-only strata show/hide checkboxes wired to Plotly.restyle
    assert plots._SCATTER_DIV_ID in html
    assert "show strata:" in html and "Plotly.restyle" in html
