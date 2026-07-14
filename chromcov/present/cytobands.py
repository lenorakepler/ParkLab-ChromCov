"""
GRCh38 cytogenetic bands, used to draw a chromosome ideogram beneath the
windowed copy-number scatter. The band table is the UCSC hg38 `cytoBand` track
restricted to the primary assembly (chr1-22, X, Y) and bundled in the package so
the plot has no network/reference dependency at render time.

Row format (tab-separated): chrom, start, end, name (band, e.g. p36.33),
Giemsa stain (gneg/gpos*/acen/gvar/stalk).
"""
from __future__ import annotations

from functools import lru_cache
from importlib.resources import files
from typing import NamedTuple

# Giemsa stain -> fill color. gpos* darken with density; acen is the centromere.
STAIN_COLORS = {
    "gneg": "#ffffff",
    "gpos25": "#c8c8c8",
    "gpos50": "#969696",
    "gpos75": "#646464",
    "gpos100": "#3c3c3c",
    "gvar": "#e0e0e0",
    "stalk": "#4f9fd9",
    "acen": "#d62728",
}
BAND_OUTLINE = "#888888"


class Band(NamedTuple):
    start: int
    end: int
    name: str
    stain: str


@lru_cache(maxsize=1)
def load_cytobands() -> dict[str, list[Band]]:
    """chrom -> bands in coordinate order. Cached; parses the bundled TSV once."""
    text = files(__package__).joinpath("data/cytoBand_hg38.tsv").read_text()
    bands: dict[str, list[Band]] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        chrom, start, end, name, stain = line.split("\t")
        bands.setdefault(chrom, []).append(Band(int(start), int(end), name, stain))
    return bands


def chrom_length(chrom: str) -> int | None:
    """Assembled length of `chrom` (last band end), or None if unknown."""
    b = load_cytobands().get(chrom)
    return b[-1].end if b else None
